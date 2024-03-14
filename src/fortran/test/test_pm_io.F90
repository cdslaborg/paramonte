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

!>  \brief This module contains tests of the module [pm_io](@ref pm_io).
!>  \author Amir Shahmoradi

module test_pm_io

    use pm_io
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none

        test = test_type(MODULE_NAME)
        call test%run(test_getRecl_1, SK_"test_getRecl_1")
        call test%run(test_getRecl_2, SK_"test_getRecl_2")
        call test%run(test_getRecl_3, SK_"test_getRecl_3")
        call test%run(test_getRecl_4, SK_"test_getRecl_4")
        call test%run(test_getForm_1, SK_"test_getForm_1")
        call test%run(test_getForm_2, SK_"test_getForm_2")
        call test%run(test_getForm_3, SK_"test_getForm_3")
        call test%run(test_getForm_4, SK_"test_getForm_4")
        call test%run(test_getName_1, SK_"test_getName_1")
        call test%run(test_getName_2, SK_"test_getName_2")
        call test%run(test_getName_3, SK_"test_getName_3")
        call test%run(test_getName_4, SK_"test_getName_4")
        call test%run(test_getBlank_1, SK_"test_getBlank_1")
        call test%run(test_getBlank_2, SK_"test_getBlank_2")
        call test%run(test_getBlank_3, SK_"test_getBlank_3")
        call test%run(test_getBlank_4, SK_"test_getBlank_4")
        call test%run(test_getDelim_1, SK_"test_getDelim_1")
        call test%run(test_getDelim_2, SK_"test_getDelim_2")
        call test%run(test_getDelim_3, SK_"test_getDelim_3")
        call test%run(test_getDelim_4, SK_"test_getDelim_4")
        call test%run(test_getAction_1, SK_"test_getAction_1")
        call test%run(test_getAction_2, SK_"test_getAction_2")
        call test%run(test_getAction_3, SK_"test_getAction_3")
        call test%run(test_getAction_4, SK_"test_getAction_4")
        call test%run(test_getAccess_1, SK_"test_getAccess_1")
        call test%run(test_getAccess_2, SK_"test_getAccess_2")
        call test%run(test_getAccess_3, SK_"test_getAccess_3")
        call test%run(test_getAccess_4, SK_"test_getAccess_4")
        call test%run(test_getNumber_1, SK_"test_getNumber_1")
        call test%run(test_getNumber_2, SK_"test_getNumber_2")
        call test%run(test_getNumber_3, SK_"test_getNumber_3")
        call test%run(test_getNumber_4, SK_"test_getNumber_4")
        call test%run(test_getPosition_1, SK_"test_getPosition_1")
        call test%run(test_getPosition_2, SK_"test_getPosition_2")
        call test%run(test_getPosition_3, SK_"test_getPosition_3")
        call test%run(test_getPosition_4, SK_"test_getPosition_4")
        call test%run(test_constructPad_1, SK_"test_constructPad_1")
        call test%run(test_constructPad_2, SK_"test_constructPad_2")
        call test%run(test_constructPad_3, SK_"test_constructPad_3")
        call test%run(test_constructPad_4, SK_"test_constructPad_4")
        call test%run(test_constructPad_5, SK_"test_constructPad_5")
        call test%run(test_constructFile_1, SK_"test_constructFile_1")
        call test%run(test_constructFile_2, SK_"test_constructFile_2")
        call test%run(test_constructFile_3, SK_"test_constructFile_3")
        call test%run(test_constructForm_1, SK_"test_constructForm_1")
        call test%run(test_constructForm_2, SK_"test_constructForm_2")
        call test%run(test_constructForm_3, SK_"test_constructForm_3")
        call test%run(test_constructForm_4, SK_"test_constructForm_4")
        call test%run(test_constructForm_5, SK_"test_constructForm_5")
        call test%run(test_constructSign_1, SK_"test_constructSign_1")
        call test%run(test_constructSign_2, SK_"test_constructSign_2")
        call test%run(test_constructSign_3, SK_"test_constructSign_3")
        call test%run(test_constructSign_4, SK_"test_constructSign_4")
        call test%run(test_constructSign_5, SK_"test_constructSign_5")
        call test%run(test_constructSign_6, SK_"test_constructSign_6")
        call test%run(test_getOpenStatus_1, SK_"test_getOpenStatus_1")
        call test%run(test_getOpenStatus_2, SK_"test_getOpenStatus_2")
        call test%run(test_getOpenStatus_3, SK_"test_getOpenStatus_3")
        call test%run(test_getOpenStatus_4, SK_"test_getOpenStatus_4")
        call test%run(test_isExtant_1, SK_"test_isExtant_1")
        call test%run(test_isExtant_2, SK_"test_isExtant_2")
        call test%run(test_isExtant_3, SK_"test_isExtant_3")
        call test%run(test_constructBlank_1, SK_"test_constructBlank_1")
        call test%run(test_constructBlank_2, SK_"test_constructBlank_2")
        call test%run(test_constructBlank_3, SK_"test_constructBlank_3")
        call test%run(test_constructBlank_4, SK_"test_constructBlank_4")
        call test%run(test_constructBlank_5, SK_"test_constructBlank_5")
        call test%run(test_constructDelim_1, SK_"test_constructDelim_1")
        call test%run(test_constructDelim_2, SK_"test_constructDelim_2")
        call test%run(test_constructDelim_3, SK_"test_constructDelim_3")
        call test%run(test_constructDelim_4, SK_"test_constructDelim_4")
        call test%run(test_constructDelim_5, SK_"test_constructDelim_5")
        call test%run(test_constructDelim_6, SK_"test_constructDelim_6")
        call test%run(test_constructRound_1, SK_"test_constructRound_1")
        call test%run(test_constructRound_2, SK_"test_constructRound_2")
        call test%run(test_constructRound_3, SK_"test_constructRound_3")
        call test%run(test_constructRound_4, SK_"test_constructRound_4")
        call test%run(test_constructRound_5, SK_"test_constructRound_5")
        call test%run(test_constructRound_6, SK_"test_constructRound_6")
        call test%run(test_constructRound_7, SK_"test_constructRound_7")
        call test%run(test_constructRound_8, SK_"test_constructRound_8")
        call test%run(test_constructRound_9, SK_"test_constructRound_9")
        call test%run(test_constructAction_1, SK_"test_constructAction_1")
        call test%run(test_constructAction_2, SK_"test_constructAction_2")
        call test%run(test_constructAction_3, SK_"test_constructAction_3")
        call test%run(test_constructAction_4, SK_"test_constructAction_4")
        call test%run(test_constructAction_5, SK_"test_constructAction_5")
        call test%run(test_constructAccess_1, SK_"test_constructAccess_1")
        call test%run(test_constructAccess_2, SK_"test_constructAccess_2")
        call test%run(test_constructAccess_3, SK_"test_constructAccess_3")
        call test%run(test_constructAccess_4, SK_"test_constructAccess_4")
        call test%run(test_constructAccess_5, SK_"test_constructAccess_5")
        call test%run(test_constructPosition_1, SK_"test_constructPosition_1")
        call test%run(test_constructPosition_2, SK_"test_constructPosition_2")
        call test%run(test_constructPosition_3, SK_"test_constructPosition_3")
        call test%run(test_constructPosition_4, SK_"test_constructPosition_4")
        call test%run(test_constructPosition_5, SK_"test_constructPosition_5")
        call test%run(test_constructPosition_6, SK_"test_constructPosition_6")
        call test%run(test_isCloseErr_1, SK_"test_isCloseErr_1")
        call test%run(test_isCloseErr_2, SK_"test_isCloseErr_2")
        call test%run(test_isCloseErr_3, SK_"test_isCloseErr_3")
        call test%run(test_getWriteErr_1, SK_"test_getWriteErr_1")
        call test%run(test_getWriteErr_2, SK_"test_getWriteErr_2")
        call test%run(test_getWriteErr_3, SK_"test_getWriteErr_3")
        call test%run(test_isOpenErr_1, SK_"test_isOpenErr_1")
        call test%run(test_isOpenErr_2, SK_"test_isOpenErr_2")
        call test%run(test_isOpenErr_3, SK_"test_isOpenErr_3")
        call test%run(test_getReadErr_1, SK_"test_getReadErr_1")
        call test%run(test_getReadErr_2, SK_"test_getReadErr_2")
        call test%run(test_getReadErr_3, SK_"test_getReadErr_3")
        call test%run(test_getReadErr_4, SK_"test_getReadErr_4")
        call test%run(test_getReadErr_5, SK_"test_getReadErr_5")
        call test%run(test_getReadErr_6, SK_"test_getReadErr_6")
        call test%run(test_isInqErr_1, SK_"test_isInqErr_1")
        call test%run(test_isInqErr_2, SK_"test_isInqErr_2")
        call test%run(test_isInqErr_3, SK_"test_isInqErr_3")
        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructFile_1() result(assertion)

        implicit none
        logical(LK) :: assertion
        type(file_type) :: file

        file = file_type()
        assertion = .not. file%err%occurred
        call test%assert(assertion, desc = "No error should occur upon `file` construction.")

        assertion = assertion .and. file%path%original  == ""
        call test%assert(assertion)
        assertion = assertion .and. file%path%modified  == ""
        call test%assert(assertion)
        assertion = assertion .and. file%path%dir       == ""
        call test%assert(assertion)
        assertion = assertion .and. file%path%name      == ""
        call test%assert(assertion)
        assertion = assertion .and. file%path%ext       == ""
        call test%assert(assertion)
        assertion = assertion .and. file%unit           == -2147483647
        call test%assert(assertion)
        assertion = assertion .and. file%number         == -2147483647
        call test%assert(assertion)
        assertion = assertion .and. file%recl           == -2147483647
        call test%assert(assertion)
        assertion = assertion .and. file%exists         .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. file%isOpen         .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. file%isNamed        .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. file%isNumbered     .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. file%status         == "unknown"
        call test%assert(assertion)
        assertion = assertion .and. file%asynchronous   == "no"
        call test%assert(assertion)
        assertion = assertion .and. file%format         == ""
        call test%assert(assertion)
        assertion = assertion .and. file%nameByCompiler == ""
        call test%assert(assertion)
        assertion = assertion .and. file%action%val   == "readwrite"
        call test%assert(assertion)
        assertion = assertion .and. file%access%val   == "sequential"
        call test%assert(assertion)
        assertion = assertion .and. file%form%val     == "formatted"
        call test%assert(assertion)
        assertion = assertion .and. file%blank%val    == "null"
        call test%assert(assertion)
        assertion = assertion .and. file%position%val == "asis"
        call test%assert(assertion)
        assertion = assertion .and. file%delim%val    == "none"
        call test%assert(assertion)
        assertion = assertion .and. file%pad%val      == "yes"
        call test%assert(assertion)
        assertion = assertion .and. file%round%val    == "processor_defined"
        call test%assert(assertion)
        assertion = assertion .and. file%sign%val     == "processor_defined"
        call test%assert(assertion)
        assertion = assertion .and. file%err%occurred   .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. file%err%msg        == ""
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))")   "file%path%original   : ", file%path%original
            write(test%disp%unit,"(*(g0))")   "file%path%modified   : ", file%path%modified
            write(test%disp%unit,"(*(g0))")   "file%path%dir        : ", file%path%dir
            write(test%disp%unit,"(*(g0))")   "file%path%name       : ", file%path%name
            write(test%disp%unit,"(*(g0))")   "file%path%ext        : ", file%path%ext
            write(test%disp%unit,"(*(g0))")   "file%path%sep        : ", file%path%sep
            write(test%disp%unit,"(*(g0))")   "file%unit            : ", file%unit
            write(test%disp%unit,"(*(g0))")   "file%number          : ", file%number
            write(test%disp%unit,"(*(g0))")   "file%recl            : ", file%recl
            write(test%disp%unit,"(*(g0))")   "file%exists          : ", file%exists
            write(test%disp%unit,"(*(g0))")   "file%isOpen          : ", file%isOpen
            write(test%disp%unit,"(*(g0))")   "file%isNamed         : ", file%isNamed
            write(test%disp%unit,"(*(g0))")   "file%isNumbered      : ", file%isNumbered
            write(test%disp%unit,"(*(g0))")   "file%status          : ", file%status
            write(test%disp%unit,"(*(g0))")   "file%asynchronous    : ", file%asynchronous
            write(test%disp%unit,"(*(g0))")   "file%format          : ", file%format
            write(test%disp%unit,"(*(g0))")   "file%nameByCompiler  : ", file%nameByCompiler
            write(test%disp%unit,"(*(g0))")   "file%action%val      : ", file%action%val
            write(test%disp%unit,"(*(g0))")   "file%access%val      : ", file%access%val
            write(test%disp%unit,"(*(g0))")   "file%form%val        : ", file%form%val
            write(test%disp%unit,"(*(g0))")   "file%blank%val       : ", file%blank%val
            write(test%disp%unit,"(*(g0))")   "file%position%val    : ", file%position%val
            write(test%disp%unit,"(*(g0))")   "file%delim%val       : ", file%delim%val
            write(test%disp%unit,"(*(g0))")   "file%pad%val         : ", file%pad%val
            write(test%disp%unit,"(*(g0))")   "file%round%val       : ", file%round%val
            write(test%disp%unit,"(*(g0))")   "file%sign%val        : ", file%sign%val
            write(test%disp%unit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructFile_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructFile_2() result(assertion)

        implicit none
        logical(LK) :: assertion
        type(file_type) :: file

        file = file_type( unit=13 &
                        , recl=9999 &
                        , path="./test_File_pmod/\test_File_pmod\-" &
                        , status="NEw" &
                        , position="APPEND" &
                        , access="direCt" &
                        , form="unformatted" &
                        , action="rEAd" &
                        , delim="quotE" &
                        , round="uP" &
                        , sign="undefined" &
                        , pad = "nO" &
                        , blank = "undefined" &
                        , format = "(A)" &
                        , asynchronous = "Yes" &
                        )

        assertion = .not. file%err%occurred
        call test%assert(assertion)

        assertion = assertion .and. file%path%original  == "./test_File_pmod/\test_File_pmod\-"

        if (file%path%sep=="/") then
            assertion = assertion .and. file%path%modified  == "./test_File_pmod/\test_File_pmod\-"
            call test%assert(assertion)
            assertion = assertion .and. file%path%dir       == "./test_File_pmod/"
            call test%assert(assertion)
            assertion = assertion .and. file%path%name      == "\test_File_pmod\-"
            call test%assert(assertion)
            assertion = assertion .and. file%path%ext       == ""
            call test%assert(assertion)
#if WINDOWS_ENABLED
        else
            assertion = assertion .and. file%path%modified  == ".\test_File_pmod\\test_File_pmod\-"
            call test%assert(assertion)
            assertion = assertion .and. file%path%dir       == ".\test_File_pmod\\test_File_pmod\"
            call test%assert(assertion)
            assertion = assertion .and. file%path%name      == "-"
            call test%assert(assertion)
            assertion = assertion .and. file%path%ext       == ""
            call test%assert(assertion)
#endif
        end if
        assertion = assertion .and. file%unit           == 13
        call test%assert(assertion)
        assertion = assertion .and. file%number         == -2147483647
        call test%assert(assertion)
        assertion = assertion .and. file%recl           == 9999
        call test%assert(assertion)
        assertion = assertion .and. file%exists         .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. file%isOpen         .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. file%isNamed        .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. file%isNumbered     .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. file%status         == "new"
        call test%assert(assertion)
        assertion = assertion .and. file%asynchronous   == "yes"
        call test%assert(assertion)
        assertion = assertion .and. file%format         == "(A)"
        call test%assert(assertion)
        assertion = assertion .and. file%nameByCompiler == ""
        call test%assert(assertion)
        assertion = assertion .and. file%action%val   == "read"
        call test%assert(assertion)
        assertion = assertion .and. file%access%val   == "direct"
        call test%assert(assertion)
        assertion = assertion .and. file%form%val     == "unformatted"
        call test%assert(assertion)
        assertion = assertion .and. file%blank%val    == "undefined"
        call test%assert(assertion)
        assertion = assertion .and. file%position%val == "append"
        call test%assert(assertion)
        assertion = assertion .and. file%delim%val    == "quote"
        call test%assert(assertion)
        assertion = assertion .and. file%pad%val      == "no"
        call test%assert(assertion)
        assertion = assertion .and. file%round%val    == "up"
        call test%assert(assertion)
        assertion = assertion .and. file%sign%val     == "undefined"
        call test%assert(assertion)
        assertion = assertion .and. file%err%occurred   .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. file%err%msg        == ""
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))")   "file%path%original  : ", file%path%original
            write(test%disp%unit,"(*(g0))")   "file%path%modified  : ", file%path%modified
            write(test%disp%unit,"(*(g0))")   "file%path%dir       : ", file%path%dir
            write(test%disp%unit,"(*(g0))")   "file%path%name      : ", file%path%name
            write(test%disp%unit,"(*(g0))")   "file%path%ext       : ", file%path%ext
            write(test%disp%unit,"(*(g0))")   "file%path%sep: ", file%path%sep
            write(test%disp%unit,"(*(g0))")   "file%unit           : ", file%unit
            write(test%disp%unit,"(*(g0))")   "file%number         : ", file%number
            write(test%disp%unit,"(*(g0))")   "file%recl           : ", file%recl
            write(test%disp%unit,"(*(g0))")   "file%exists         : ", file%exists
            write(test%disp%unit,"(*(g0))")   "file%isOpen         : ", file%isOpen
            write(test%disp%unit,"(*(g0))")   "file%isNamed        : ", file%isNamed
            write(test%disp%unit,"(*(g0))")   "file%isNumbered     : ", file%isNumbered
            write(test%disp%unit,"(*(g0))")   "file%status         : ", file%status
            write(test%disp%unit,"(*(g0))")   "file%asynchronous   : ", file%asynchronous
            write(test%disp%unit,"(*(g0))")   "file%format         : ", file%format
            write(test%disp%unit,"(*(g0))")   "file%nameByCompiler : ", file%nameByCompiler
            write(test%disp%unit,"(*(g0))")   "file%action%val   : ", file%action%val
            write(test%disp%unit,"(*(g0))")   "file%access%val   : ", file%access%val
            write(test%disp%unit,"(*(g0))")   "file%form%val     : ", file%form%val
            write(test%disp%unit,"(*(g0))")   "file%blank%val    : ", file%blank%val
            write(test%disp%unit,"(*(g0))")   "file%position%val : ", file%position%val
            write(test%disp%unit,"(*(g0))")   "file%delim%val    : ", file%delim%val
            write(test%disp%unit,"(*(g0))")   "file%pad%val      : ", file%pad%val
            write(test%disp%unit,"(*(g0))")   "file%round%val    : ", file%round%val
            write(test%disp%unit,"(*(g0))")   "file%sign%val     : ", file%sign%val
            write(test%disp%unit,"(*(g0))")   "file%err%occurred   : ", file%err%occurred
            write(test%disp%unit,"(*(g0))")   "file%err%msg        : ", file%err%msg
            write(test%disp%unit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructFile_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Test the effects of missing input argument `form`.
    function test_constructFile_3() result(assertion)

        implicit none
        logical(LK) :: assertion
        type(file_type) :: file

        file = file_type( unit=13 &
                        , recl=9999 &
                        , path="./test_File_pmod/\test_File_pmod\-" &
                        , status="NEw" &
                        , position="APPEND" &
                        , access="direCt" &
                        , action="rEAd" &
                        , delim="quotE" &
                        , round="uP" &
                        , sign="undefined" &
                        , pad = "nO" &
                        , blank = "undefined" &
                        , format = "(A)" &
                        , asynchronous = "Yes" &
                        )

        assertion = .not. file%err%occurred
        call test%assert(assertion)

        assertion = assertion .and. file%path%original  == "./test_File_pmod/\test_File_pmod\-"
        call test%assert(assertion)

        if (file%path%sep=="/") then
            assertion = assertion .and. file%path%modified  == "./test_File_pmod/\test_File_pmod\-"
            call test%assert(assertion)
            assertion = assertion .and. file%path%dir       == "./test_File_pmod/"
            call test%assert(assertion)
            assertion = assertion .and. file%path%name      == "\test_File_pmod\-"
            call test%assert(assertion)
            assertion = assertion .and. file%path%ext       == ""
            call test%assert(assertion)
#if WINDOWS_ENABLED
        else
            assertion = assertion .and. file%path%modified  == ".\test_File_pmod\\test_File_pmod\-"
            call test%assert(assertion)
            assertion = assertion .and. file%path%dir       == ".\test_File_pmod\\test_File_pmod\"
            call test%assert(assertion)
            assertion = assertion .and. file%path%name      == "-"
            call test%assert(assertion)
            assertion = assertion .and. file%path%ext       == ""
            call test%assert(assertion)
#endif
        end if

        assertion = assertion .and. file%unit           == 13
        call test%assert(assertion)
        assertion = assertion .and. file%number         == -2147483647
        call test%assert(assertion)
        assertion = assertion .and. file%recl           == 9999
        call test%assert(assertion)
        assertion = assertion .and. file%exists         .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. file%isOpen         .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. file%isNamed        .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. file%isNumbered     .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. file%status         == "new"
        call test%assert(assertion)
        assertion = assertion .and. file%asynchronous   == "yes"
        call test%assert(assertion)
        assertion = assertion .and. file%format         == "(A)"
        call test%assert(assertion)
        assertion = assertion .and. file%nameByCompiler == ""
        call test%assert(assertion)
        assertion = assertion .and. file%action%val   == "read"
        call test%assert(assertion)
        assertion = assertion .and. file%access%val   == "direct"
        call test%assert(assertion)
        assertion = assertion .and. file%form%val     == "unformatted"
        call test%assert(assertion)
        assertion = assertion .and. file%blank%val    == "undefined"
        call test%assert(assertion)
        assertion = assertion .and. file%position%val == "append"
        call test%assert(assertion)
        assertion = assertion .and. file%delim%val    == "quote"
        call test%assert(assertion)
        assertion = assertion .and. file%pad%val      == "no"
        call test%assert(assertion)
        assertion = assertion .and. file%round%val    == "up"
        call test%assert(assertion)
        assertion = assertion .and. file%sign%val     == "undefined"
        call test%assert(assertion)
        assertion = assertion .and. file%err%occurred   .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. file%err%msg        == ""
        call test%assert(assertion)

        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))")   "file%path%original  : ", file%path%original
            write(test%disp%unit,"(*(g0))")   "file%path%modified  : ", file%path%modified
            write(test%disp%unit,"(*(g0))")   "file%path%dir       : ", file%path%dir
            write(test%disp%unit,"(*(g0))")   "file%path%name      : ", file%path%name
            write(test%disp%unit,"(*(g0))")   "file%path%ext       : ", file%path%ext
            write(test%disp%unit,"(*(g0))")   "file%path%sep: ", file%path%sep
            write(test%disp%unit,"(*(g0))")   "file%unit           : ", file%unit
            write(test%disp%unit,"(*(g0))")   "file%number         : ", file%number
            write(test%disp%unit,"(*(g0))")   "file%recl           : ", file%recl
            write(test%disp%unit,"(*(g0))")   "file%exists         : ", file%exists
            write(test%disp%unit,"(*(g0))")   "file%isOpen         : ", file%isOpen
            write(test%disp%unit,"(*(g0))")   "file%isNamed        : ", file%isNamed
            write(test%disp%unit,"(*(g0))")   "file%isNumbered     : ", file%isNumbered
            write(test%disp%unit,"(*(g0))")   "file%status         : ", file%status
            write(test%disp%unit,"(*(g0))")   "file%asynchronous   : ", file%asynchronous
            write(test%disp%unit,"(*(g0))")   "file%format         : ", file%format
            write(test%disp%unit,"(*(g0))")   "file%nameByCompiler : ", file%nameByCompiler
            write(test%disp%unit,"(*(g0))")   "file%action%val   : ", file%action%val
            write(test%disp%unit,"(*(g0))")   "file%access%val   : ", file%access%val
            write(test%disp%unit,"(*(g0))")   "file%form%val     : ", file%form%val
            write(test%disp%unit,"(*(g0))")   "file%blank%val    : ", file%blank%val
            write(test%disp%unit,"(*(g0))")   "file%position%val : ", file%position%val
            write(test%disp%unit,"(*(g0))")   "file%delim%val    : ", file%delim%val
            write(test%disp%unit,"(*(g0))")   "file%pad%val      : ", file%pad%val
            write(test%disp%unit,"(*(g0))")   "file%round%val    : ", file%round%val
            write(test%disp%unit,"(*(g0))")   "file%sign%val     : ", file%sign%val
            write(test%disp%unit,"(*(g0))")   "file%err%occurred   : ", file%err%occurred
            write(test%disp%unit,"(*(g0))")   "file%err%msg        : ", file%err%msg
            write(test%disp%unit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructFile_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructPad_1() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Pad_type) :: pad
        pad = Pad_type()
        assertion = .not. pad%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. pad%val == "yes"
        call test%assert(assertion)
        assertion = assertion .and. pad%isPadded .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. pad%isNotPadded .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. pad%isUndefined .eqv. .false._LK
    end function test_constructPad_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Testing `Pad_type( value = 'undefined' )`
    function test_constructPad_2() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Pad_type) :: pad
        pad = Pad_type(value='undefined')
        assertion = .not. pad%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. pad%val == "undefined"
        call test%assert(assertion)
        assertion = assertion .and. pad%isPadded .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. pad%isNotPadded .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. pad%isUndefined .eqv. .true._LK
    end function test_constructPad_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Testing `Pad_type( value = 'yes' )`
    function test_constructPad_3() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Pad_type) :: pad
        pad = Pad_type(value='yes')
        assertion = .not. pad%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. pad%val == "yes"
        call test%assert(assertion)
        assertion = assertion .and. pad%isPadded .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. pad%isNotPadded .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. pad%isUndefined .eqv. .false._LK
    end function test_constructPad_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Testing `Pad_type( value = 'no' )`
    function test_constructPad_4() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Pad_type) :: pad
        pad = Pad_type(value='no')
        assertion = .not. pad%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. pad%val == "no"
        call test%assert(assertion)
        assertion = assertion .and. pad%isPadded .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. pad%isNotPadded .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. pad%isUndefined .eqv. .false._LK
    end function test_constructPad_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Testing `Pad_type( value = 'nonsense' )`
    function test_constructPad_5() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Pad_type) :: pad
        pad = Pad_type(value="nonsense")
        assertion = pad%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. pad%val == ""
        call test%assert(assertion)
        assertion = assertion .and. pad%isPadded .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. pad%isNotPadded .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. pad%isUndefined .eqv. .false._LK
    end function test_constructPad_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAction_1() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Action_type) :: action
        action = Action_type()
        assertion = .not. action%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. action%val == "readwrite"
        call test%assert(assertion)
        assertion = assertion .and. action%isRead .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. action%isWrite .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. action%isReadWrite .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. action%isUndefined .eqv. .false._LK
    end function test_constructAction_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAction_2() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Action_type) :: action
        action = Action_type(value="undefined")
        assertion = .not. action%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. action%val == "undefined"
        call test%assert(assertion)
        assertion = assertion .and. action%isRead .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. action%isWrite .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. action%isReadWrite .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. action%isUndefined .eqv. .true._LK
    end function test_constructAction_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAction_3() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Action_type) :: action
        action = Action_type(value="read")
        assertion = .not. action%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. action%val == "read"
        call test%assert(assertion)
        assertion = assertion .and. action%isRead .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. action%isWrite .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. action%isReadWrite .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. action%isUndefined .eqv. .false._LK
    end function test_constructAction_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAction_4() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Action_type) :: action
        action = Action_type(value="write")
        assertion = .not. action%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. action%val == "write"
        call test%assert(assertion)
        assertion = assertion .and. action%isRead .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. action%isWrite .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. action%isReadWrite .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. action%isUndefined .eqv. .false._LK
    end function test_constructAction_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAction_5() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Action_type) :: action
        action = Action_type(value="nonsense")
        assertion = action%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. action%val == ""
        call test%assert(assertion)
        assertion = assertion .and. action%isRead .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. action%isWrite .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. action%isReadWrite .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. action%isUndefined .eqv. .false._LK
    end function test_constructAction_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAccess_1() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Access_type) :: access
        access = Access_type()
        assertion = .not. access%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. access%val == "sequential"
        call test%assert(assertion)
        assertion = assertion .and. access%isSequential .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. access%isDirect .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. access%isUndefined .eqv. .false._LK
    end function test_constructAccess_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAccess_2() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Access_type) :: access
        access = Access_type("undefined")
        assertion = .not. access%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. access%val == "undefined"
        call test%assert(assertion)
        assertion = assertion .and. access%isSequential .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. access%isDirect .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. access%isUndefined .eqv. .true._LK
    end function test_constructAccess_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAccess_3() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Access_type) :: access
        access = Access_type("sequential")
        assertion = .not. access%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. access%val == "sequential"
        call test%assert(assertion)
        assertion = assertion .and. access%isSequential .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. access%isDirect .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. access%isUndefined .eqv. .false._LK
    end function test_constructAccess_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAccess_4() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Access_type) :: access
        access = Access_type("direct")
        assertion = .not. access%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. access%val == "direct"
        call test%assert(assertion)
        assertion = assertion .and. access%isSequential .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. access%isDirect .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. access%isUndefined .eqv. .false._LK
    end function test_constructAccess_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAccess_5() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Access_type) :: access
        access = Access_type("nonsense")
        assertion = access%err%occurred
       !call test%assert(assertion)
       !assertion = assertion .and. access%val == "nonsense"
       !call test%assert(assertion)
       !assertion = assertion .and. access%isSequential .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. access%isDirect .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. access%isUndefined .eqv. .false._LK
    end function test_constructAccess_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructForm_1() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Form_type) :: form
        form = Form_type()
        assertion = .not. form%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. form%val == "formatted"
        call test%assert(assertion)
        assertion = assertion .and. form%isFormatted .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. form%isUnformatted .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. form%isUndefined .eqv. .false._LK
    end function test_constructForm_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructForm_2() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Form_type) :: form
        form = Form_type("undefined")
        assertion = .not. form%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. form%val == "undefined"
        call test%assert(assertion)
        assertion = assertion .and. form%isFormatted .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. form%isUnformatted .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. form%isUndefined .eqv. .true._LK
    end function test_constructForm_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructForm_3() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Form_type) :: form
        form = Form_type("formatted")
        assertion = .not. form%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. form%val == "formatted"
        call test%assert(assertion)
        call test%assert(assertion)
        assertion = assertion .and. form%isFormatted .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. form%isUnformatted .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. form%isUndefined .eqv. .false._LK
    end function test_constructForm_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructForm_4() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Form_type) :: form
        form = Form_type("unformatted")
        assertion = .not. form%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. form%val == "unformatted"
        call test%assert(assertion)
        assertion = assertion .and. form%isFormatted .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. form%isUnformatted .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. form%isUndefined .eqv. .false._LK
    end function test_constructForm_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructForm_5() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Form_type) :: form
        form = Form_type("nonsense")
        assertion = form%err%occurred
       !call test%assert(assertion)
       !assertion = assertion .and. form%val == "nonsense"
       !call test%assert(assertion)
       !assertion = assertion .and. form%isFormatted .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. form%isUnformatted .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. form%isUndefined .eqv. .false._LK
    end function test_constructForm_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructBlank_1() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Blank_type) :: blank
        blank = Blank_type()
        assertion = .not. blank%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. blank%val == "null"
        call test%assert(assertion)
        assertion = assertion .and. blank%isNull .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. blank%isZero .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. blank%isUndefined .eqv. .false._LK
    end function test_constructBlank_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructBlank_2() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Blank_type) :: blank
        blank = Blank_type("undefined")
        assertion = .not. blank%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. blank%val == "undefined"
        call test%assert(assertion)
        assertion = assertion .and. blank%isNull .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. blank%isZero .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. blank%isUndefined .eqv. .true._LK
    end function test_constructBlank_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructBlank_3() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Blank_type) :: blank
        blank = Blank_type("null")
        assertion = .not. blank%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. blank%val == "null"
        call test%assert(assertion)
        assertion = assertion .and. blank%isNull .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. blank%isZero .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. blank%isUndefined .eqv. .false._LK
    end function test_constructBlank_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructBlank_4() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Blank_type) :: blank
        blank = Blank_type("zero")
        assertion = .not. blank%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. blank%val == "zero"
        call test%assert(assertion)
        assertion = assertion .and. blank%isNull .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. blank%isZero .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. blank%isUndefined .eqv. .false._LK
    end function test_constructBlank_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructBlank_5() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Blank_type) :: blank
        blank = Blank_type("nonsense")
        assertion = blank%err%occurred
       !call test%assert(assertion)
       !assertion = assertion .and. blank%val == "nonsense"
       !call test%assert(assertion)
       !assertion = assertion .and. blank%isNull .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. blank%isZero .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. blank%isUndefined .eqv. .false._LK
    end function test_constructBlank_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructPosition_1() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Position_type) :: position
        position = Position_type()
        assertion = .not. position%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. position%val == "asis"
        call test%assert(assertion)
        assertion = assertion .and. position%isRewind .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. position%isAppend .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. position%isAsis .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. position%isUndefined .eqv. .false._LK
    end function test_constructPosition_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructPosition_2() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Position_type) :: position
        position = Position_type("undefined")
        assertion = .not. position%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. position%val == "undefined"
        call test%assert(assertion)
        assertion = assertion .and. position%isRewind .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. position%isAppend .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. position%isAsis .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. position%isUndefined .eqv. .true._LK
    end function test_constructPosition_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructPosition_3() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Position_type) :: position
        position = Position_type("rewind")
        assertion = .not. position%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. position%val == "rewind"
        call test%assert(assertion)
        assertion = assertion .and. position%isRewind .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. position%isAppend .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. position%isAsis .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. position%isUndefined .eqv. .false._LK
    end function test_constructPosition_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructPosition_4() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Position_type) :: position
        position = Position_type("APPEND")
        assertion = .not. position%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. position%val == "append"
        call test%assert(assertion)
        assertion = assertion .and. position%isRewind .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. position%isAppend .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. position%isAsis .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. position%isUndefined .eqv. .false._LK
    end function test_constructPosition_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructPosition_5() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Position_type) :: position
        position = Position_type("ASIS")
        assertion = .not. position%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. position%val == "asis"
        call test%assert(assertion)
        assertion = assertion .and. position%isRewind .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. position%isAppend .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. position%isAsis .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. position%isUndefined .eqv. .false._LK
    end function test_constructPosition_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructPosition_6() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Position_type) :: position
        position = Position_type("nonsense")
        assertion = position%err%occurred
       !call test%assert(assertion)
       !assertion = assertion .and. position%val == ""
       !call test%assert(assertion)
       !assertion = assertion .and. position%isRewind .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. position%isAppend .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. position%isAsis .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. position%isUndefined .eqv. .false._LK
    end function test_constructPosition_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructDelim_1() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Delim_type) :: delim
        delim = Delim_type()
        assertion = .not. delim%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. delim%val == "none"
        call test%assert(assertion)
        assertion = assertion .and. delim%isQuote .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. delim%isApostrophe .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. delim%isNone .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. delim%isUndefined .eqv. .false._LK
    end function test_constructDelim_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructDelim_2() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Delim_type) :: delim
        delim = Delim_type("Undefined")
        assertion = .not. delim%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. delim%val == "undefined"
        call test%assert(assertion)
        assertion = assertion .and. delim%isQuote .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. delim%isApostrophe .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. delim%isNone .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. delim%isUndefined .eqv. .true._LK
    end function test_constructDelim_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructDelim_3() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Delim_type) :: delim
        delim = Delim_type("Quote")
        assertion = .not. delim%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. delim%val == "quote"
        call test%assert(assertion)
        assertion = assertion .and. delim%isQuote .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. delim%isApostrophe .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. delim%isNone .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. delim%isUndefined .eqv. .false._LK
    end function test_constructDelim_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructDelim_4() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Delim_type) :: delim
        delim = Delim_type("Apostrophe")
        assertion = .not. delim%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. delim%val == "apostrophe"
        call test%assert(assertion)
        assertion = assertion .and. delim%isQuote .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. delim%isApostrophe .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. delim%isNone .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. delim%isUndefined .eqv. .false._LK
    end function test_constructDelim_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructDelim_5() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Delim_type) :: delim
        delim = Delim_type("None")
        assertion = .not. delim%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. delim%val == "none"
        call test%assert(assertion)
        assertion = assertion .and. delim%isQuote .eqv. .false._LK
        call test%assert(assertion)
        call test%assert(assertion)
        assertion = assertion .and. delim%isApostrophe .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. delim%isNone .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. delim%isUndefined .eqv. .false._LK
    end function test_constructDelim_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructDelim_6() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Delim_type) :: delim
        delim = Delim_type("nonsense")
        assertion = delim%err%occurred
       !call test%assert(assertion)
       !assertion = assertion .and. delim%val == ""
       !call test%assert(assertion)
       !assertion = assertion .and. delim%isQuote .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. delim%isApostrophe .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. delim%isNone .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. delim%isUndefined .eqv. .false._LK
    end function test_constructDelim_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructRound_1() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Round_type) :: round
        round = Round_type()
        assertion = .not. round%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. round%val == "processor_defined"
        call test%assert(assertion)
        assertion = assertion .and. round%isUp .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isDown .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isZero .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isNearest .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isCompatible .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isProcessDefined .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isUndefined .eqv. .false._LK
    end function test_constructRound_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructRound_2() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Round_type) :: round
        round = Round_type("UP")
        assertion = .not. round%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. round%val == "up"
        call test%assert(assertion)
        assertion = assertion .and. round%isUp .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isDown .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isZero .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isNearest .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isCompatible .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isProcessDefined .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isUndefined .eqv. .false._LK
    end function test_constructRound_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructRound_3() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Round_type) :: round
        round = Round_type("Down")
        assertion = .not. round%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. round%val == "down"
        call test%assert(assertion)
        assertion = assertion .and. round%isUp .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isDown .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isZero .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isNearest .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isCompatible .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isProcessDefined .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isUndefined .eqv. .false._LK
    end function test_constructRound_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructRound_4() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Round_type) :: round
        round = Round_type("Zero")
        assertion = .not. round%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. round%val == "zero"
        call test%assert(assertion)
        assertion = assertion .and. round%isUp .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isDown .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isZero .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isNearest .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isCompatible .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isProcessDefined .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isUndefined .eqv. .false._LK
    end function test_constructRound_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructRound_5() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Round_type) :: round
        round = Round_type("Nearest")
        assertion = .not. round%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. round%val == "nearest"
        call test%assert(assertion)
        assertion = assertion .and. round%isUp .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isDown .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isZero .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isNearest .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isCompatible .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isProcessDefined .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isUndefined .eqv. .false._LK
    end function test_constructRound_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructRound_6() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Round_type) :: round
        round = Round_type("Nearest")
        assertion = .not. round%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. round%val == "nearest"
        call test%assert(assertion)
        assertion = assertion .and. round%isUp .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isDown .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isZero .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isNearest .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isCompatible .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isProcessDefined .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isUndefined .eqv. .false._LK
    end function test_constructRound_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructRound_7() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Round_type) :: round
        round = Round_type("Compatible")
        assertion = .not. round%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. round%val == "compatible"
        call test%assert(assertion)
        assertion = assertion .and. round%isUp .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isDown .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isZero .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isNearest .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isCompatible .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isProcessDefined .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isUndefined .eqv. .false._LK
    end function test_constructRound_7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructRound_8() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Round_type) :: round
        round = Round_type("Processor_defined")
        assertion = .not. round%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. round%val == "processor_defined"
        call test%assert(assertion)
        assertion = assertion .and. round%isUp .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isDown .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isZero .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isNearest .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isCompatible .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isProcessDefined .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. round%isUndefined .eqv. .false._LK
    end function test_constructRound_8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructRound_9() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Round_type) :: round
        round = Round_type("nonsense")
        assertion = round%err%occurred
       !call test%assert(assertion)
       !assertion = assertion .and. round%val == ""
       !call test%assert(assertion)
       !assertion = assertion .and. round%isUp .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. round%isDown .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. round%isZero .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. round%isNearest .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. round%isCompatible .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. round%isProcessDefined .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. round%isUndefined .eqv. .false._LK
    end function test_constructRound_9

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructSign_1() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Sign_type) :: sign
        sign = Sign_type()
        assertion = .not. sign%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. sign%val == "processor_defined"
        call test%assert(assertion)
        assertion = assertion .and. sign%isSuppress .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. sign%isPlus .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. sign%isProcessDefined .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. sign%isUndefined .eqv. .false._LK
    end function test_constructSign_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructSign_2() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Sign_type) :: sign
        sign = Sign_type("Undefined")
        assertion = .not. sign%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. sign%val == "undefined"
        call test%assert(assertion)
        assertion = assertion .and. sign%isSuppress .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. sign%isPlus .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. sign%isProcessDefined .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. sign%isUndefined .eqv. .true._LK
    end function test_constructSign_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructSign_3() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Sign_type) :: sign
        sign = Sign_type("Suppress")
        assertion = .not. sign%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. sign%val == "suppress"
        call test%assert(assertion)
        assertion = assertion .and. sign%isSuppress .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. sign%isPlus .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. sign%isProcessDefined .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. sign%isUndefined .eqv. .false._LK
    end function test_constructSign_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructSign_4() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Sign_type) :: sign
        sign = Sign_type("Plus")
        assertion = .not. sign%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. sign%val == "plus"
        call test%assert(assertion)
        assertion = assertion .and. sign%isSuppress .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. sign%isPlus .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. sign%isProcessDefined .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. sign%isUndefined .eqv. .false._LK
    end function test_constructSign_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructSign_5() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Sign_type) :: sign
        sign = Sign_type("Processor_defined")
        assertion = .not. sign%err%occurred
        call test%assert(assertion)
        assertion = assertion .and. sign%val == "processor_defined"
        call test%assert(assertion)
        assertion = assertion .and. sign%isSuppress .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. sign%isPlus .eqv. .false._LK
        call test%assert(assertion)
        assertion = assertion .and. sign%isProcessDefined .eqv. .true._LK
        call test%assert(assertion)
        assertion = assertion .and. sign%isUndefined .eqv. .false._LK
    end function test_constructSign_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructSign_6() result(assertion)
        implicit none
        logical(LK) :: assertion
        type(Sign_type) :: sign
        sign = Sign_type("nonsense")
        assertion = sign%err%occurred
       !call test%assert(assertion)
       !assertion = assertion .and. sign%err%occurred .eqv. .true._LK
       !call test%assert(assertion)
       !assertion = assertion .and. sign%val == ""
       !call test%assert(assertion)
       !assertion = assertion .and. sign%isSuppress .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. sign%isPlus .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. sign%isProcessDefined .eqv. .false._LK
       !call test%assert(assertion)
       !assertion = assertion .and. sign%isUndefined .eqv. .false._LK
    end function test_constructSign_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getWriteErr_1() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        integer :: stat
        stat = 1
        Err = getWriteErr(stat) ! LCOV_EXCL_LINE
        assertion = err%occurred
    end function test_getWriteErr_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getWriteErr_2() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        integer :: stat
        stat = 0
        Err = getWriteErr(stat) ! LCOV_EXCL_LINE
        assertion = .not. err%occurred
    end function test_getWriteErr_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getWriteErr_3() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        integer :: stat
        stat = -1
        Err = getWriteErr(stat) ! LCOV_EXCL_LINE
        assertion = err%occurred
    end function test_getWriteErr_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getReadErr_1() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        integer :: stat
        stat = 1
        Err = getReadErr(stat,"./path") ! LCOV_EXCL_LINE
        assertion = err%occurred
    end function test_getReadErr_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getReadErr_2() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        integer :: stat
        stat = 1
        Err = getReadErr(stat) ! LCOV_EXCL_LINE
        assertion = err%occurred
    end function test_getReadErr_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getReadErr_3() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        integer :: stat
        stat = 0
        Err = getReadErr(stat,"./path") ! LCOV_EXCL_LINE
        assertion = .not. err%occurred
    end function test_getReadErr_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getReadErr_4() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        integer :: stat
        stat = 0
        Err = getReadErr(stat) ! LCOV_EXCL_LINE
        assertion = .not. err%occurred
    end function test_getReadErr_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getReadErr_5() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        integer :: stat
        stat = -1
        Err = getReadErr(stat,"./path") ! LCOV_EXCL_LINE
        assertion = err%occurred
    end function test_getReadErr_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getReadErr_6() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        integer :: stat
        stat = -1
        Err = getReadErr(stat) ! LCOV_EXCL_LINE
        assertion = err%occurred
    end function test_getReadErr_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isCloseErr_1() result(assertion)
        logical(LK) :: assertion
        assertion = .not. isCloseErr(0_IK)
    end function test_isCloseErr_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isCloseErr_2() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        integer :: stat
        stat = 1
        Err = isCloseErr(stat) ! LCOV_EXCL_LINE
        assertion = err%occurred
    end function test_isCloseErr_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isCloseErr_3() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        Err = isCloseErr(-1) ! LCOV_EXCL_LINE
        assertion = .not. err%occurred
    end function test_isCloseErr_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isOpenErr_1() result(assertion)
        logical(LK) :: assertion
        assertion = .not. isOpenErr(0_IK)
    end function test_isOpenErr_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isOpenErr_2() result(assertion)
        logical(LK) :: assertion
        assertion = isOpenErr(1_IK)
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getOpenErr_3() result(assertion)
        logical(LK) :: assertion
        assertion = isOpenErr(-1_IK)
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isInqErr_1() result(assertion)
        logical(LK) :: assertion
        assertion = .not. isInqErr(0_IK)
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isInqErr_2() result(assertion)
        logical(LK) :: assertion
        assertion = isInqErr(1_IK)
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isInqErr_3() result(assertion)
        logical(LK) :: assertion
        assertion = isInqErr(-1_IK)
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_isExtant_1() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        logical(LK) :: exists
        type(err_type) :: Err
        call isExtant(exists, Err, unit = -1, file = "")
        assertion = err%occurred
    end

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_isExtant_2() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        logical(LK) :: exists
        type(err_type) :: Err
        call isExtant(exists, Err)
        assertion = err%occurred
    end function test_isExtant_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_isExtant_3() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        logical(LK) :: exists
        type(err_type) :: Err
        call isExtant(exists, Err, unit = 123)
        assertion = .not. err%occurred
    end function test_isExtant_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getOpenStatus_1() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        logical(LK) :: isOpen
        type(err_type) :: Err
        call getOpenStatus(isOpen, Err, unit = 124, file = "")
        assertion = err%occurred
    end function test_getOpenStatus_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getOpenStatus_2() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        logical(LK) :: isOpen
        type(err_type) :: Err
        call getOpenStatus(isOpen, Err)
        assertion = err%occurred
    end function test_getOpenStatus_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getOpenStatus_3() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        logical(LK) :: isOpen
        type(err_type) :: Err
        call getOpenStatus(isOpen, Err, unit = 123)
        assertion = .not. err%occurred
    end function test_getOpenStatus_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Check inquire with an input file name.
    function test_getOpenStatus_4() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        logical(LK) :: isOpen
        type(err_type) :: Err
        call getOpenStatus(isOpen, Err, file = "nonexisting.file")
        assertion = .not. isOpen .and. .not. err%occurred
    end function test_getOpenStatus_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getNumber_1() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        logical(LK) :: isNumbered
        integer :: number
        call getNumber(isNumbered, number, Err, unit = 124, file = "")
        assertion = err%occurred
    end function test_getNumber_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getNumber_2() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        logical(LK) :: isNumbered
        integer :: number
        call getNumber(isNumbered, number, Err)
        assertion = err%occurred
    end function test_getNumber_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getNumber_3() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        logical(LK) :: isNumbered
        integer :: number
        call getNumber(isNumbered, number, Err, unit = 124)
        assertion = .not. err%occurred .and. number == -1 .and. .not. isNumbered
    end function test_getNumber_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Check inquire with an input file name.
    function test_getNumber_4() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        logical(LK) :: isNumbered
        integer :: number
        call getNumber(isNumbered, number, Err, file = "nonexisting.file")
        assertion = .not. err%occurred .and. number == -1 .and. .not. isNumbered
    end function test_getNumber_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getName_1() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        logical(LK) :: isNamed
        character(:, SK), allocatable :: nameByCompiler
        call getName(isNamed, nameByCompiler, Err, unit = 124, file = "")
        assertion = err%occurred
    end function test_getName_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getName_2() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        logical(LK) :: isNamed
        character(:, SK), allocatable :: nameByCompiler
        call getName(isNamed, nameByCompiler, Err)
        assertion = err%occurred
    end function test_getName_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getName_3() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        logical(LK) :: isNamed
        character(:, SK), allocatable :: nameByCompiler
        call getName(isNamed, nameByCompiler, Err, unit = 124)
        assertion = .not. err%occurred .and. .not. isNamed
    end function test_getName_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Check inquire with an input file name.
    function test_getName_4() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        logical(LK) :: isNamed
        character(:, SK), allocatable :: nameByCompiler
        call getName(isNamed, nameByCompiler, Err, file = "nonexisting.file")
        assertion = .not. err%occurred .and. isNamed
#if !CODECOV_ENABLED
        assertion = .true._LK
#endif
    end function test_getName_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getAccess_1() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: access
        call getAccess(access, Err, unit = 124, file = "")
        assertion = err%occurred
    end function test_getAccess_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getAccess_2() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: access
        call getAccess(access, Err)
        assertion = err%occurred
    end function test_getAccess_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getAccess_3() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: access
        call getAccess(access, Err, unit = 124)
        assertion = .not. err%occurred
    end function test_getAccess_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Check inquire with an input file name.
    function test_getAccess_4() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: access
        call getAccess(access, Err, file = "nonexisting.file")
        assertion = .not. err%occurred
    end function test_getAccess_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getForm_1() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: form
        call getForm(form, Err, unit = 124, file = "")
        assertion = err%occurred
    end function test_getForm_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getForm_2() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: form
        call getForm(form, Err)
        assertion = err%occurred
    end function test_getForm_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getForm_3() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: form
        call getForm(form, Err, unit = 124)
        assertion = .not. err%occurred
    end function test_getForm_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Check inquire with an input file name.
    function test_getForm_4() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: form
        call getForm(form, Err, file = "nonexisting.file")
        assertion = .not. err%occurred
    end function test_getForm_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getRecl_1() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        integer :: recl
        call getRecl(recl, Err, unit = 124, file = "")
        assertion = err%occurred
    end function test_getRecl_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getRecl_2() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        integer :: recl
        call getRecl(recl, Err)
        assertion = err%occurred
    end function test_getRecl_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getRecl_3() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        integer :: recl
        call getRecl(recl, Err, unit = 124)
        assertion = .not. err%occurred
    end function test_getRecl_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Check inquire with an input file name.
    function test_getRecl_4() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        integer :: recl
        call getRecl(recl, Err, file = "nonexisting.file")
        assertion = .not. err%occurred
    end function test_getRecl_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getBlank_1() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: blank
        call getBlank(blank, Err, unit = 124, file = "")
        assertion = err%occurred
    end function test_getBlank_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getBlank_2() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: blank
        call getBlank(blank, Err)
        assertion = err%occurred
    end function test_getBlank_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getBlank_3() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: blank
        call getBlank(blank, Err, unit = 124)
        assertion = .not. err%occurred .and. blank == "undefined"
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))")   "err%occurred : ", err%occurred
            write(test%disp%unit,"(*(g0))")   "blank : ", blank
            write(test%disp%unit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getBlank_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Check inquire with an input file name.
    function test_getBlank_4() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: blank
        call getBlank(blank, Err, file = "nonexisting.file")
        assertion = .not. err%occurred .and. blank == "undefined"
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))")   "err%occurred : ", err%occurred
            write(test%disp%unit,"(*(g0))")   "blank : ", blank
            write(test%disp%unit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getBlank_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getPosition_1() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: position
        call getPosition(position, Err, unit = 124, file = "")
        assertion = err%occurred
    end function test_getPosition_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getPosition_2() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: position
        call getPosition(position, Err)
        assertion = err%occurred
    end function test_getPosition_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getPosition_3() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: position
        call getPosition(position, Err, unit = 124)
        assertion = .not. err%occurred .and. position == "undefined"
        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))")   "err%occurred : ", err%occurred
            write(test%disp%unit,"(*(g0))")   "position : ", position
            write(test%disp%unit,"(*(g0))")
            ! LCOV_EXCL_STOP
        end if
    end function test_getPosition_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Check inquire with an input file name.
    function test_getPosition_4() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: position
        call getPosition(position, Err, file = "nonexisting.file")
        assertion = .not. err%occurred .and. position == "undefined"
        if (test%traceable .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))")   "err%occurred : ", err%occurred
            write(test%disp%unit,"(*(g0))")   "position : ", position
            write(test%disp%unit,"(*(g0))")
            ! LCOV_EXCL_STOP
        end if
    end function test_getPosition_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getAction_1() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: action
        call getAction(action, Err, unit = 124, file = "")
        assertion = err%occurred
    end function test_getAction_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getAction_2() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: action
        call getAction(action, Err)
        assertion = err%occurred
    end function test_getAction_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getAction_3() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: action
        call getAction(action, Err, unit = 124)
        assertion = .not. err%occurred .and. action == "undefined"
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))")   "err%occurred : ", err%occurred
            write(test%disp%unit,"(*(g0))")   "action : ", action
            write(test%disp%unit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getAction_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Check inquire with an input file name.
    function test_getAction_4() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: action
        call getAction(action, Err, file = "nonexisting.file")
        assertion = .not. err%occurred !.and. action == "undefined"
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))")   "err%occurred : ", err%occurred
            write(test%disp%unit,"(*(g0))")   "action : ", action
            write(test%disp%unit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getAction_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getDelim_1() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: delim
        call getDelim(delim, Err, unit = 124, file = "")
        assertion = err%occurred
    end function test_getDelim_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getDelim_2() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: delim
        call getDelim(delim, Err)
        assertion = err%occurred
    end function test_getDelim_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getDelim_3() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: delim
        call getDelim(delim, Err, unit = 124)
        assertion = .not. err%occurred .and. delim == "undefined"
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))")   "err%occurred : ", err%occurred
            write(test%disp%unit,"(*(g0))")   "delim : ", delim
            write(test%disp%unit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getDelim_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !> Check inquire with an input file name.
    function test_getDelim_4() result(assertion)
        use pm_err, only: err_type
        implicit none
        logical(LK) :: assertion
        type(err_type) :: Err
        character(:, SK), allocatable :: delim
        call getDelim(delim, Err, file = "nonexisting.file")
        assertion = .not. err%occurred !.and. delim == "undefined"
        if (test%traceable .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(test%disp%unit,"(*(g0))")
            write(test%disp%unit,"(*(g0))")   "err%occurred : ", err%occurred
            write(test%disp%unit,"(*(g0))")   "delim : ", delim
            write(test%disp%unit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getDelim_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_io ! LCOV_EXCL_LINE