!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!   MIT License
!!!!
!!!!   ParaMonte: plain powerful parallel Monte Carlo library.
!!!!
!!!!   Copyright (C) 2012-present, The Computational Data Science Lab
!!!!
!!!!   This file is part of the ParaMonte library.
!!!!
!!!!   Permission is hereby granted, free of charge, to any person obtaining a
!!!!   copy of this software and associated documentation files (the "Software"),
!!!!   to deal in the Software without restriction, including without limitation
!!!!   the rights to use, copy, modify, merge, publish, distribute, sublicense,
!!!!   and/or sell copies of the Software, and to permit persons to whom the
!!!!   Software is furnished to do so, subject to the following conditions:
!!!!
!!!!   The above copyright notice and this permission notice shall be
!!!!   included in all copies or substantial portions of the Software.
!!!!
!!!!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!!!!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!!!!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!!!!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
!!!!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
!!!!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
!!!!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!!!!
!!!!   ACKNOWLEDGMENT
!!!!
!!!!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!!!!   As per the ParaMonte library license agreement terms, if you use any parts of
!!!!   this library for any purposes, kindly acknowledge the use of ParaMonte in your
!!!!   work (education/research/industry/development/...) by citing the ParaMonte
!!!!   library as described on this page:
!!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This module contains tests of the module [File_mod](@ref file_mod).
!>  \author Amir Shahmoradi

module Test_File_mod

    use File_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_File

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_File()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_getRecl_1, "test_getRecl_1")
        call Test%run(test_getRecl_2, "test_getRecl_2")
        call Test%run(test_getRecl_3, "test_getRecl_3")
        call Test%run(test_getRecl_4, "test_getRecl_4")
        call Test%run(test_getForm_1, "test_getForm_1")
        call Test%run(test_getForm_2, "test_getForm_2")
        call Test%run(test_getForm_3, "test_getForm_3")
        call Test%run(test_getForm_4, "test_getForm_4")
        call Test%run(test_getName_1, "test_getName_1")
        call Test%run(test_getName_2, "test_getName_2")
        call Test%run(test_getName_3, "test_getName_3")
        call Test%run(test_getName_4, "test_getName_4")
        call Test%run(test_getBlank_1, "test_getBlank_1")
        call Test%run(test_getBlank_2, "test_getBlank_2")
        call Test%run(test_getBlank_3, "test_getBlank_3")
        call Test%run(test_getBlank_4, "test_getBlank_4")
        call Test%run(test_getDelim_1, "test_getDelim_1")
        call Test%run(test_getDelim_2, "test_getDelim_2")
        call Test%run(test_getDelim_3, "test_getDelim_3")
        call Test%run(test_getDelim_4, "test_getDelim_4")
        call Test%run(test_getAction_1, "test_getAction_1")
        call Test%run(test_getAction_2, "test_getAction_2")
        call Test%run(test_getAction_3, "test_getAction_3")
        call Test%run(test_getAction_4, "test_getAction_4")
        call Test%run(test_getAccess_1, "test_getAccess_1")
        call Test%run(test_getAccess_2, "test_getAccess_2")
        call Test%run(test_getAccess_3, "test_getAccess_3")
        call Test%run(test_getAccess_4, "test_getAccess_4")
        call Test%run(test_getNumber_1, "test_getNumber_1")
        call Test%run(test_getNumber_2, "test_getNumber_2")
        call Test%run(test_getNumber_3, "test_getNumber_3")
        call Test%run(test_getNumber_4, "test_getNumber_4")
        call Test%run(test_getPosition_1, "test_getPosition_1")
        call Test%run(test_getPosition_2, "test_getPosition_2")
        call Test%run(test_getPosition_3, "test_getPosition_3")
        call Test%run(test_getPosition_4, "test_getPosition_4")
        call Test%run(test_constructPad_1, "test_constructPad_1")
        call Test%run(test_constructPad_2, "test_constructPad_2")
        call Test%run(test_constructPad_3, "test_constructPad_3")
        call Test%run(test_constructPad_4, "test_constructPad_4")
        call Test%run(test_constructPad_5, "test_constructPad_5")
        call Test%run(test_constructFile_1, "test_constructFile_1")
        call Test%run(test_constructFile_2, "test_constructFile_2")
        call Test%run(test_constructFile_3, "test_constructFile_3")
        call Test%run(test_constructForm_1, "test_constructForm_1")
        call Test%run(test_constructForm_2, "test_constructForm_2")
        call Test%run(test_constructForm_3, "test_constructForm_3")
        call Test%run(test_constructForm_4, "test_constructForm_4")
        call Test%run(test_constructForm_5, "test_constructForm_5")
        call Test%run(test_constructSign_1, "test_constructSign_1")
        call Test%run(test_constructSign_2, "test_constructSign_2")
        call Test%run(test_constructSign_3, "test_constructSign_3")
        call Test%run(test_constructSign_4, "test_constructSign_4")
        call Test%run(test_constructSign_5, "test_constructSign_5")
        call Test%run(test_constructSign_6, "test_constructSign_6")
        call Test%run(test_getOpenStatus_1, "test_getOpenStatus_1")
        call Test%run(test_getOpenStatus_2, "test_getOpenStatus_2")
        call Test%run(test_getOpenStatus_3, "test_getOpenStatus_3")
        call Test%run(test_getOpenStatus_4, "test_getOpenStatus_4")
        call Test%run(test_getExistStatus_1, "test_getExistStatus_1")
        call Test%run(test_getExistStatus_2, "test_getExistStatus_2")
        call Test%run(test_getExistStatus_3, "test_getExistStatus_3")
        call Test%run(test_constructBlank_1, "test_constructBlank_1")
        call Test%run(test_constructBlank_2, "test_constructBlank_2")
        call Test%run(test_constructBlank_3, "test_constructBlank_3")
        call Test%run(test_constructBlank_4, "test_constructBlank_4")
        call Test%run(test_constructBlank_5, "test_constructBlank_5")
        call Test%run(test_constructDelim_1, "test_constructDelim_1")
        call Test%run(test_constructDelim_2, "test_constructDelim_2")
        call Test%run(test_constructDelim_3, "test_constructDelim_3")
        call Test%run(test_constructDelim_4, "test_constructDelim_4")
        call Test%run(test_constructDelim_5, "test_constructDelim_5")
        call Test%run(test_constructDelim_6, "test_constructDelim_6")
        call Test%run(test_constructRound_1, "test_constructRound_1")
        call Test%run(test_constructRound_2, "test_constructRound_2")
        call Test%run(test_constructRound_3, "test_constructRound_3")
        call Test%run(test_constructRound_4, "test_constructRound_4")
        call Test%run(test_constructRound_5, "test_constructRound_5")
        call Test%run(test_constructRound_6, "test_constructRound_6")
        call Test%run(test_constructRound_7, "test_constructRound_7")
        call Test%run(test_constructRound_8, "test_constructRound_8")
        call Test%run(test_constructRound_9, "test_constructRound_9")
        call Test%run(test_constructAction_1, "test_constructAction_1")
        call Test%run(test_constructAction_2, "test_constructAction_2")
        call Test%run(test_constructAction_3, "test_constructAction_3")
        call Test%run(test_constructAction_4, "test_constructAction_4")
        call Test%run(test_constructAction_5, "test_constructAction_5")
        call Test%run(test_constructAccess_1, "test_constructAccess_1")
        call Test%run(test_constructAccess_2, "test_constructAccess_2")
        call Test%run(test_constructAccess_3, "test_constructAccess_3")
        call Test%run(test_constructAccess_4, "test_constructAccess_4")
        call Test%run(test_constructAccess_5, "test_constructAccess_5")
        call Test%run(test_constructPosition_1, "test_constructPosition_1")
        call Test%run(test_constructPosition_2, "test_constructPosition_2")
        call Test%run(test_constructPosition_3, "test_constructPosition_3")
        call Test%run(test_constructPosition_4, "test_constructPosition_4")
        call Test%run(test_constructPosition_5, "test_constructPosition_5")
        call Test%run(test_constructPosition_6, "test_constructPosition_6")
        call Test%run(test_getCloseErr_1, "test_getCloseErr_1")
        call Test%run(test_getCloseErr_2, "test_getCloseErr_2")
        call Test%run(test_getCloseErr_3, "test_getCloseErr_3")
        call Test%run(test_getWriteErr_1, "test_getWriteErr_1")
        call Test%run(test_getWriteErr_2, "test_getWriteErr_2")
        call Test%run(test_getWriteErr_3, "test_getWriteErr_3")
        call Test%run(test_getOpenErr_1, "test_getOpenErr_1")
        call Test%run(test_getOpenErr_2, "test_getOpenErr_2")
        call Test%run(test_getOpenErr_3, "test_getOpenErr_3")
        call Test%run(test_getReadErr_1, "test_getReadErr_1")
        call Test%run(test_getReadErr_2, "test_getReadErr_2")
        call Test%run(test_getReadErr_3, "test_getReadErr_3")
        call Test%run(test_getReadErr_4, "test_getReadErr_4")
        call Test%run(test_getReadErr_5, "test_getReadErr_5")
        call Test%run(test_getReadErr_6, "test_getReadErr_6")
        call Test%run(test_getInqErr_1, "test_getInqErr_1")
        call Test%run(test_getInqErr_2, "test_getInqErr_2")
        call Test%run(test_getInqErr_3, "test_getInqErr_3")
        call Test%finalize()

    end subroutine test_File

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructFile_1() result(assertion)

        implicit none
        logical :: assertion
        type(File_type) :: File

        File = File_type()
        assertion = .not. File%Err%occurred; if (.not. assertion) return

        assertion = assertion .and. File%Path%original  == ""
        assertion = assertion .and. File%Path%modified  == ""
        assertion = assertion .and. File%Path%dir       == ""
        assertion = assertion .and. File%Path%name      == ""
        assertion = assertion .and. File%Path%ext       == ""
        assertion = assertion .and. File%unit           == -2147483647
        assertion = assertion .and. File%number         == -2147483647
        assertion = assertion .and. File%recl           == -2147483647
        assertion = assertion .and. File%exists         .eqv. .false.
        assertion = assertion .and. File%isOpen         .eqv. .false.
        assertion = assertion .and. File%isNamed        .eqv. .false.
        assertion = assertion .and. File%isNumbered     .eqv. .false.
        assertion = assertion .and. File%status         == "unknown"
        assertion = assertion .and. File%asynchronous   == "no"
        assertion = assertion .and. File%format         == ""
        assertion = assertion .and. File%nameByCompiler == ""
        assertion = assertion .and. File%Action%value   == "readwrite"
        assertion = assertion .and. File%Access%value   == "sequential"
        assertion = assertion .and. File%Form%value     == "formatted"
        assertion = assertion .and. File%Blank%value    == "null"
        assertion = assertion .and. File%Position%value == "asis"
        assertion = assertion .and. File%Delim%value    == "none"
        assertion = assertion .and. File%Pad%value      == "yes"
        assertion = assertion .and. File%Round%value    == "processor_defined"
        assertion = assertion .and. File%Sign%value     == "processor_defined"
        assertion = assertion .and. File%Err%occurred   .eqv. .false.
        assertion = assertion .and. File%Err%msg        == ""

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "File%Path%original  : ", File%Path%original
            write(Test%outputUnit,"(*(g0))")   "File%Path%modified  : ", File%Path%modified
            write(Test%outputUnit,"(*(g0))")   "File%Path%dir       : ", File%Path%dir
            write(Test%outputUnit,"(*(g0))")   "File%Path%name      : ", File%Path%name
            write(Test%outputUnit,"(*(g0))")   "File%Path%ext       : ", File%Path%ext
            write(Test%outputUnit,"(*(g0))")   "File%Path%shellSlash: ", File%Path%shellSlash
            write(Test%outputUnit,"(*(g0))")   "File%unit           : ", File%unit
            write(Test%outputUnit,"(*(g0))")   "File%number         : ", File%number
            write(Test%outputUnit,"(*(g0))")   "File%recl           : ", File%recl
            write(Test%outputUnit,"(*(g0))")   "File%exists         : ", File%exists
            write(Test%outputUnit,"(*(g0))")   "File%isOpen         : ", File%isOpen
            write(Test%outputUnit,"(*(g0))")   "File%isNamed        : ", File%isNamed
            write(Test%outputUnit,"(*(g0))")   "File%isNumbered     : ", File%isNumbered
            write(Test%outputUnit,"(*(g0))")   "File%status         : ", File%status
            write(Test%outputUnit,"(*(g0))")   "File%asynchronous   : ", File%asynchronous
            write(Test%outputUnit,"(*(g0))")   "File%format         : ", File%format
            write(Test%outputUnit,"(*(g0))")   "File%nameByCompiler : ", File%nameByCompiler
            write(Test%outputUnit,"(*(g0))")   "File%Action%value   : ", File%Action%value
            write(Test%outputUnit,"(*(g0))")   "File%Access%value   : ", File%Access%value
            write(Test%outputUnit,"(*(g0))")   "File%Form%value     : ", File%Form%value
            write(Test%outputUnit,"(*(g0))")   "File%Blank%value    : ", File%Blank%value
            write(Test%outputUnit,"(*(g0))")   "File%Position%value : ", File%Position%value
            write(Test%outputUnit,"(*(g0))")   "File%Delim%value    : ", File%Delim%value
            write(Test%outputUnit,"(*(g0))")   "File%Pad%value      : ", File%Pad%value
            write(Test%outputUnit,"(*(g0))")   "File%Round%value    : ", File%Round%value
            write(Test%outputUnit,"(*(g0))")   "File%Sign%value     : ", File%Sign%value
            write(Test%outputUnit,"(*(g0))")   "File%Err%occurred   : ", File%Err%occurred
            write(Test%outputUnit,"(*(g0))")   "File%Err%msg        : ", File%Err%msg
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructFile_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructFile_2() result(assertion)

        implicit none
        logical :: assertion
        type(File_type) :: File

        File = File_type( unit=13 &
                        , recl=9999 &
                        , path="./test_File_mod/\test_File_mod\-" &
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

        assertion = .not. File%Err%occurred; if (.not. assertion) return

        assertion = assertion .and. File%Path%original  == "./test_File_mod/\test_File_mod\-"

        if (File%Path%shellSlash=="/") then
            assertion = assertion .and. File%Path%modified  == "./test_File_mod/\test_File_mod\-"
            assertion = assertion .and. File%Path%dir       == "./test_File_mod/"
            assertion = assertion .and. File%Path%name      == "\test_File_mod\-"
            assertion = assertion .and. File%Path%ext       == ""
#if defined OS_IS_WINDOWS
        else
            assertion = assertion .and. File%Path%modified  == ".\test_File_mod\\test_File_mod\-"
            assertion = assertion .and. File%Path%dir       == ".\test_File_mod\\test_File_mod\"
            assertion = assertion .and. File%Path%name      == "-"
            assertion = assertion .and. File%Path%ext       == ""
#endif
        end if
        assertion = assertion .and. File%unit           == 13
        assertion = assertion .and. File%number         == -2147483647
        assertion = assertion .and. File%recl           == 9999
        assertion = assertion .and. File%exists         .eqv. .false.
        assertion = assertion .and. File%isOpen         .eqv. .false.
        assertion = assertion .and. File%isNamed        .eqv. .false.
        assertion = assertion .and. File%isNumbered     .eqv. .false.
        assertion = assertion .and. File%status         == "new"
        assertion = assertion .and. File%asynchronous   == "yes"
        assertion = assertion .and. File%format         == "(A)"
        assertion = assertion .and. File%nameByCompiler == ""
        assertion = assertion .and. File%Action%value   == "read"
        assertion = assertion .and. File%Access%value   == "direct"
        assertion = assertion .and. File%Form%value     == "unformatted"
        assertion = assertion .and. File%Blank%value    == "undefined"
        assertion = assertion .and. File%Position%value == "append"
        assertion = assertion .and. File%Delim%value    == "quote"
        assertion = assertion .and. File%Pad%value      == "no"
        assertion = assertion .and. File%Round%value    == "up"
        assertion = assertion .and. File%Sign%value     == "undefined"
        assertion = assertion .and. File%Err%occurred   .eqv. .false.
        assertion = assertion .and. File%Err%msg        == ""

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "File%Path%original  : ", File%Path%original
            write(Test%outputUnit,"(*(g0))")   "File%Path%modified  : ", File%Path%modified
            write(Test%outputUnit,"(*(g0))")   "File%Path%dir       : ", File%Path%dir
            write(Test%outputUnit,"(*(g0))")   "File%Path%name      : ", File%Path%name
            write(Test%outputUnit,"(*(g0))")   "File%Path%ext       : ", File%Path%ext
            write(Test%outputUnit,"(*(g0))")   "File%Path%shellSlash: ", File%Path%shellSlash
            write(Test%outputUnit,"(*(g0))")   "File%unit           : ", File%unit
            write(Test%outputUnit,"(*(g0))")   "File%number         : ", File%number
            write(Test%outputUnit,"(*(g0))")   "File%recl           : ", File%recl
            write(Test%outputUnit,"(*(g0))")   "File%exists         : ", File%exists
            write(Test%outputUnit,"(*(g0))")   "File%isOpen         : ", File%isOpen
            write(Test%outputUnit,"(*(g0))")   "File%isNamed        : ", File%isNamed
            write(Test%outputUnit,"(*(g0))")   "File%isNumbered     : ", File%isNumbered
            write(Test%outputUnit,"(*(g0))")   "File%status         : ", File%status
            write(Test%outputUnit,"(*(g0))")   "File%asynchronous   : ", File%asynchronous
            write(Test%outputUnit,"(*(g0))")   "File%format         : ", File%format
            write(Test%outputUnit,"(*(g0))")   "File%nameByCompiler : ", File%nameByCompiler
            write(Test%outputUnit,"(*(g0))")   "File%Action%value   : ", File%Action%value
            write(Test%outputUnit,"(*(g0))")   "File%Access%value   : ", File%Access%value
            write(Test%outputUnit,"(*(g0))")   "File%Form%value     : ", File%Form%value
            write(Test%outputUnit,"(*(g0))")   "File%Blank%value    : ", File%Blank%value
            write(Test%outputUnit,"(*(g0))")   "File%Position%value : ", File%Position%value
            write(Test%outputUnit,"(*(g0))")   "File%Delim%value    : ", File%Delim%value
            write(Test%outputUnit,"(*(g0))")   "File%Pad%value      : ", File%Pad%value
            write(Test%outputUnit,"(*(g0))")   "File%Round%value    : ", File%Round%value
            write(Test%outputUnit,"(*(g0))")   "File%Sign%value     : ", File%Sign%value
            write(Test%outputUnit,"(*(g0))")   "File%Err%occurred   : ", File%Err%occurred
            write(Test%outputUnit,"(*(g0))")   "File%Err%msg        : ", File%Err%msg
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructFile_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Test the effects of missing input argument `form`.
    function test_constructFile_3() result(assertion)

        implicit none
        logical :: assertion
        type(File_type) :: File

        File = File_type( unit=13 &
                        , recl=9999 &
                        , path="./test_File_mod/\test_File_mod\-" &
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

        assertion = .not. File%Err%occurred; if (.not. assertion) return

        assertion = assertion .and. File%Path%original  == "./test_File_mod/\test_File_mod\-"

        if (File%Path%shellSlash=="/") then
            assertion = assertion .and. File%Path%modified  == "./test_File_mod/\test_File_mod\-"
            assertion = assertion .and. File%Path%dir       == "./test_File_mod/"
            assertion = assertion .and. File%Path%name      == "\test_File_mod\-"
            assertion = assertion .and. File%Path%ext       == ""
#if defined OS_IS_WINDOWS
        else
            assertion = assertion .and. File%Path%modified  == ".\test_File_mod\\test_File_mod\-"
            assertion = assertion .and. File%Path%dir       == ".\test_File_mod\\test_File_mod\"
            assertion = assertion .and. File%Path%name      == "-"
            assertion = assertion .and. File%Path%ext       == ""
#endif
        end if

        assertion = assertion .and. File%unit           == 13
        assertion = assertion .and. File%number         == -2147483647
        assertion = assertion .and. File%recl           == 9999
        assertion = assertion .and. File%exists         .eqv. .false.
        assertion = assertion .and. File%isOpen         .eqv. .false.
        assertion = assertion .and. File%isNamed        .eqv. .false.
        assertion = assertion .and. File%isNumbered     .eqv. .false.
        assertion = assertion .and. File%status         == "new"
        assertion = assertion .and. File%asynchronous   == "yes"
        assertion = assertion .and. File%format         == "(A)"
        assertion = assertion .and. File%nameByCompiler == ""
        assertion = assertion .and. File%Action%value   == "read"
        assertion = assertion .and. File%Access%value   == "direct"
        assertion = assertion .and. File%Form%value     == "unformatted"
        assertion = assertion .and. File%Blank%value    == "undefined"
        assertion = assertion .and. File%Position%value == "append"
        assertion = assertion .and. File%Delim%value    == "quote"
        assertion = assertion .and. File%Pad%value      == "no"
        assertion = assertion .and. File%Round%value    == "up"
        assertion = assertion .and. File%Sign%value     == "undefined"
        assertion = assertion .and. File%Err%occurred   .eqv. .false.
        assertion = assertion .and. File%Err%msg        == ""

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "File%Path%original  : ", File%Path%original
            write(Test%outputUnit,"(*(g0))")   "File%Path%modified  : ", File%Path%modified
            write(Test%outputUnit,"(*(g0))")   "File%Path%dir       : ", File%Path%dir
            write(Test%outputUnit,"(*(g0))")   "File%Path%name      : ", File%Path%name
            write(Test%outputUnit,"(*(g0))")   "File%Path%ext       : ", File%Path%ext
            write(Test%outputUnit,"(*(g0))")   "File%Path%shellSlash: ", File%Path%shellSlash
            write(Test%outputUnit,"(*(g0))")   "File%unit           : ", File%unit
            write(Test%outputUnit,"(*(g0))")   "File%number         : ", File%number
            write(Test%outputUnit,"(*(g0))")   "File%recl           : ", File%recl
            write(Test%outputUnit,"(*(g0))")   "File%exists         : ", File%exists
            write(Test%outputUnit,"(*(g0))")   "File%isOpen         : ", File%isOpen
            write(Test%outputUnit,"(*(g0))")   "File%isNamed        : ", File%isNamed
            write(Test%outputUnit,"(*(g0))")   "File%isNumbered     : ", File%isNumbered
            write(Test%outputUnit,"(*(g0))")   "File%status         : ", File%status
            write(Test%outputUnit,"(*(g0))")   "File%asynchronous   : ", File%asynchronous
            write(Test%outputUnit,"(*(g0))")   "File%format         : ", File%format
            write(Test%outputUnit,"(*(g0))")   "File%nameByCompiler : ", File%nameByCompiler
            write(Test%outputUnit,"(*(g0))")   "File%Action%value   : ", File%Action%value
            write(Test%outputUnit,"(*(g0))")   "File%Access%value   : ", File%Access%value
            write(Test%outputUnit,"(*(g0))")   "File%Form%value     : ", File%Form%value
            write(Test%outputUnit,"(*(g0))")   "File%Blank%value    : ", File%Blank%value
            write(Test%outputUnit,"(*(g0))")   "File%Position%value : ", File%Position%value
            write(Test%outputUnit,"(*(g0))")   "File%Delim%value    : ", File%Delim%value
            write(Test%outputUnit,"(*(g0))")   "File%Pad%value      : ", File%Pad%value
            write(Test%outputUnit,"(*(g0))")   "File%Round%value    : ", File%Round%value
            write(Test%outputUnit,"(*(g0))")   "File%Sign%value     : ", File%Sign%value
            write(Test%outputUnit,"(*(g0))")   "File%Err%occurred   : ", File%Err%occurred
            write(Test%outputUnit,"(*(g0))")   "File%Err%msg        : ", File%Err%msg
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_constructFile_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructPad_1() result(assertion)
        implicit none
        logical :: assertion
        type(Pad_type) :: Pad
        Pad = Pad_type()
        assertion = .not. Pad%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Pad%value == "yes"
        assertion = assertion .and. Pad%isPadded .eqv. .true.
        assertion = assertion .and. Pad%isNotPadded .eqv. .false.
        assertion = assertion .and. Pad%isUndefined .eqv. .false.
    end function test_constructPad_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Testing `Pad_type( value = 'undefined' )`
    function test_constructPad_2() result(assertion)
        implicit none
        logical :: assertion
        type(Pad_type) :: Pad
        Pad = Pad_type(value='undefined')
        assertion = .not. Pad%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Pad%value == "undefined"
        assertion = assertion .and. Pad%isPadded .eqv. .false.
        assertion = assertion .and. Pad%isNotPadded .eqv. .false.
        assertion = assertion .and. Pad%isUndefined .eqv. .true.
    end function test_constructPad_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Testing `Pad_type( value = 'yes' )`
    function test_constructPad_3() result(assertion)
        implicit none
        logical :: assertion
        type(Pad_type) :: Pad
        Pad = Pad_type(value='yes')
        assertion = .not. Pad%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Pad%value == "yes"
        assertion = assertion .and. Pad%isPadded .eqv. .true.
        assertion = assertion .and. Pad%isNotPadded .eqv. .false.
        assertion = assertion .and. Pad%isUndefined .eqv. .false.
    end function test_constructPad_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Testing `Pad_type( value = 'no' )`
    function test_constructPad_4() result(assertion)
        implicit none
        logical :: assertion
        type(Pad_type) :: Pad
        Pad = Pad_type(value='no')
        assertion = .not. Pad%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Pad%value == "no"
        assertion = assertion .and. Pad%isPadded .eqv. .false.
        assertion = assertion .and. Pad%isNotPadded .eqv. .true.
        assertion = assertion .and. Pad%isUndefined .eqv. .false.
    end function test_constructPad_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Testing `Pad_type( value = 'nonsense' )`
    function test_constructPad_5() result(assertion)
        implicit none
        logical :: assertion
        type(Pad_type) :: Pad
        Pad = Pad_type(value="nonsense")
        assertion = Pad%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Pad%value == ""
        assertion = assertion .and. Pad%isPadded .eqv. .false.
        assertion = assertion .and. Pad%isNotPadded .eqv. .false.
        assertion = assertion .and. Pad%isUndefined .eqv. .false.
    end function test_constructPad_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAction_1() result(assertion)
        implicit none
        logical :: assertion
        type(Action_type) :: Action
        Action = Action_type()
        assertion = .not. Action%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Action%value == "readwrite"
        assertion = assertion .and. Action%isRead .eqv. .false.
        assertion = assertion .and. Action%isWrite .eqv. .false.
        assertion = assertion .and. Action%isReadWrite .eqv. .true.
        assertion = assertion .and. Action%isUndefined .eqv. .false.
    end function test_constructAction_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAction_2() result(assertion)
        implicit none
        logical :: assertion
        type(Action_type) :: Action
        Action = Action_type(value="undefined")
        assertion = .not. Action%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Action%value == "readwrite"
        assertion = assertion .and. Action%isRead .eqv. .false.
        assertion = assertion .and. Action%isWrite .eqv. .false.
        assertion = assertion .and. Action%isReadWrite .eqv. .false.
        assertion = assertion .and. Action%isUndefined .eqv. .true.
    end function test_constructAction_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAction_3() result(assertion)
        implicit none
        logical :: assertion
        type(Action_type) :: Action
        Action = Action_type(value="read")
        assertion = .not. Action%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Action%value == "read"
        assertion = assertion .and. Action%isRead .eqv. .true.
        assertion = assertion .and. Action%isWrite .eqv. .false.
        assertion = assertion .and. Action%isReadWrite .eqv. .false.
        assertion = assertion .and. Action%isUndefined .eqv. .false.
    end function test_constructAction_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAction_4() result(assertion)
        implicit none
        logical :: assertion
        type(Action_type) :: Action
        Action = Action_type(value="write")
        assertion = .not. Action%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Action%value == "write"
        assertion = assertion .and. Action%isRead .eqv. .false.
        assertion = assertion .and. Action%isWrite .eqv. .true.
        assertion = assertion .and. Action%isReadWrite .eqv. .false.
        assertion = assertion .and. Action%isUndefined .eqv. .false.
    end function test_constructAction_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAction_5() result(assertion)
        implicit none
        logical :: assertion
        type(Action_type) :: Action
        Action = Action_type(value="nonsense")
        assertion = Action%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Action%value == ""
        assertion = assertion .and. Action%isRead .eqv. .false.
        assertion = assertion .and. Action%isWrite .eqv. .false.
        assertion = assertion .and. Action%isReadWrite .eqv. .false.
        assertion = assertion .and. Action%isUndefined .eqv. .false.
    end function test_constructAction_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAccess_1() result(assertion)
        implicit none
        logical :: assertion
        type(Access_type) :: Access
        Access = Access_type()
        assertion = .not. Access%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Access%value == "sequential"
        assertion = assertion .and. Access%isSequential .eqv. .true.
        assertion = assertion .and. Access%isDirect .eqv. .false.
        assertion = assertion .and. Access%isUndefined .eqv. .false.
    end function test_constructAccess_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAccess_2() result(assertion)
        implicit none
        logical :: assertion
        type(Access_type) :: Access
        Access = Access_type("undefined")
        assertion = .not. Access%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Access%value == "undefined"
        assertion = assertion .and. Access%isSequential .eqv. .false.
        assertion = assertion .and. Access%isDirect .eqv. .false.
        assertion = assertion .and. Access%isUndefined .eqv. .true.
    end function test_constructAccess_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAccess_3() result(assertion)
        implicit none
        logical :: assertion
        type(Access_type) :: Access
        Access = Access_type("sequential")
        assertion = .not. Access%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Access%value == "sequential"
        assertion = assertion .and. Access%isSequential .eqv. .true.
        assertion = assertion .and. Access%isDirect .eqv. .false.
        assertion = assertion .and. Access%isUndefined .eqv. .false.
    end function test_constructAccess_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAccess_4() result(assertion)
        implicit none
        logical :: assertion
        type(Access_type) :: Access
        Access = Access_type("direct")
        assertion = .not. Access%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Access%value == "direct"
        assertion = assertion .and. Access%isSequential .eqv. .false.
        assertion = assertion .and. Access%isDirect .eqv. .true.
        assertion = assertion .and. Access%isUndefined .eqv. .false.
    end function test_constructAccess_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructAccess_5() result(assertion)
        implicit none
        logical :: assertion
        type(Access_type) :: Access
        Access = Access_type("nonsense")
        assertion = Access%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Access%value == "nonsense"
        assertion = assertion .and. Access%isSequential .eqv. .false.
        assertion = assertion .and. Access%isDirect .eqv. .false.
        assertion = assertion .and. Access%isUndefined .eqv. .false.
    end function test_constructAccess_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructForm_1() result(assertion)
        implicit none
        logical :: assertion
        type(Form_type) :: Form
        Form = Form_type()
        assertion = .not. Form%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Form%value == "formatted"
        assertion = assertion .and. Form%isFormatted .eqv. .true.
        assertion = assertion .and. Form%isUnformatted .eqv. .false.
        assertion = assertion .and. Form%isUndefined .eqv. .false.
    end function test_constructForm_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructForm_2() result(assertion)
        implicit none
        logical :: assertion
        type(Form_type) :: Form
        Form = Form_type("undefined")
        assertion = .not. Form%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Form%value == "undefined"
        assertion = assertion .and. Form%isFormatted .eqv. .false.
        assertion = assertion .and. Form%isUnformatted .eqv. .false.
        assertion = assertion .and. Form%isUndefined .eqv. .true.
    end function test_constructForm_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructForm_3() result(assertion)
        implicit none
        logical :: assertion
        type(Form_type) :: Form
        Form = Form_type("formatted")
        assertion = .not. Form%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Form%value == "formatted"
        assertion = assertion .and. Form%isFormatted .eqv. .true.
        assertion = assertion .and. Form%isUnformatted .eqv. .false.
        assertion = assertion .and. Form%isUndefined .eqv. .false.
    end function test_constructForm_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructForm_4() result(assertion)
        implicit none
        logical :: assertion
        type(Form_type) :: Form
        Form = Form_type("unformatted")
        assertion = .not. Form%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Form%value == "unformatted"
        assertion = assertion .and. Form%isFormatted .eqv. .false.
        assertion = assertion .and. Form%isUnformatted .eqv. .true.
        assertion = assertion .and. Form%isUndefined .eqv. .false.
    end function test_constructForm_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructForm_5() result(assertion)
        implicit none
        logical :: assertion
        type(Form_type) :: Form
        Form = Form_type("nonsense")
        assertion = Form%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Form%value == "nonsense"
        assertion = assertion .and. Form%isFormatted .eqv. .false.
        assertion = assertion .and. Form%isUnformatted .eqv. .false.
        assertion = assertion .and. Form%isUndefined .eqv. .false.
    end function test_constructForm_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructBlank_1() result(assertion)
        implicit none
        logical :: assertion
        type(Blank_type) :: Blank
        Blank = Blank_type()
        assertion = .not. Blank%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Blank%value == "null"
        assertion = assertion .and. Blank%isNull .eqv. .true.
        assertion = assertion .and. Blank%isZero .eqv. .false.
        assertion = assertion .and. Blank%isUndefined .eqv. .false.
    end function test_constructBlank_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructBlank_2() result(assertion)
        implicit none
        logical :: assertion
        type(Blank_type) :: Blank
        Blank = Blank_type("undefined")
        assertion = .not. Blank%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Blank%value == "undefined"
        assertion = assertion .and. Blank%isNull .eqv. .false.
        assertion = assertion .and. Blank%isZero .eqv. .false.
        assertion = assertion .and. Blank%isUndefined .eqv. .true.
    end function test_constructBlank_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructBlank_3() result(assertion)
        implicit none
        logical :: assertion
        type(Blank_type) :: Blank
        Blank = Blank_type("null")
        assertion = .not. Blank%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Blank%value == "null"
        assertion = assertion .and. Blank%isNull .eqv. .true.
        assertion = assertion .and. Blank%isZero .eqv. .false.
        assertion = assertion .and. Blank%isUndefined .eqv. .false.
    end function test_constructBlank_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructBlank_4() result(assertion)
        implicit none
        logical :: assertion
        type(Blank_type) :: Blank
        Blank = Blank_type("zero")
        assertion = .not. Blank%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Blank%value == "zero"
        assertion = assertion .and. Blank%isNull .eqv. .false.
        assertion = assertion .and. Blank%isZero .eqv. .true.
        assertion = assertion .and. Blank%isUndefined .eqv. .false.
    end function test_constructBlank_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructBlank_5() result(assertion)
        implicit none
        logical :: assertion
        type(Blank_type) :: Blank
        Blank = Blank_type("nonsense")
        assertion = Blank%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Blank%value == "nonsense"
        assertion = assertion .and. Blank%isNull .eqv. .false.
        assertion = assertion .and. Blank%isZero .eqv. .false.
        assertion = assertion .and. Blank%isUndefined .eqv. .false.
    end function test_constructBlank_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructPosition_1() result(assertion)
        implicit none
        logical :: assertion
        type(Position_type) :: Position
        Position = Position_type()
        assertion = .not. Position%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Position%value == "asis"
        assertion = assertion .and. Position%isRewind .eqv. .false.
        assertion = assertion .and. Position%isAppend .eqv. .false.
        assertion = assertion .and. Position%isAsis .eqv. .true.
        assertion = assertion .and. Position%isUndefined .eqv. .false.
    end function test_constructPosition_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructPosition_2() result(assertion)
        implicit none
        logical :: assertion
        type(Position_type) :: Position
        Position = Position_type("undefined")
        assertion = .not. Position%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Position%value == "undefined"
        assertion = assertion .and. Position%isRewind .eqv. .false.
        assertion = assertion .and. Position%isAppend .eqv. .false.
        assertion = assertion .and. Position%isAsis .eqv. .false.
        assertion = assertion .and. Position%isUndefined .eqv. .true.
    end function test_constructPosition_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructPosition_3() result(assertion)
        implicit none
        logical :: assertion
        type(Position_type) :: Position
        Position = Position_type("rewind")
        assertion = .not. Position%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Position%value == "rewind"
        assertion = assertion .and. Position%isRewind .eqv. .true.
        assertion = assertion .and. Position%isAppend .eqv. .false.
        assertion = assertion .and. Position%isAsis .eqv. .false.
        assertion = assertion .and. Position%isUndefined .eqv. .false.
    end function test_constructPosition_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructPosition_4() result(assertion)
        implicit none
        logical :: assertion
        type(Position_type) :: Position
        Position = Position_type("APPEND")
        assertion = .not. Position%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Position%value == "append"
        assertion = assertion .and. Position%isRewind .eqv. .false.
        assertion = assertion .and. Position%isAppend .eqv. .true.
        assertion = assertion .and. Position%isAsis .eqv. .false.
        assertion = assertion .and. Position%isUndefined .eqv. .false.
    end function test_constructPosition_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructPosition_5() result(assertion)
        implicit none
        logical :: assertion
        type(Position_type) :: Position
        Position = Position_type("ASIS")
        assertion = .not. Position%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Position%value == "asis"
        assertion = assertion .and. Position%isRewind .eqv. .false.
        assertion = assertion .and. Position%isAppend .eqv. .false.
        assertion = assertion .and. Position%isAsis .eqv. .true.
        assertion = assertion .and. Position%isUndefined .eqv. .false.
    end function test_constructPosition_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructPosition_6() result(assertion)
        implicit none
        logical :: assertion
        type(Position_type) :: Position
        Position = Position_type("nonsense")
        assertion = Position%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Position%value == ""
        assertion = assertion .and. Position%isRewind .eqv. .false.
        assertion = assertion .and. Position%isAppend .eqv. .false.
        assertion = assertion .and. Position%isAsis .eqv. .false.
        assertion = assertion .and. Position%isUndefined .eqv. .false.
    end function test_constructPosition_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructDelim_1() result(assertion)
        implicit none
        logical :: assertion
        type(Delim_type) :: Delim
        Delim = Delim_type()
        assertion = .not. Delim%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Delim%value == "none"
        assertion = assertion .and. Delim%isQuote .eqv. .false.
        assertion = assertion .and. Delim%isApostrophe .eqv. .false.
        assertion = assertion .and. Delim%isNone .eqv. .true.
        assertion = assertion .and. Delim%isUndefined .eqv. .false.
    end function test_constructDelim_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructDelim_2() result(assertion)
        implicit none
        logical :: assertion
        type(Delim_type) :: Delim
        Delim = Delim_type("Undefined")
        assertion = .not. Delim%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Delim%value == "undefined"
        assertion = assertion .and. Delim%isQuote .eqv. .false.
        assertion = assertion .and. Delim%isApostrophe .eqv. .false.
        assertion = assertion .and. Delim%isNone .eqv. .false.
        assertion = assertion .and. Delim%isUndefined .eqv. .true.
    end function test_constructDelim_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructDelim_3() result(assertion)
        implicit none
        logical :: assertion
        type(Delim_type) :: Delim
        Delim = Delim_type("Quote")
        assertion = .not. Delim%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Delim%value == "quote"
        assertion = assertion .and. Delim%isQuote .eqv. .true.
        assertion = assertion .and. Delim%isApostrophe .eqv. .false.
        assertion = assertion .and. Delim%isNone .eqv. .false.
        assertion = assertion .and. Delim%isUndefined .eqv. .false.
    end function test_constructDelim_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructDelim_4() result(assertion)
        implicit none
        logical :: assertion
        type(Delim_type) :: Delim
        Delim = Delim_type("Apostrophe")
        assertion = .not. Delim%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Delim%value == "apostrophe"
        assertion = assertion .and. Delim%isQuote .eqv. .false.
        assertion = assertion .and. Delim%isApostrophe .eqv. .true.
        assertion = assertion .and. Delim%isNone .eqv. .false.
        assertion = assertion .and. Delim%isUndefined .eqv. .false.
    end function test_constructDelim_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructDelim_5() result(assertion)
        implicit none
        logical :: assertion
        type(Delim_type) :: Delim
        Delim = Delim_type("None")
        assertion = .not. Delim%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Delim%value == "none"
        assertion = assertion .and. Delim%isQuote .eqv. .false.
        assertion = assertion .and. Delim%isApostrophe .eqv. .false.
        assertion = assertion .and. Delim%isNone .eqv. .true.
        assertion = assertion .and. Delim%isUndefined .eqv. .false.
    end function test_constructDelim_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructDelim_6() result(assertion)
        implicit none
        logical :: assertion
        type(Delim_type) :: Delim
        Delim = Delim_type("nonsense")
        assertion = Delim%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Delim%value == ""
        assertion = assertion .and. Delim%isQuote .eqv. .false.
        assertion = assertion .and. Delim%isApostrophe .eqv. .false.
        assertion = assertion .and. Delim%isNone .eqv. .false.
        assertion = assertion .and. Delim%isUndefined .eqv. .false.
    end function test_constructDelim_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructRound_1() result(assertion)
        implicit none
        logical :: assertion
        type(Round_type) :: Round
        Round = Round_type()
        assertion = .not. Round%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Round%value == "processor_defined"
        assertion = assertion .and. Round%isUp .eqv. .false.
        assertion = assertion .and. Round%isDown .eqv. .false.
        assertion = assertion .and. Round%isZero .eqv. .false.
        assertion = assertion .and. Round%isNearest .eqv. .false.
        assertion = assertion .and. Round%isCompatible .eqv. .false.
        assertion = assertion .and. Round%isProcessDefined .eqv. .true.
        assertion = assertion .and. Round%isUndefined .eqv. .false.
    end function test_constructRound_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructRound_2() result(assertion)
        implicit none
        logical :: assertion
        type(Round_type) :: Round
        Round = Round_type("UP")
        assertion = .not. Round%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Round%value == "up"
        assertion = assertion .and. Round%isUp .eqv. .true.
        assertion = assertion .and. Round%isDown .eqv. .false.
        assertion = assertion .and. Round%isZero .eqv. .false.
        assertion = assertion .and. Round%isNearest .eqv. .false.
        assertion = assertion .and. Round%isCompatible .eqv. .false.
        assertion = assertion .and. Round%isProcessDefined .eqv. .false.
        assertion = assertion .and. Round%isUndefined .eqv. .false.
    end function test_constructRound_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructRound_3() result(assertion)
        implicit none
        logical :: assertion
        type(Round_type) :: Round
        Round = Round_type("Down")
        assertion = .not. Round%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Round%value == "down"
        assertion = assertion .and. Round%isUp .eqv. .false.
        assertion = assertion .and. Round%isDown .eqv. .true.
        assertion = assertion .and. Round%isZero .eqv. .false.
        assertion = assertion .and. Round%isNearest .eqv. .false.
        assertion = assertion .and. Round%isCompatible .eqv. .false.
        assertion = assertion .and. Round%isProcessDefined .eqv. .false.
        assertion = assertion .and. Round%isUndefined .eqv. .false.
    end function test_constructRound_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructRound_4() result(assertion)
        implicit none
        logical :: assertion
        type(Round_type) :: Round
        Round = Round_type("Zero")
        assertion = .not. Round%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Round%value == "zero"
        assertion = assertion .and. Round%isUp .eqv. .false.
        assertion = assertion .and. Round%isDown .eqv. .false.
        assertion = assertion .and. Round%isZero .eqv. .true.
        assertion = assertion .and. Round%isNearest .eqv. .false.
        assertion = assertion .and. Round%isCompatible .eqv. .false.
        assertion = assertion .and. Round%isProcessDefined .eqv. .false.
        assertion = assertion .and. Round%isUndefined .eqv. .false.
    end function test_constructRound_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructRound_5() result(assertion)
        implicit none
        logical :: assertion
        type(Round_type) :: Round
        Round = Round_type("Nearest")
        assertion = .not. Round%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Round%value == "nearest"
        assertion = assertion .and. Round%isUp .eqv. .false.
        assertion = assertion .and. Round%isDown .eqv. .false.
        assertion = assertion .and. Round%isZero .eqv. .false.
        assertion = assertion .and. Round%isNearest .eqv. .true.
        assertion = assertion .and. Round%isCompatible .eqv. .false.
        assertion = assertion .and. Round%isProcessDefined .eqv. .false.
        assertion = assertion .and. Round%isUndefined .eqv. .false.
    end function test_constructRound_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructRound_6() result(assertion)
        implicit none
        logical :: assertion
        type(Round_type) :: Round
        Round = Round_type("Nearest")
        assertion = .not. Round%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Round%value == "nearest"
        assertion = assertion .and. Round%isUp .eqv. .false.
        assertion = assertion .and. Round%isDown .eqv. .false.
        assertion = assertion .and. Round%isZero .eqv. .false.
        assertion = assertion .and. Round%isNearest .eqv. .true.
        assertion = assertion .and. Round%isCompatible .eqv. .false.
        assertion = assertion .and. Round%isProcessDefined .eqv. .false.
        assertion = assertion .and. Round%isUndefined .eqv. .false.
    end function test_constructRound_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructRound_7() result(assertion)
        implicit none
        logical :: assertion
        type(Round_type) :: Round
        Round = Round_type("Compatible")
        assertion = .not. Round%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Round%value == "compatible"
        assertion = assertion .and. Round%isUp .eqv. .false.
        assertion = assertion .and. Round%isDown .eqv. .false.
        assertion = assertion .and. Round%isZero .eqv. .false.
        assertion = assertion .and. Round%isNearest .eqv. .false.
        assertion = assertion .and. Round%isCompatible .eqv. .true.
        assertion = assertion .and. Round%isProcessDefined .eqv. .false.
        assertion = assertion .and. Round%isUndefined .eqv. .false.
    end function test_constructRound_7

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructRound_8() result(assertion)
        implicit none
        logical :: assertion
        type(Round_type) :: Round
        Round = Round_type("Processor_defined")
        assertion = .not. Round%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Round%value == "processor_defined"
        assertion = assertion .and. Round%isUp .eqv. .false.
        assertion = assertion .and. Round%isDown .eqv. .false.
        assertion = assertion .and. Round%isZero .eqv. .false.
        assertion = assertion .and. Round%isNearest .eqv. .false.
        assertion = assertion .and. Round%isCompatible .eqv. .false.
        assertion = assertion .and. Round%isProcessDefined .eqv. .true.
        assertion = assertion .and. Round%isUndefined .eqv. .false.
    end function test_constructRound_8

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructRound_9() result(assertion)
        implicit none
        logical :: assertion
        type(Round_type) :: Round
        Round = Round_type("nonsense")
        assertion = Round%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Round%value == ""
        assertion = assertion .and. Round%isUp .eqv. .false.
        assertion = assertion .and. Round%isDown .eqv. .false.
        assertion = assertion .and. Round%isZero .eqv. .false.
        assertion = assertion .and. Round%isNearest .eqv. .false.
        assertion = assertion .and. Round%isCompatible .eqv. .false.
        assertion = assertion .and. Round%isProcessDefined .eqv. .false.
        assertion = assertion .and. Round%isUndefined .eqv. .false.
    end function test_constructRound_9

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructSign_1() result(assertion)
        implicit none
        logical :: assertion
        type(Sign_type) :: Sign
        Sign = Sign_type()
        assertion = .not. Sign%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Sign%value == "processor_defined"
        assertion = assertion .and. Sign%isSuppress .eqv. .false.
        assertion = assertion .and. Sign%isPlus .eqv. .false.
        assertion = assertion .and. Sign%isProcessDefined .eqv. .true.
        assertion = assertion .and. Sign%isUndefined .eqv. .false.
    end function test_constructSign_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructSign_2() result(assertion)
        implicit none
        logical :: assertion
        type(Sign_type) :: Sign
        Sign = Sign_type("Undefined")
        assertion = .not. Sign%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Sign%value == "undefined"
        assertion = assertion .and. Sign%isSuppress .eqv. .false.
        assertion = assertion .and. Sign%isPlus .eqv. .false.
        assertion = assertion .and. Sign%isProcessDefined .eqv. .false.
        assertion = assertion .and. Sign%isUndefined .eqv. .true.
    end function test_constructSign_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructSign_3() result(assertion)
        implicit none
        logical :: assertion
        type(Sign_type) :: Sign
        Sign = Sign_type("Suppress")
        assertion = .not. Sign%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Sign%value == "suppress"
        assertion = assertion .and. Sign%isSuppress .eqv. .true.
        assertion = assertion .and. Sign%isPlus .eqv. .false.
        assertion = assertion .and. Sign%isProcessDefined .eqv. .false.
        assertion = assertion .and. Sign%isUndefined .eqv. .false.
    end function test_constructSign_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructSign_4() result(assertion)
        implicit none
        logical :: assertion
        type(Sign_type) :: Sign
        Sign = Sign_type("Plus")
        assertion = .not. Sign%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Sign%value == "plus"
        assertion = assertion .and. Sign%isSuppress .eqv. .false.
        assertion = assertion .and. Sign%isPlus .eqv. .true.
        assertion = assertion .and. Sign%isProcessDefined .eqv. .false.
        assertion = assertion .and. Sign%isUndefined .eqv. .false.
    end function test_constructSign_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructSign_5() result(assertion)
        implicit none
        logical :: assertion
        type(Sign_type) :: Sign
        Sign = Sign_type("Processor_defined")
        assertion = .not. Sign%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Sign%value == "processor_defined"
        assertion = assertion .and. Sign%isSuppress .eqv. .false.
        assertion = assertion .and. Sign%isPlus .eqv. .false.
        assertion = assertion .and. Sign%isProcessDefined .eqv. .true.
        assertion = assertion .and. Sign%isUndefined .eqv. .false.
    end function test_constructSign_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_constructSign_6() result(assertion)
        implicit none
        logical :: assertion
        type(Sign_type) :: Sign
        Sign = Sign_type("nonsense")
        assertion = Sign%Err%occurred; if (.not. assertion) return
        assertion = assertion .and. Sign%Err%occurred .eqv. .true.
        assertion = assertion .and. Sign%value == ""
        assertion = assertion .and. Sign%isSuppress .eqv. .false.
        assertion = assertion .and. Sign%isPlus .eqv. .false.
        assertion = assertion .and. Sign%isProcessDefined .eqv. .false.
        assertion = assertion .and. Sign%isUndefined .eqv. .false.
    end function test_constructSign_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getWriteErr_1() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        integer :: stat
        stat = 1
        Err = getWriteErr(stat) ! LCOV_EXCL_LINE
        assertion = Err%occurred
    end function test_getWriteErr_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getWriteErr_2() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        integer :: stat
        stat = 0
        Err = getWriteErr(stat) ! LCOV_EXCL_LINE
        assertion = .not. Err%occurred
    end function test_getWriteErr_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getWriteErr_3() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        integer :: stat
        stat = -1
        Err = getWriteErr(stat) ! LCOV_EXCL_LINE
        assertion = Err%occurred
    end function test_getWriteErr_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getReadErr_1() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        integer :: stat
        stat = 1
        Err = getReadErr(stat,"./path") ! LCOV_EXCL_LINE
        assertion = Err%occurred
    end function test_getReadErr_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getReadErr_2() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        integer :: stat
        stat = 1
        Err = getReadErr(stat) ! LCOV_EXCL_LINE
        assertion = Err%occurred
    end function test_getReadErr_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getReadErr_3() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        integer :: stat
        stat = 0
        Err = getReadErr(stat,"./path") ! LCOV_EXCL_LINE
        assertion = .not. Err%occurred
    end function test_getReadErr_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getReadErr_4() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        integer :: stat
        stat = 0
        Err = getReadErr(stat) ! LCOV_EXCL_LINE
        assertion = .not. Err%occurred
    end function test_getReadErr_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getReadErr_5() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        integer :: stat
        stat = -1
        Err = getReadErr(stat,"./path") ! LCOV_EXCL_LINE
        assertion = Err%occurred
    end function test_getReadErr_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getReadErr_6() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        integer :: stat
        stat = -1
        Err = getReadErr(stat) ! LCOV_EXCL_LINE
        assertion = Err%occurred
    end function test_getReadErr_6

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getCloseErr_1() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        integer :: stat
        stat = 0
        Err = getCloseErr(stat) ! LCOV_EXCL_LINE
        assertion = .not. Err%occurred
    end function test_getCloseErr_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getCloseErr_2() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        integer :: stat
        stat = 1
        Err = getCloseErr(stat) ! LCOV_EXCL_LINE
        assertion = Err%occurred
    end function test_getCloseErr_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getCloseErr_3() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        Err = getCloseErr(-1) ! LCOV_EXCL_LINE
        assertion = .not. Err%occurred
    end function test_getCloseErr_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getOpenErr_1() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        Err = getOpenErr(0) ! LCOV_EXCL_LINE
        assertion = .not. Err%occurred
    end function test_getOpenErr_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getOpenErr_2() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        Err = getOpenErr(1) ! LCOV_EXCL_LINE
        assertion = Err%occurred
    end function test_getOpenErr_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getOpenErr_3() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        Err = getOpenErr(-1) ! LCOV_EXCL_LINE
        assertion = .not. Err%occurred
    end function test_getOpenErr_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getInqErr_1() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        Err = getInqErr(0) ! LCOV_EXCL_LINE
        assertion = .not. Err%occurred
    end function test_getInqErr_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getInqErr_2() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        Err = getInqErr(1) ! LCOV_EXCL_LINE
        assertion = Err%occurred
    end function test_getInqErr_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getInqErr_3() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        Err = getInqErr(-1) ! LCOV_EXCL_LINE
        assertion = Err%occurred
    end function test_getInqErr_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getExistStatus_1() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        logical :: exists
        type(Err_type) :: Err
        call getExistStatus(exists, Err, unit = -1, file = "")
        assertion = Err%occurred
    end function test_getExistStatus_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getExistStatus_2() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        logical :: exists
        type(Err_type) :: Err
        call getExistStatus(exists, Err)
        assertion = Err%occurred
    end function test_getExistStatus_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getExistStatus_3() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        logical :: exists
        type(Err_type) :: Err
        call getExistStatus(exists, Err, unit = 123)
        assertion = .not. Err%occurred
    end function test_getExistStatus_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getOpenStatus_1() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        logical :: isOpen
        type(Err_type) :: Err
        call getOpenStatus(isOpen, Err, unit = 124, file = "")
        assertion = Err%occurred
    end function test_getOpenStatus_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getOpenStatus_2() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        logical :: isOpen
        type(Err_type) :: Err
        call getOpenStatus(isOpen, Err)
        assertion = Err%occurred
    end function test_getOpenStatus_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getOpenStatus_3() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        logical :: isOpen
        type(Err_type) :: Err
        call getOpenStatus(isOpen, Err, unit = 123)
        assertion = .not. Err%occurred
    end function test_getOpenStatus_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Check inquire with an input file name.
    function test_getOpenStatus_4() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        logical :: isOpen
        type(Err_type) :: Err
        call getOpenStatus(isOpen, Err, file = "nonexisting.file")
        assertion = .not. isOpen .and. .not. Err%occurred
    end function test_getOpenStatus_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getNumber_1() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        logical :: isNumbered
        integer :: number
        call getNumber(isNumbered, number, Err, unit = 124, file = "")
        assertion = Err%occurred
    end function test_getNumber_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getNumber_2() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        logical :: isNumbered
        integer :: number
        call getNumber(isNumbered, number, Err)
        assertion = Err%occurred
    end function test_getNumber_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getNumber_3() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        logical :: isNumbered
        integer :: number
        call getNumber(isNumbered, number, Err, unit = 124)
        assertion = .not. Err%occurred .and. number == -1 .and. .not. isNumbered
    end function test_getNumber_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Check inquire with an input file name.
    function test_getNumber_4() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        logical :: isNumbered
        integer :: number
        call getNumber(isNumbered, number, Err, file = "nonexisting.file")
        assertion = .not. Err%occurred .and. number == -1 .and. .not. isNumbered
    end function test_getNumber_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getName_1() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        logical :: isNamed
        character(:), allocatable :: nameByCompiler
        call getName(isNamed, nameByCompiler, Err, unit = 124, file = "")
        assertion = Err%occurred
    end function test_getName_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getName_2() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        logical :: isNamed
        character(:), allocatable :: nameByCompiler
        call getName(isNamed, nameByCompiler, Err)
        assertion = Err%occurred
    end function test_getName_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getName_3() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        logical :: isNamed
        character(:), allocatable :: nameByCompiler
        call getName(isNamed, nameByCompiler, Err, unit = 124)
        assertion = .not. Err%occurred .and. .not. isNamed
    end function test_getName_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Check inquire with an input file name.
    function test_getName_4() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        logical :: isNamed
        character(:), allocatable :: nameByCompiler
        call getName(isNamed, nameByCompiler, Err, file = "nonexisting.file")
        assertion = .not. Err%occurred .and. isNamed
#if !defined CODECOV_ENABLED
        assertion = .true.
#endif
    end function test_getName_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getAccess_1() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: access
        call getAccess(access, Err, unit = 124, file = "")
        assertion = Err%occurred
    end function test_getAccess_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getAccess_2() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: access
        call getAccess(access, Err)
        assertion = Err%occurred
    end function test_getAccess_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getAccess_3() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: access
        call getAccess(access, Err, unit = 124)
        assertion = .not. Err%occurred
    end function test_getAccess_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Check inquire with an input file name.
    function test_getAccess_4() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: access
        call getAccess(access, Err, file = "nonexisting.file")
        assertion = .not. Err%occurred
    end function test_getAccess_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getForm_1() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: form
        call getForm(form, Err, unit = 124, file = "")
        assertion = Err%occurred
    end function test_getForm_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getForm_2() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: form
        call getForm(form, Err)
        assertion = Err%occurred
    end function test_getForm_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getForm_3() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: form
        call getForm(form, Err, unit = 124)
        assertion = .not. Err%occurred
    end function test_getForm_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Check inquire with an input file name.
    function test_getForm_4() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: form
        call getForm(form, Err, file = "nonexisting.file")
        assertion = .not. Err%occurred
    end function test_getForm_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getRecl_1() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        integer :: recl
        call getRecl(recl, Err, unit = 124, file = "")
        assertion = Err%occurred
    end function test_getRecl_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getRecl_2() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        integer :: recl
        call getRecl(recl, Err)
        assertion = Err%occurred
    end function test_getRecl_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getRecl_3() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        integer :: recl
        call getRecl(recl, Err, unit = 124)
        assertion = .not. Err%occurred
    end function test_getRecl_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Check inquire with an input file name.
    function test_getRecl_4() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        integer :: recl
        call getRecl(recl, Err, file = "nonexisting.file")
        assertion = .not. Err%occurred
    end function test_getRecl_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getBlank_1() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: blank
        call getBlank(blank, Err, unit = 124, file = "")
        assertion = Err%occurred
    end function test_getBlank_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getBlank_2() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: blank
        call getBlank(blank, Err)
        assertion = Err%occurred
    end function test_getBlank_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getBlank_3() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: blank
        call getBlank(blank, Err, unit = 124)
        assertion = .not. Err%occurred .and. blank == "undefined"
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "Err%occurred : ", Err%occurred
            write(Test%outputUnit,"(*(g0))")   "blank : ", blank
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getBlank_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Check inquire with an input file name.
    function test_getBlank_4() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: blank
        call getBlank(blank, Err, file = "nonexisting.file")
        assertion = .not. Err%occurred .and. blank == "undefined"
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "Err%occurred : ", Err%occurred
            write(Test%outputUnit,"(*(g0))")   "blank : ", blank
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getBlank_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getPosition_1() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: position
        call getPosition(position, Err, unit = 124, file = "")
        assertion = Err%occurred
    end function test_getPosition_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getPosition_2() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: position
        call getPosition(position, Err)
        assertion = Err%occurred
    end function test_getPosition_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getPosition_3() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: position
        call getPosition(position, Err, unit = 124)
        assertion = .not. Err%occurred .and. position == "undefined"
        if (Test%isDebugMode .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "Err%occurred : ", Err%occurred
            write(Test%outputUnit,"(*(g0))")   "position : ", position
            write(Test%outputUnit,"(*(g0))")
            ! LCOV_EXCL_STOP
        end if
    end function test_getPosition_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Check inquire with an input file name.
    function test_getPosition_4() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: position
        call getPosition(position, Err, file = "nonexisting.file")
        assertion = .not. Err%occurred .and. position == "undefined"
        if (Test%isDebugMode .and. .not. assertion) then
            ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "Err%occurred : ", Err%occurred
            write(Test%outputUnit,"(*(g0))")   "position : ", position
            write(Test%outputUnit,"(*(g0))")
            ! LCOV_EXCL_STOP
        end if
    end function test_getPosition_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getAction_1() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: action
        call getAction(action, Err, unit = 124, file = "")
        assertion = Err%occurred
    end function test_getAction_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getAction_2() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: action
        call getAction(action, Err)
        assertion = Err%occurred
    end function test_getAction_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getAction_3() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: action
        call getAction(action, Err, unit = 124)
        assertion = .not. Err%occurred .and. action == "undefined"
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "Err%occurred : ", Err%occurred
            write(Test%outputUnit,"(*(g0))")   "action : ", action
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getAction_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Check inquire with an input file name.
    function test_getAction_4() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: action
        call getAction(action, Err, file = "nonexisting.file")
        assertion = .not. Err%occurred !.and. action == "undefined"
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "Err%occurred : ", Err%occurred
            write(Test%outputUnit,"(*(g0))")   "action : ", action
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getAction_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input arguments `unit` and `file` must NOT be present simultaneously.
    function test_getDelim_1() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: delim
        call getDelim(delim, Err, unit = 124, file = "")
        assertion = Err%occurred
    end function test_getDelim_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> At least one of the two identifiers (`unit` or `file`) must be present.
    function test_getDelim_2() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: delim
        call getDelim(delim, Err)
        assertion = Err%occurred
    end function test_getDelim_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The input `unit` can point to any existing or non-existing opened or closed file at runtime.
    function test_getDelim_3() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: delim
        call getDelim(delim, Err, unit = 124)
        assertion = .not. Err%occurred .and. delim == "undefined"
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "Err%occurred : ", Err%occurred
            write(Test%outputUnit,"(*(g0))")   "delim : ", delim
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getDelim_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Check inquire with an input file name.
    function test_getDelim_4() result(assertion)
        use Err_mod, only: Err_type
        implicit none
        logical :: assertion
        type(Err_type) :: Err
        character(:), allocatable :: delim
        call getDelim(delim, Err, file = "nonexisting.file")
        assertion = .not. Err%occurred !.and. delim == "undefined"
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "Err%occurred : ", Err%occurred
            write(Test%outputUnit,"(*(g0))")   "delim : ", delim
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_getDelim_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_File_mod ! LCOV_EXCL_LINE