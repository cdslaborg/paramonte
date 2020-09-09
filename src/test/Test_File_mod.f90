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

module Test_File_mod

    use File_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_File

    type(Test_type) :: Test

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_File()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)

        call test_File_type()
        call test_Pad_type()
        call test_Action_type()
        call test_Access_type()
        call test_Form_type()
        call test_Blank_type()
        call test_Position_type()
        call test_Delim_type()
        call test_Round_type()
        call test_Sign_type()
        call Test%finalize()
        call test_getFileExistStatus()
        call test_getFileOpenStatus()
        call test_getFileNumber()
        call test_getFileName()
        call test_getFileAccess()
        call test_getFileForm()
        call test_getFileRecl()
        call test_getFileBlank()
        call test_getFilePosition()
        call test_getFileAction()
        call test_getFileDelim()
        call test_getWriteErr()
        call test_getReadErr()
        call test_getCloseErr()
        call test_getOpenErr()
        call test_getInqErr()
#ifdef CAF_ENABLED
        sync all
#endif


    end subroutine test_File

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_File_type()
        
        !use Constants_mod, only: RK
        implicit none
        type(File_type) :: File


        if (Test%Image%isFirst) call Test%testing("File_type default values")

        File = File_type()
        call Test%checkForErr(File%Err)

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "File%Path%original  : ", File%Path%original
            write(Test%outputUnit,"(*(g0))")   "File%Path%modified  : ", File%Path%modified
            write(Test%outputUnit,"(*(g0))")   "File%Path%dir       : ", File%Path%dir
            write(Test%outputUnit,"(*(g0))")   "File%Path%name      : ", File%Path%name
            write(Test%outputUnit,"(*(g0))")   "File%Path%ext       : ", File%Path%ext
            write(Test%outputUnit,"(*(g0))")   "File%Path%slashOS   : ", File%Path%slashOS
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
                            Test%assertion = File%Path%original  == ""
        call Test%verify(); Test%assertion = File%Path%modified  == ""
        call Test%verify(); Test%assertion = File%Path%dir       == ""
        call Test%verify(); Test%assertion = File%Path%name      == ""
        call Test%verify(); Test%assertion = File%Path%ext       == ""
        call Test%verify(); Test%assertion = File%unit           == -2147483647
        call Test%verify(); Test%assertion = File%number         == -2147483647
        call Test%verify(); Test%assertion = File%recl           == -2147483647
        call Test%verify(); Test%assertion = File%exists         .eqv. .false.
        call Test%verify(); Test%assertion = File%isOpen         .eqv. .false.
        call Test%verify(); Test%assertion = File%isNamed        .eqv. .false.
        call Test%verify(); Test%assertion = File%isNumbered     .eqv. .false.
        call Test%verify(); Test%assertion = File%status         == "unknown"
        call Test%verify(); Test%assertion = File%asynchronous   == "no"
        call Test%verify(); Test%assertion = File%format         == ""
        call Test%verify(); Test%assertion = File%nameByCompiler == ""
        call Test%verify(); Test%assertion = File%Action%value   == "readwrite"
        call Test%verify(); Test%assertion = File%Access%value   == "sequential"
        call Test%verify(); Test%assertion = File%Form%value     == "formatted"
        call Test%verify(); Test%assertion = File%Blank%value    == "null"
        call Test%verify(); Test%assertion = File%Position%value == "asis"
        call Test%verify(); Test%assertion = File%Delim%value    == "none"
        call Test%verify(); Test%assertion = File%Pad%value      == "yes"
        call Test%verify(); Test%assertion = File%Round%value    == "processor_defined"
        call Test%verify(); Test%assertion = File%Sign%value     == "processor_defined"
        call Test%verify(); Test%assertion = File%Err%occurred   .eqv. .false.
        call Test%verify(); Test%assertion = File%Err%msg        == ""
        call Test%verify()




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
        call Test%checkForErr(File%Err)

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))")   "File%Path%original  : ", File%Path%original
            write(Test%outputUnit,"(*(g0))")   "File%Path%modified  : ", File%Path%modified
            write(Test%outputUnit,"(*(g0))")   "File%Path%dir       : ", File%Path%dir
            write(Test%outputUnit,"(*(g0))")   "File%Path%name      : ", File%Path%name
            write(Test%outputUnit,"(*(g0))")   "File%Path%ext       : ", File%Path%ext
            write(Test%outputUnit,"(*(g0))")   "File%Path%slashOS   : ", File%Path%slashOS
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
        Test%assertion = File%Path%original  == "./test_File_mod/\test_File_mod\-"
        if (File%Path%slashOS=="/") then
            call Test%verify(); Test%assertion = File%Path%modified  == "./test_File_mod/\test_File_mod\-"
            call Test%verify(); Test%assertion = File%Path%dir       == "./test_File_mod/"
            call Test%verify(); Test%assertion = File%Path%name      == "\test_File_mod\-"
            call Test%verify(); Test%assertion = File%Path%ext       == ""
        else
            call Test%verify(); Test%assertion = File%Path%modified  == ".\test_File_mod\\test_File_mod\-"
            call Test%verify(); Test%assertion = File%Path%dir       == ".\test_File_mod\\test_File_mod\"
            call Test%verify(); Test%assertion = File%Path%name      == "-"
            call Test%verify(); Test%assertion = File%Path%ext       == ""
        end if
        call Test%verify(); Test%assertion = File%unit           == 13
        call Test%verify(); Test%assertion = File%number         == -2147483647
        call Test%verify(); Test%assertion = File%recl           == 9999
        call Test%verify(); Test%assertion = File%exists         .eqv. .false.
        call Test%verify(); Test%assertion = File%isOpen         .eqv. .false.
        call Test%verify(); Test%assertion = File%isNamed        .eqv. .false.
        call Test%verify(); Test%assertion = File%isNumbered     .eqv. .false.
        call Test%verify(); Test%assertion = File%status         == "new"
        call Test%verify(); Test%assertion = File%asynchronous   == "yes"
        call Test%verify(); Test%assertion = File%format         == "(A)"
        call Test%verify(); Test%assertion = File%nameByCompiler == ""
        call Test%verify(); Test%assertion = File%Action%value   == "read"
        call Test%verify(); Test%assertion = File%Access%value   == "direct"
        call Test%verify(); Test%assertion = File%Form%value     == "unformatted"
        call Test%verify(); Test%assertion = File%Blank%value    == "undefined"
        call Test%verify(); Test%assertion = File%Position%value == "append"
        call Test%verify(); Test%assertion = File%Delim%value    == "quote"
        call Test%verify(); Test%assertion = File%Pad%value      == "no"
        call Test%verify(); Test%assertion = File%Round%value    == "up"
        call Test%verify(); Test%assertion = File%Sign%value     == "undefined"
        call Test%verify(); Test%assertion = File%Err%occurred   .eqv. .false.
        call Test%verify(); Test%assertion = File%Err%msg        == ""
        call Test%verify()

    end subroutine test_File_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_Pad_type

        type(Pad_type) :: Pad

        if (Test%Image%isFirst) call Test%testing("Pad_type()")
        Pad = Pad_type()
        call Test%checkForErr(Pad%Err)
        Test%assertion = Pad%value == "yes"; call Test%verify()
        Test%assertion = Pad%isPadded .eqv. .true.; call Test%verify()
        Test%assertion = Pad%isNotPadded .eqv. .false.; call Test%verify()
        Test%assertion = Pad%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Pad_type(value='undefined')")
        Pad = Pad_type(value='undefined')
        call Test%checkForErr(Pad%Err)
        Test%assertion = Pad%value == "undefined"; call Test%verify()
        Test%assertion = Pad%isPadded .eqv. .false.; call Test%verify()
        Test%assertion = Pad%isNotPadded .eqv. .false.; call Test%verify()
        Test%assertion = Pad%isUndefined .eqv. .true.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Pad_type(value='yes')")
        Pad = Pad_type(value='yes')
        call Test%checkForErr(Pad%Err)
        Test%assertion = Pad%value == "yes"; call Test%verify()
        Test%assertion = Pad%isPadded .eqv. .true.; call Test%verify()
        Test%assertion = Pad%isNotPadded .eqv. .false.; call Test%verify()
        Test%assertion = Pad%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Pad_type(value='no')")
        Pad = Pad_type(value='no')
        call Test%checkForErr(Pad%Err)
        Test%assertion = Pad%value == "no"; call Test%verify()
        Test%assertion = Pad%isPadded .eqv. .false.; call Test%verify()
        Test%assertion = Pad%isNotPadded .eqv. .true.; call Test%verify()
        Test%assertion = Pad%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Pad_type(value='nonsense')")
        Pad = Pad_type(value='nonsense')
        Test%assertion = Pad%Err%occurred .eqv. .true.; call Test%verify()
        Test%assertion = Pad%value == ""; call Test%verify()
        Test%assertion = Pad%isPadded .eqv. .false.; call Test%verify()
        Test%assertion = Pad%isNotPadded .eqv. .false.; call Test%verify()
        Test%assertion = Pad%isUndefined .eqv. .false.; call Test%verify()

    end subroutine test_Pad_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_Action_type

        type(Action_type) :: Action

        if (Test%Image%isFirst) call Test%testing("Action_type()")
        Action = Action_type()
        call Test%checkForErr(Action%Err)
        Test%assertion = Action%value == "readwrite"; call Test%verify()
        Test%assertion = Action%isRead .eqv. .false.; call Test%verify()
        Test%assertion = Action%isWrite .eqv. .false.; call Test%verify()
        Test%assertion = Action%isReadWrite .eqv. .true.; call Test%verify()
        Test%assertion = Action%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Action_type(value='undefined')")
        Action = Action_type(value='undefined')
        call Test%checkForErr(Action%Err)
        Test%assertion = Action%value == "undefined"; call Test%verify()
        Test%assertion = Action%isRead .eqv. .false.; call Test%verify()
        Test%assertion = Action%isWrite .eqv. .false.; call Test%verify()
        Test%assertion = Action%isReadWrite .eqv. .false.; call Test%verify()
        Test%assertion = Action%isUndefined .eqv. .true.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Action_type(value='read')")
        Action = Action_type(value='read')
        call Test%checkForErr(Action%Err)
        Test%assertion = Action%value == "read"; call Test%verify()
        Test%assertion = Action%isRead .eqv. .true.; call Test%verify()
        Test%assertion = Action%isWrite .eqv. .false.; call Test%verify()
        Test%assertion = Action%isReadWrite .eqv. .false.; call Test%verify()
        Test%assertion = Action%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Action_type(value='write')")
        Action = Action_type(value='write')
        call Test%checkForErr(Action%Err)
        Test%assertion = Action%value == "write"; call Test%verify()
        Test%assertion = Action%isRead .eqv. .false.; call Test%verify()
        Test%assertion = Action%isWrite .eqv. .true.; call Test%verify()
        Test%assertion = Action%isReadWrite .eqv. .false.; call Test%verify()
        Test%assertion = Action%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Action_type(value='readwrite')")
        Action = Action_type(value='readwrite')
        call Test%checkForErr(Action%Err)
        Test%assertion = Action%value == "readwrite"; call Test%verify()
        Test%assertion = Action%isRead .eqv. .false.; call Test%verify()
        Test%assertion = Action%isWrite .eqv. .false.; call Test%verify()
        Test%assertion = Action%isReadWrite .eqv. .true.; call Test%verify()
        Test%assertion = Action%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Action_type(value='nonsense')")
        Action = Action_type(value='nonsense')
        Test%assertion = Action%Err%occurred .eqv. .true.; call Test%verify()
        Test%assertion = Action%value == ""; call Test%verify()
        Test%assertion = Action%isRead .eqv. .false.; call Test%verify()
        Test%assertion = Action%isWrite .eqv. .false.; call Test%verify()
        Test%assertion = Action%isReadWrite .eqv. .false.; call Test%verify()
        Test%assertion = Action%isUndefined .eqv. .false.; call Test%verify()

    end subroutine test_Action_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_Access_type

        type(Access_type) :: Access

        if (Test%Image%isFirst) call Test%testing("Access_type()")
        Access = Access_type()
        call Test%checkForErr(Access%Err)
        Test%assertion = Access%value == "sequential"; call Test%verify()
        Test%assertion = Access%isSequential .eqv. .true.; call Test%verify()
        Test%assertion = Access%isDirect .eqv. .false.; call Test%verify()
        Test%assertion = Access%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Access_type(value='undefined')")
        Access = Access_type(value='undefined')
        call Test%checkForErr(Access%Err)
        Test%assertion = Access%value == "undefined"; call Test%verify()
        Test%assertion = Access%isSequential .eqv. .false.; call Test%verify()
        Test%assertion = Access%isDirect .eqv. .false.; call Test%verify()
        Test%assertion = Access%isUndefined .eqv. .true.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Access_type(value='Sequential')")
        Access = Access_type(value='Sequential')
        call Test%checkForErr(Access%Err)
        Test%assertion = Access%value == "sequential"; call Test%verify()
        Test%assertion = Access%isSequential .eqv. .true.; call Test%verify()
        Test%assertion = Access%isDirect .eqv. .false.; call Test%verify()
        Test%assertion = Access%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Access_type(value='Direct')")
        Access = Access_type(value='Direct')
        call Test%checkForErr(Access%Err)
        Test%assertion = Access%value == "direct"; call Test%verify()
        Test%assertion = Access%isSequential .eqv. .false.; call Test%verify()
        Test%assertion = Access%isDirect .eqv. .true.; call Test%verify()
        Test%assertion = Access%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Access_type(value='nonsense')")
        Access = Access_type(value='nonsense')
        Test%assertion = Access%Err%occurred .eqv. .true.; call Test%verify()
        Test%assertion = Access%value == ""; call Test%verify()
        Test%assertion = Access%isSequential .eqv. .false.; call Test%verify()
        Test%assertion = Access%isDirect .eqv. .false.; call Test%verify()
        Test%assertion = Access%isUndefined .eqv. .false.; call Test%verify()

    end subroutine test_Access_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_Form_type

        type(Form_type) :: Form

        if (Test%Image%isFirst) call Test%testing("Form_type()")
        Form = Form_type()
        call Test%checkForErr(Form%Err)
        Test%assertion = Form%value == "formatted"; call Test%verify()
        Test%assertion = Form%isFormatted .eqv. .true.; call Test%verify()
        Test%assertion = Form%isUnformatted .eqv. .false.; call Test%verify()
        Test%assertion = Form%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Form_type(value='undefined')")
        Form = Form_type(value='undefined')
        call Test%checkForErr(Form%Err)
        Test%assertion = Form%value == "undefined"; call Test%verify()
        Test%assertion = Form%isFormatted .eqv. .false.; call Test%verify()
        Test%assertion = Form%isUnformatted .eqv. .false.; call Test%verify()
        Test%assertion = Form%isUndefined .eqv. .true.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Form_type(value='Formatted')")
        Form = Form_type(value='Formatted')
        call Test%checkForErr(Form%Err)
        Test%assertion = Form%value == "formatted"; call Test%verify()
        Test%assertion = Form%isFormatted .eqv. .true.; call Test%verify()
        Test%assertion = Form%isUnformatted .eqv. .false.; call Test%verify()
        Test%assertion = Form%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Form_type(value='Unformatted')")
        Form = Form_type(value='Unformatted')
        call Test%checkForErr(Form%Err)
        Test%assertion = Form%value == "unformatted"; call Test%verify()
        Test%assertion = Form%isFormatted .eqv. .false.; call Test%verify()
        Test%assertion = Form%isUnformatted .eqv. .true.; call Test%verify()
        Test%assertion = Form%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Form_type(value='nonsense')")
        Form = Form_type(value='nonsense')
        Test%assertion = Form%Err%occurred .eqv. .true.; call Test%verify()
        Test%assertion = Form%value == ""; call Test%verify()
        Test%assertion = Form%isFormatted .eqv. .false.; call Test%verify()
        Test%assertion = Form%isUnformatted .eqv. .false.; call Test%verify()
        Test%assertion = Form%isUndefined .eqv. .false.; call Test%verify()

    end subroutine test_Form_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_Blank_type

        type(Blank_type) :: Blank

        if (Test%Image%isFirst) call Test%testing("Blank_type()")
        Blank = Blank_type()
        call Test%checkForErr(Blank%Err)
        Test%assertion = Blank%value == "null"; call Test%verify()
        Test%assertion = Blank%isNull .eqv. .true.; call Test%verify()
        Test%assertion = Blank%isZero .eqv. .false.; call Test%verify()
        Test%assertion = Blank%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Blank_type(value='undefined')")
        Blank = Blank_type(value='undefined')
        call Test%checkForErr(Blank%Err)
        Test%assertion = Blank%value == "undefined"; call Test%verify()
        Test%assertion = Blank%isNull .eqv. .false.; call Test%verify()
        Test%assertion = Blank%isZero .eqv. .false.; call Test%verify()
        Test%assertion = Blank%isUndefined .eqv. .true.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Blank_type(value='Null')")
        Blank = Blank_type(value='Null')
        call Test%checkForErr(Blank%Err)
        Test%assertion = Blank%value == "null"; call Test%verify()
        Test%assertion = Blank%isNull .eqv. .true.; call Test%verify()
        Test%assertion = Blank%isZero .eqv. .false.; call Test%verify()
        Test%assertion = Blank%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Blank_type(value='Zero')")
        Blank = Blank_type(value='Zero')
        call Test%checkForErr(Blank%Err)
        Test%assertion = Blank%value == "zero"; call Test%verify()
        Test%assertion = Blank%isNull .eqv. .false.; call Test%verify()
        Test%assertion = Blank%isZero .eqv. .true.; call Test%verify()
        Test%assertion = Blank%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Blank_type(value='nonsense')")
        Blank = Blank_type(value='nonsense')
        Test%assertion = Blank%Err%occurred .eqv. .true.; call Test%verify()
        Test%assertion = Blank%value == ""; call Test%verify()
        Test%assertion = Blank%isNull .eqv. .false.; call Test%verify()
        Test%assertion = Blank%isZero .eqv. .false.; call Test%verify()
        Test%assertion = Blank%isUndefined .eqv. .false.; call Test%verify()

    end subroutine test_Blank_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_Position_type

        type(Position_type) :: Position

        if (Test%Image%isFirst) call Test%testing("Position_type()")
        Position = Position_type()
        call Test%checkForErr(Position%Err)
        Test%assertion = Position%value == "asis"; call Test%verify()
        Test%assertion = Position%isRewind .eqv. .false.; call Test%verify()
        Test%assertion = Position%isAppend .eqv. .false.; call Test%verify()
        Test%assertion = Position%isAsis .eqv. .true.; call Test%verify()
        Test%assertion = Position%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Position_type(value='Undefined')")
        Position = Position_type(value='Undefined')
        call Test%checkForErr(Position%Err)
        Test%assertion = Position%value == "undefined"; call Test%verify()
        Test%assertion = Position%isRewind .eqv. .false.; call Test%verify()
        Test%assertion = Position%isAppend .eqv. .false.; call Test%verify()
        Test%assertion = Position%isAsis .eqv. .false.; call Test%verify()
        Test%assertion = Position%isUndefined .eqv. .true.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Position_type(value='Rewind')")
        Position = Position_type(value='Rewind')
        call Test%checkForErr(Position%Err)
        Test%assertion = Position%value == "rewind"; call Test%verify()
        Test%assertion = Position%isRewind .eqv. .true.; call Test%verify()
        Test%assertion = Position%isAppend .eqv. .false.; call Test%verify()
        Test%assertion = Position%isAsis .eqv. .false.; call Test%verify()
        Test%assertion = Position%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Position_type(value='Append')")
        Position = Position_type(value='Append')
        call Test%checkForErr(Position%Err)
        Test%assertion = Position%value == "append"; call Test%verify()
        Test%assertion = Position%isRewind .eqv. .false.; call Test%verify()
        Test%assertion = Position%isAppend .eqv. .true.; call Test%verify()
        Test%assertion = Position%isAsis .eqv. .false.; call Test%verify()
        Test%assertion = Position%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Position_type(value='Asis')")
        Position = Position_type(value='Asis')
        call Test%checkForErr(Position%Err)
        Test%assertion = Position%value == "asis"; call Test%verify()
        Test%assertion = Position%isRewind .eqv. .false.; call Test%verify()
        Test%assertion = Position%isAppend .eqv. .false.; call Test%verify()
        Test%assertion = Position%isAsis .eqv. .true.; call Test%verify()
        Test%assertion = Position%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Position_type(value='nonsense')")
        Position = Position_type(value='nonsense')
        Test%assertion = Position%Err%occurred .eqv. .true.; call Test%verify()
        Test%assertion = Position%value == ""; call Test%verify()
        Test%assertion = Position%isRewind .eqv. .false.; call Test%verify()
        Test%assertion = Position%isAppend .eqv. .false.; call Test%verify()
        Test%assertion = Position%isAsis .eqv. .false.; call Test%verify()
        Test%assertion = Position%isUndefined .eqv. .false.; call Test%verify()

    end subroutine test_Position_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_Delim_type

        type(Delim_type) :: Delim

        if (Test%Image%isFirst) call Test%testing("Delim_type()")
        Delim = Delim_type()
        call Test%checkForErr(Delim%Err)
        Test%assertion = Delim%value == "none"; call Test%verify()
        Test%assertion = Delim%isQuote .eqv. .false.; call Test%verify()
        Test%assertion = Delim%isApostrophe .eqv. .false.; call Test%verify()
        Test%assertion = Delim%isNone .eqv. .true.; call Test%verify()
        Test%assertion = Delim%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Delim_type(value='Undefined')")
        Delim = Delim_type(value='Undefined')
        call Test%checkForErr(Delim%Err)
        Test%assertion = Delim%value == "undefined"; call Test%verify()
        Test%assertion = Delim%isQuote .eqv. .false.; call Test%verify()
        Test%assertion = Delim%isApostrophe .eqv. .false.; call Test%verify()
        Test%assertion = Delim%isNone .eqv. .false.; call Test%verify()
        Test%assertion = Delim%isUndefined .eqv. .true.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Delim_type(value='Quote')")
        Delim = Delim_type(value='Quote')
        call Test%checkForErr(Delim%Err)
        Test%assertion = Delim%value == "quote"; call Test%verify()
        Test%assertion = Delim%isQuote .eqv. .true.; call Test%verify()
        Test%assertion = Delim%isApostrophe .eqv. .false.; call Test%verify()
        Test%assertion = Delim%isNone .eqv. .false.; call Test%verify()
        Test%assertion = Delim%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Delim_type(value='Apostrophe')")
        Delim = Delim_type(value='Apostrophe')
        call Test%checkForErr(Delim%Err)
        Test%assertion = Delim%value == "apostrophe"; call Test%verify()
        Test%assertion = Delim%isQuote .eqv. .false.; call Test%verify()
        Test%assertion = Delim%isApostrophe .eqv. .true.; call Test%verify()
        Test%assertion = Delim%isNone .eqv. .false.; call Test%verify()
        Test%assertion = Delim%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Delim_type(value='None')")
        Delim = Delim_type(value='None')
        call Test%checkForErr(Delim%Err)
        Test%assertion = Delim%value == "none"; call Test%verify()
        Test%assertion = Delim%isQuote .eqv. .false.; call Test%verify()
        Test%assertion = Delim%isApostrophe .eqv. .false.; call Test%verify()
        Test%assertion = Delim%isNone .eqv. .true.; call Test%verify()
        Test%assertion = Delim%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Delim_type(value='nonsense')")
        Delim = Delim_type(value='nonsense')
        Test%assertion = Delim%Err%occurred .eqv. .true.; call Test%verify()
        Test%assertion = Delim%value == ""; call Test%verify()
        Test%assertion = Delim%isQuote .eqv. .false.; call Test%verify()
        Test%assertion = Delim%isApostrophe .eqv. .false.; call Test%verify()
        Test%assertion = Delim%isNone .eqv. .false.; call Test%verify()
        Test%assertion = Delim%isUndefined .eqv. .false.; call Test%verify()

    end subroutine test_Delim_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_Round_type

        type(Round_type) :: Round

        if (Test%Image%isFirst) call Test%testing("Round_type()")
        Round = Round_type()
        call Test%checkForErr(Round%Err)
        Test%assertion = Round%value == "processor_defined"; call Test%verify()
        Test%assertion = Round%isUp .eqv. .false.; call Test%verify()
        Test%assertion = Round%isDown .eqv. .false.; call Test%verify()
        Test%assertion = Round%isZero .eqv. .false.; call Test%verify()
        Test%assertion = Round%isNearest .eqv. .false.; call Test%verify()
        Test%assertion = Round%isCompatible .eqv. .false.; call Test%verify()
        Test%assertion = Round%isProcessDefined .eqv. .true.; call Test%verify()
        Test%assertion = Round%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Round_type(value='Undefined')")
        Round = Round_type(value='Undefined')
        call Test%checkForErr(Round%Err)
        Test%assertion = Round%value == "undefined"; call Test%verify()
        Test%assertion = Round%isUp .eqv. .false.; call Test%verify()
        Test%assertion = Round%isDown .eqv. .false.; call Test%verify()
        Test%assertion = Round%isZero .eqv. .false.; call Test%verify()
        Test%assertion = Round%isNearest .eqv. .false.; call Test%verify()
        Test%assertion = Round%isCompatible .eqv. .false.; call Test%verify()
        Test%assertion = Round%isProcessDefined .eqv. .false.; call Test%verify()
        Test%assertion = Round%isUndefined .eqv. .true.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Round_type(value='UP')")
        Round = Round_type(value='UP')
        call Test%checkForErr(Round%Err)
        Test%assertion = Round%value == "up"; call Test%verify()
        Test%assertion = Round%isUp .eqv. .true.; call Test%verify()
        Test%assertion = Round%isDown .eqv. .false.; call Test%verify()
        Test%assertion = Round%isZero .eqv. .false.; call Test%verify()
        Test%assertion = Round%isNearest .eqv. .false.; call Test%verify()
        Test%assertion = Round%isCompatible .eqv. .false.; call Test%verify()
        Test%assertion = Round%isProcessDefined .eqv. .false.; call Test%verify()
        Test%assertion = Round%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Round_type(value='Down')")
        Round = Round_type(value='Down')
        call Test%checkForErr(Round%Err)
        Test%assertion = Round%value == "down"; call Test%verify()
        Test%assertion = Round%isUp .eqv. .false.; call Test%verify()
        Test%assertion = Round%isDown .eqv. .true.; call Test%verify()
        Test%assertion = Round%isZero .eqv. .false.; call Test%verify()
        Test%assertion = Round%isNearest .eqv. .false.; call Test%verify()
        Test%assertion = Round%isCompatible .eqv. .false.; call Test%verify()
        Test%assertion = Round%isProcessDefined .eqv. .false.; call Test%verify()
        Test%assertion = Round%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Round_type(value='Zero')")
        Round = Round_type(value='Zero')
        call Test%checkForErr(Round%Err)
        Test%assertion = Round%value == "zero"; call Test%verify()
        Test%assertion = Round%isUp .eqv. .false.; call Test%verify()
        Test%assertion = Round%isDown .eqv. .false.; call Test%verify()
        Test%assertion = Round%isZero .eqv. .true.; call Test%verify()
        Test%assertion = Round%isNearest .eqv. .false.; call Test%verify()
        Test%assertion = Round%isCompatible .eqv. .false.; call Test%verify()
        Test%assertion = Round%isProcessDefined .eqv. .false.; call Test%verify()
        Test%assertion = Round%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Round_type(value='Nearest')")
        Round = Round_type(value='Nearest')
        call Test%checkForErr(Round%Err)
        Test%assertion = Round%value == "nearest"; call Test%verify()
        Test%assertion = Round%isUp .eqv. .false.; call Test%verify()
        Test%assertion = Round%isDown .eqv. .false.; call Test%verify()
        Test%assertion = Round%isZero .eqv. .false.; call Test%verify()
        Test%assertion = Round%isNearest .eqv. .true.; call Test%verify()
        Test%assertion = Round%isCompatible .eqv. .false.; call Test%verify()
        Test%assertion = Round%isProcessDefined .eqv. .false.; call Test%verify()
        Test%assertion = Round%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Round_type(value='Compatible')")
        Round = Round_type(value='Compatible')
        call Test%checkForErr(Round%Err)
        Test%assertion = Round%value == "compatible"; call Test%verify()
        Test%assertion = Round%isUp .eqv. .false.; call Test%verify()
        Test%assertion = Round%isDown .eqv. .false.; call Test%verify()
        Test%assertion = Round%isZero .eqv. .false.; call Test%verify()
        Test%assertion = Round%isNearest .eqv. .false.; call Test%verify()
        Test%assertion = Round%isCompatible .eqv. .true.; call Test%verify()
        Test%assertion = Round%isProcessDefined .eqv. .false.; call Test%verify()
        Test%assertion = Round%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Round_type(value='Processor_defined')")
        Round = Round_type(value='Processor_defined')
        call Test%checkForErr(Round%Err)
        Test%assertion = Round%value == "processor_defined"; call Test%verify()
        Test%assertion = Round%isUp .eqv. .false.; call Test%verify()
        Test%assertion = Round%isDown .eqv. .false.; call Test%verify()
        Test%assertion = Round%isZero .eqv. .false.; call Test%verify()
        Test%assertion = Round%isNearest .eqv. .false.; call Test%verify()
        Test%assertion = Round%isCompatible .eqv. .false.; call Test%verify()
        Test%assertion = Round%isProcessDefined .eqv. .true.; call Test%verify()
        Test%assertion = Round%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Round_type(value='nonsense')")
        Round = Round_type(value='nonsense')
        Test%assertion = Round%Err%occurred .eqv. .true.; call Test%verify()
        Test%assertion = Round%value == ""; call Test%verify()
        Test%assertion = Round%isUp .eqv. .false.; call Test%verify()
        Test%assertion = Round%isDown .eqv. .false.; call Test%verify()
        Test%assertion = Round%isZero .eqv. .false.; call Test%verify()
        Test%assertion = Round%isNearest .eqv. .false.; call Test%verify()
        Test%assertion = Round%isCompatible .eqv. .false.; call Test%verify()
        Test%assertion = Round%isProcessDefined .eqv. .false.; call Test%verify()
        Test%assertion = Round%isUndefined .eqv. .false.; call Test%verify()

    end subroutine test_Round_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_Sign_type

        type(Sign_type) :: Sign

        if (Test%Image%isFirst) call Test%testing("Sign_type()")
        Sign = Sign_type()
        call Test%checkForErr(Sign%Err)
        Test%assertion = Sign%value == "processor_defined"; call Test%verify()
        Test%assertion = Sign%isSuppress .eqv. .false.; call Test%verify()
        Test%assertion = Sign%isPlus .eqv. .false.; call Test%verify()
        Test%assertion = Sign%isProcessDefined .eqv. .true.; call Test%verify()
        Test%assertion = Sign%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Sign_type(value='Undefined')")
        Sign = Sign_type(value='Undefined')
        call Test%checkForErr(Sign%Err)
        Test%assertion = Sign%value == "undefined"; call Test%verify()
        Test%assertion = Sign%isSuppress .eqv. .false.; call Test%verify()
        Test%assertion = Sign%isPlus .eqv. .false.; call Test%verify()
        Test%assertion = Sign%isProcessDefined .eqv. .false.; call Test%verify()
        Test%assertion = Sign%isUndefined .eqv. .true.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Sign_type(value='Suppress')")
        Sign = Sign_type(value='Suppress')
        call Test%checkForErr(Sign%Err)
        Test%assertion = Sign%value == "suppress"; call Test%verify()
        Test%assertion = Sign%isSuppress .eqv. .true.; call Test%verify()
        Test%assertion = Sign%isPlus .eqv. .false.; call Test%verify()
        Test%assertion = Sign%isProcessDefined .eqv. .false.; call Test%verify()
        Test%assertion = Sign%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Sign_type(value='Plus')")
        Sign = Sign_type(value='Plus')
        call Test%checkForErr(Sign%Err)
        Test%assertion = Sign%value == "plus"; call Test%verify()
        Test%assertion = Sign%isSuppress .eqv. .false.; call Test%verify()
        Test%assertion = Sign%isPlus .eqv. .true.; call Test%verify()
        Test%assertion = Sign%isProcessDefined .eqv. .false.; call Test%verify()
        Test%assertion = Sign%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Sign_type(value='Processor_defined')")
        Sign = Sign_type(value='Processor_defined')
        call Test%checkForErr(Sign%Err)
        Test%assertion = Sign%value == "processor_defined"; call Test%verify()
        Test%assertion = Sign%isSuppress .eqv. .false.; call Test%verify()
        Test%assertion = Sign%isPlus .eqv. .false.; call Test%verify()
        Test%assertion = Sign%isProcessDefined .eqv. .true.; call Test%verify()
        Test%assertion = Sign%isUndefined .eqv. .false.; call Test%verify()

        if (Test%Image%isFirst) call Test%testing("Sign_type(value='nonsense')")
        Sign = Sign_type(value='nonsense')
        Test%assertion = Sign%Err%occurred .eqv. .true.; call Test%verify()
        Test%assertion = Sign%value == ""; call Test%verify()
        Test%assertion = Sign%isSuppress .eqv. .false.; call Test%verify()
        Test%assertion = Sign%isPlus .eqv. .false.; call Test%verify()
        Test%assertion = Sign%isProcessDefined .eqv. .false.; call Test%verify()
        Test%assertion = Sign%isUndefined .eqv. .false.; call Test%verify()

    end subroutine test_Sign_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getFileExistStatus()
        if (Test%Image%isFirst) call Test%testing("@getFileExistStatus()")
        call Test%skipping()
    end subroutine test_getFileExistStatus

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getFileOpenStatus()
        if (Test%Image%isFirst) call Test%testing("@getFileOpenStatus()")
        call Test%skipping()
    end subroutine test_getFileOpenStatus

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getFileNumber()
        if (Test%Image%isFirst) call Test%testing("@getFileNumber()")
        call Test%skipping()
    end subroutine test_getFileNumber

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getFileName()
        if (Test%Image%isFirst) call Test%testing("@test_getFileName()")
        call Test%skipping()
    end subroutine test_getFileName

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getFileAccess()
        if (Test%Image%isFirst) call Test%testing("@test_getFileAccess()")
        call Test%skipping()
    end subroutine test_getFileAccess

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getFileForm()
        if (Test%Image%isFirst) call Test%testing("@test_getFileForm()")
        call Test%skipping()
    end subroutine test_getFileForm

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getFileRecl()
        if (Test%Image%isFirst) call Test%testing("@test_getFileRecl()")
        call Test%skipping()
    end subroutine test_getFileRecl

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getFileBlank()
        if (Test%Image%isFirst) call Test%testing("@test_getFileBlank()")
        call Test%skipping()
    end subroutine test_getFileBlank

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getFilePosition()
        if (Test%Image%isFirst) call Test%testing("@test_getFilePosition()")
        call Test%skipping()
    end subroutine test_getFilePosition

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getFileAction()
        if (Test%Image%isFirst) call Test%testing("@test_getFileAction()")
        call Test%skipping()
    end subroutine test_getFileAction

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getFileDelim()
        if (Test%Image%isFirst) call Test%testing("@test_getFileDelim()")
        call Test%skipping()
    end subroutine test_getFileDelim

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getWriteErr()
        if (Test%Image%isFirst) call Test%testing("@test_getWriteErr()")
        call Test%skipping()
    end subroutine test_getWriteErr

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getReadErr()
        if (Test%Image%isFirst) call Test%testing("@test_getReadErr()")
        call Test%skipping()
    end subroutine test_getReadErr

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getCloseErr()
        if (Test%Image%isFirst) call Test%testing("@test_getCloseErr()")
        call Test%skipping()
    end subroutine test_getCloseErr

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getOpenErr()
        if (Test%Image%isFirst) call Test%testing("@test_getOpenErr()")
        call Test%skipping()
    end subroutine test_getOpenErr

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getInqErr()
        if (Test%Image%isFirst) call Test%testing("@test_getInqErr()")
        call Test%skipping()
    end subroutine test_getInqErr

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_File_mod