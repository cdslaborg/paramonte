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
!!!!       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This module contains tests of the module [String_mod](@ref string_mod).
!>  \author Amir Shahmoradi

module Test_String_mod

    use String_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_String

    type(Test_type) :: Test

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine test_String()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)
        call Test%run(test_split_1, "test_split_1")
        call Test%run(test_split_2, "test_split_2")
        call Test%run(test_split_3, "test_split_3")
        call Test%run(test_split_4, "test_split_4")
        call Test%run(test_split_5, "test_split_5")
        call Test%run(test_num2str_1, "test_num2str_1")
        call Test%run(test_num2str_2, "test_num2str_2")
        call Test%run(test_isDigit_1, "test_isDigit_1")
        call Test%run(test_str2int_1, "test_str2int_1")
        call Test%run(test_str2int_2, "test_str2int_2")
        call Test%run(test_log2str_1, "test_log2str_1")
        call Test%run(test_log2str_2, "test_log2str_2")
        call Test%run(test_str2real_1, "test_str2real_1")
        call Test%run(test_str2real_2, "test_str2real_2")
        call Test%run(test_splitStr_1, "test_splitStr_1")
        call Test%run(test_splitStr_2, "test_splitStr_2")
        call Test%run(test_splitStr_3, "test_splitStr_3")
        call Test%run(test_splitStr_4, "test_splitStr_4")
        call Test%run(test_splitStr_5, "test_splitStr_5")
        call Test%run(test_padString_1, "test_padString_1")
        call Test%run(test_padString_2, "test_padString_2")
        call Test%run(test_padString_3, "test_padString_3")
        call Test%run(test_isInteger_1, "test_isInteger_1")
        call Test%run(test_str2int32_1, "test_str2int32_1")
        call Test%run(test_str2int32_2, "test_str2int32_2")
        call Test%run(test_str2int64_1, "test_str2int64_1")
        call Test%run(test_str2int64_2, "test_str2int64_2")
        call Test%run(test_int322str_1, "test_int322str_1")
        call Test%run(test_int322str_2, "test_int322str_2")
        call Test%run(test_int322str_3, "test_int322str_3")
        call Test%run(test_int642str_1, "test_int642str_1")
        call Test%run(test_int642str_2, "test_int642str_2")
        call Test%run(test_int642str_3, "test_int642str_3")
        call Test%run(test_real322str_1, "test_real322str_1")
        call Test%run(test_real322str_2, "test_real322str_2")
        call Test%run(test_real322str_3, "test_real322str_3")
        call Test%run(test_real322str_4, "test_real322str_4")
        call Test%run(test_real642str_1, "test_real642str_1")
        call Test%run(test_real642str_2, "test_real642str_2")
        call Test%run(test_real642str_3, "test_real642str_3")
        call Test%run(test_real642str_4, "test_real642str_4")
        call Test%run(test_str2real32_1, "test_str2real32_1")
        call Test%run(test_str2real32_2, "test_str2real32_2")
        call Test%run(test_str2real64_1, "test_str2real64_1")
        call Test%run(test_str2real64_2, "test_str2real64_2")
        call Test%run(test_replaceStr_1, "test_replaceStr_1")
        call Test%run(test_replaceStr_2, "test_replaceStr_2")
        call Test%run(test_IntStr_type_1, "test_IntStr_type_1")
        call Test%run(test_IntStr_type_2, "test_IntStr_type_2")
        call Test%run(test_getLowerCase_1, "test_getLowerCase_1")
        call Test%run(test_getLowerCase_2, "test_getLowerCase_2")
        call Test%run(test_getUpperCase_1, "test_getUpperCase_1")
        call Test%run(test_getUpperCase_2, "test_getUpperCase_2")
        call Test%run(test_real642str_1D_1, "test_real642str_1D_1")
        call Test%run(test_real642str_1D_2, "test_real642str_1D_2")
        call Test%run(test_real642str_1D_3, "test_real642str_1D_3")
        call Test%run(test_real642str_1D_4, "test_real642str_1D_4")
        call Test%run(test_real642str_2D_1, "test_real642str_2D_1")
        call Test%run(test_real642str_2D_2, "test_real642str_2D_2")
        call Test%run(test_real642str_2D_3, "test_real642str_2D_3")
        call Test%run(test_real642str_2D_4, "test_real642str_2D_4")
        call Test%finalize()

    end subroutine test_String

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_split_1() result(assertion)

        implicit none
        logical             :: assertion
        type(String_type)   :: String
        integer             :: i

        String%value = "  StringString "
        String%Parts  = split(string=String%value,delim="String")
        assertion = String%Parts(1)%record == "  " .and. String%Parts(2)%record == " "

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "String%value = '", String%value, "'"
            write(Test%outputUnit,"(*(g0))") "split(string=String%value,delim='String') = ", (String%Parts(i)%record,i=1,size(String%Parts))
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_split_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_split_2() result(assertion)

        implicit none
        logical             :: assertion
        type(String_type)   :: String
        integer             :: i

        String%value = "  Stringstring "
        String%Parts  = split(string=String%value,delim='str')
        assertion = String%Parts(1)%record == "  String" .and. String%Parts(2)%record == "ing "

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "String%value = '", String%value, "'"
            write(Test%outputUnit,"(*(g0))") "split(string=String%value,delim='str') = ", (String%Parts(i)%record,i=1,size(String%Parts))
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_split_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_split_3() result(assertion)

        implicit none
        logical             :: assertion
        type(String_type)   :: String
        integer             :: i

        String%value = "   "
        String%Parts  = split(string=String%value,delim=' ')
        assertion = String%Parts(1)%record == ""

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "String%value = '", String%value, "'"
            write(Test%outputUnit,"(*(g0))") "split(string=String%value,delim=' ') = ", ("'"//String%Parts(i)%record//"'",i=1,size(String%Parts))
            write(Test%outputUnit,"(A,*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_split_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_split_4() result(assertion)

        implicit none
        logical             :: assertion
        type(String_type)   :: String
        integer             :: i

        String%value = "   "
        String%Parts  = split(string = String%value, delim = " ", nPart = String%nPart)
        assertion = String%nPart == 4_IK

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "String%value = '", String%value, "'"
            write(Test%outputUnit,"(*(g0))") "split(string=String%value,delim=' ') = ", ("'"//String%Parts(i)%record//"'",i=1,size(String%Parts))
            write(Test%outputUnit,"(*(g0))") "String%nPart = ", String%nPart
            write(Test%outputUnit,"(A,*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_split_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> When the input delim is empty, the whole string, split character by character, should be returned.
    function test_split_5() result(assertion)

        implicit none
        logical             :: assertion
        type(String_type)   :: String
        integer             :: i

        String%value = "ParaMonte"
        String%Parts  = split(string = String%value, delim = "", nPart = String%nPart)
        assertion = String%nPart == len(String%value)
        do i = 1, len(String%value)
            assertion = assertion .and. string%parts(i)%record == string%value(i:i)
        end do

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "String%value = '", String%value, "'"
            write(Test%outputUnit,"(*(g0))") "split(string=String%value,delim=' ') = ", ("'"//String%Parts(i)%record//"'",i=1,size(String%Parts))
            write(Test%outputUnit,"(*(g0))") "String%nPart = ", String%nPart
            write(Test%outputUnit,"(A,*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_split_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_splitStr_1() result(assertion)

        implicit none
        logical             :: assertion
        type(String_type)   :: String
        integer             :: i

        String%value = "  StringString "
        String%Parts  = splitStr(string=String%value,delimiter="String")
        assertion = String%Parts(1)%record == "  " .and. String%Parts(2)%record == " "

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "String%value = '", String%value, "'"
            write(Test%outputUnit,"(*(g0))") "splitStr(string=String%value,delimiter='String') = ", (String%Parts(i)%record,i=1,size(String%Parts))
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_splitStr_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_splitStr_2() result(assertion)

        implicit none
        logical             :: assertion
        type(String_type)   :: String
        integer             :: i

        String%value = "  Stringstring "
        String%Parts  = splitStr(string=String%value,delimiter='str')
        assertion = String%Parts(1)%record == "  String" .and. String%Parts(2)%record == "ing "

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "String%value = '", String%value, "'"
            write(Test%outputUnit,"(*(g0))") "splitStr(string=String%value,delimiter='str') = ", (String%Parts(i)%record,i=1,size(String%Parts))
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_splitStr_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_splitStr_3() result(assertion)

        implicit none
        logical             :: assertion
        type(String_type)   :: String
        integer             :: i

        String%value = "   "
        String%Parts  = splitStr(string=String%value,delimiter=' ')
        assertion = String%Parts(1)%record == ""

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "String%value = '", String%value, "'"
            write(Test%outputUnit,"(*(g0))") "splitStr(string=String%value,delimiter=' ') = ", ("'"//String%Parts(i)%record//"'",i=1,size(String%Parts))
            write(Test%outputUnit,"(A,*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_splitStr_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_splitStr_4() result(assertion)

        implicit none
        logical             :: assertion
        type(String_type)   :: String
        integer             :: i

        String%value = "   "
        String%Parts  = splitStr(string = String%value, delimiter = " ", nPart = String%nPart)
        assertion = String%nPart == 1_IK

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "String%value = '", String%value, "'"
            write(Test%outputUnit,"(*(g0))") "splitStr(string=String%value,delimiter=' ') = ", ("'"//String%Parts(i)%record//"'",i=1,size(String%Parts))
            write(Test%outputUnit,"(*(g0))") "String%nPart = ", String%nPart
            write(Test%outputUnit,"(A,*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_splitStr_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> When the input delim is empty, the whole string should be returned.
    function test_splitStr_5() result(assertion)

        implicit none
        logical             :: assertion
        type(String_type)   :: String
        integer             :: i

        String%value = "ParaMonte"
        String%Parts  = splitStr(string = String%value, delimiter = "", nPart = String%nPart)
        assertion = String%nPart == 1_IK .and. string%parts(1)%record == string%value

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "String%value = '", String%value, "'"
            write(Test%outputUnit,"(*(g0))") "splitStr(string=String%value,delim=' ') = ", ("'"//String%Parts(i)%record//"'",i=1,size(String%Parts))
            write(Test%outputUnit,"(*(g0))") "String%nPart = ", String%nPart
            write(Test%outputUnit,"(A,*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_splitStr_5

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_str2int_1() result(assertion)

        use Constants_mod, only: IK
        implicit none
        logical                 :: assertion
        character(*), parameter :: this_str = "123456"
        integer(IK) , parameter :: this_ref = 123456_IK
        integer(IK)             :: this

        this = str2int(this_str)
        assertion = this == this_ref

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_str  = '", this_str, "'"
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_str2int_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_str2int_2() result(assertion)

        use Constants_mod, only: IK
        implicit none
        logical                 :: assertion
        character(*), parameter :: this_str = "123456"
        integer(IK) , parameter :: this_ref = 123456_IK
        integer(IK)             :: this
        integer                 :: iostat

        this = str2int(this_str, iostat)
        assertion = iostat /= 0 .or. (iostat == 0 .and. this == this_ref)

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_str  = '", this_str, "'"
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_str2int_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_str2int32_1() result(assertion)

        use, intrinsic :: iso_fortran_env, only: IK => int32
        implicit none
        logical                 :: assertion
        character(*), parameter :: this_str = "123456"
        integer(IK) , parameter :: this_ref = 123456_IK
        integer(IK)             :: this

        this = str2int32(this_str)
        assertion = this == this_ref

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_str  = '", this_str, "'"
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_str2int32_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_str2int32_2() result(assertion)

        use, intrinsic :: iso_fortran_env, only: IK => int32
        implicit none
        logical                 :: assertion
        character(*), parameter :: this_str = "123456"
        integer(IK) , parameter :: this_ref = 123456_IK
        integer(IK)             :: this
        integer                 :: iostat

        this = str2int32(this_str, iostat)
        assertion = iostat /= 0 .or. (iostat == 0 .and. this == this_ref)

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_str  = '", this_str, "'"
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_str2int32_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_str2int64_1() result(assertion)

        use, intrinsic :: iso_fortran_env, only: IK => int64
        implicit none
        logical                 :: assertion
        character(*), parameter :: this_str = "123456"
        integer(IK) , parameter :: this_ref = 123456_IK
        integer(IK)             :: this

        this = str2int64(this_str)
        assertion = this == this_ref

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_str  = '", this_str, "'"
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_str2int64_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_str2int64_2() result(assertion)

        use, intrinsic :: iso_fortran_env, only: IK => int64
        implicit none
        logical                 :: assertion
        character(*), parameter :: this_str = "1234561284676"
        integer(IK) , parameter :: this_ref = 1234561284676_IK
        integer(IK)             :: this
        integer                 :: iostat

        this = str2int64(this_str, iostat)
        assertion = iostat /= 0 .or. (iostat == 0 .and. this == this_ref)

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_str  = '", this_str, "'"
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_str2int64_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_str2real_1() result(assertion)

        use Constants_mod, only: RK
        implicit none
        logical                 :: assertion
        character(*), parameter :: this_str = "1.1234567890123456789e-30"
        real(RK)    , parameter :: this_ref = 1.1234567890123456789e-30_RK
        real(RK)                :: this

        this = str2real(this_str)
        assertion = this == this_ref

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_str  = '", this_str, "'"
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_str2real_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_str2real_2() result(assertion)

        use Constants_mod, only: RK
        implicit none
        logical                 :: assertion
        character(*), parameter :: this_str = "1.1234567890123456789e-30"
        real(RK)    , parameter :: this_ref = 1.1234567890123456789e-30_RK
        real(RK)                :: this
        integer                 :: iostat

        this = str2real(this_str, iostat)
        assertion = iostat /= 0 .or. (iostat == 0 .and. this == this_ref)

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_str  = '", this_str, "'"
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_str2real_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_str2real32_1() result(assertion)

        use, intrinsic :: iso_fortran_env, only: RK => real32
        implicit none
        logical                 :: assertion
        character(*), parameter :: this_str = "1.1234567890123456789e-30"
        real(RK)    , parameter :: this_ref = 1.1234567890123456789e-30_RK
        real(RK)                :: this

        this = str2real32(this_str)
        assertion = this == this_ref

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_str  = '", this_str, "'"
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_str2real32_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_str2real32_2() result(assertion)

        use, intrinsic :: iso_fortran_env, only: RK => real32
        implicit none
        logical                 :: assertion
        character(*), parameter :: this_str = "1.1234567890123456789e-30"
        real(RK)    , parameter :: this_ref = 1.1234567890123456789e-30_RK
        real(RK)                :: this
        integer                 :: iostat

        this = str2real32(this_str, iostat)
        assertion = iostat /= 0 .or. (iostat == 0 .and. this == this_ref)

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_str  = '", this_str, "'"
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_str2real32_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_str2real64_1() result(assertion)

        use, intrinsic :: iso_fortran_env, only: RK => real64
        implicit none
        logical                 :: assertion
        character(*), parameter :: this_str = "1.1234567890123456789e-30"
        real(RK)    , parameter :: this_ref = 1.1234567890123456789e-30_RK
        real(RK)                :: this

        this = str2real64(this_str)
        assertion = this == this_ref

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_str  = '", this_str, "'"
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_str2real64_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_str2real64_2() result(assertion)

        use, intrinsic :: iso_fortran_env, only: RK => real64
        implicit none
        logical                 :: assertion
        character(*), parameter :: this_str = "1.1234567890123456789e-30"
        real(RK)    , parameter :: this_ref = 1.1234567890123456789e-30_RK
        real(RK)                :: this
        integer                 :: iostat

        this = str2real64(this_str, iostat)
        assertion = iostat /= 0 .or. (iostat == 0 .and. this == this_ref)

        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_str  = '", this_str, "'"
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP

    end function test_str2real64_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_replaceStr_1() result(assertion)
        implicit none
        logical                     :: assertion
        type(String_type)           :: String
        character(:), allocatable   :: strout
        String%value = "  StringString "
        strout = String%replaceStr(string=String%value, search="Str", substitute="")
        assertion = strout == "  inging "
    end function test_replaceStr_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_replaceStr_2() result(assertion)

        implicit none
        logical                     :: assertion
        type(String_type)           :: String
        character(:), allocatable   :: strout
        String%value = ""
        strout = String%replaceStr(string=String%value, search="", substitute="")
        assertion = strout == ""
    end function test_replaceStr_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLowerCase_1()  result(assertion)
        implicit none
        logical                     :: assertion
        type(String_type)           :: String
        character(:), allocatable   :: strout
        String%value = "  StringString !@#$"
        strout = String%getLowerCase(string=String%value)
        assertion = strout == "  stringstring !@#$"
    end function test_getLowerCase_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLowerCase_2()  result(assertion)
        implicit none
        logical                     :: assertion
        type(String_type)           :: String
        character(:), allocatable   :: strout
        String%value = ""
        strout = String%getLowerCase(string=String%value)
        assertion = strout == ""
    end function test_getLowerCase_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getUpperCase_1() result(assertion)
        implicit none
        logical                     :: assertion
        type(String_type)           :: String
        character(:), allocatable   :: strout
        String%value = "  StringString !@#$"
        strout = String%getUpperCase(string=String%value)
        assertion = strout == "  STRINGSTRING !@#$"
    end function test_getUpperCase_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getUpperCase_2() result(assertion)
        implicit none
        logical                     :: assertion
        type(String_type)           :: String
        character(:), allocatable   :: strout
        String%value = ""
        strout = String%getUpperCase(string=String%value)
        assertion = strout == ""
    end function test_getUpperCase_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isInteger_1() result(assertion)
        implicit none
        logical                     :: assertion
        type(String_type)           :: String
        assertion = String%isInteger("1283985")
    end function test_isInteger_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isDigit_1() result(assertion)
        implicit none
        logical                     :: assertion
        type(String_type)           :: String
        integer                     :: i
        assertion = .true.
        do i = 0, 9
            assertion = assertion .and. String%isDigit(num2str(i))
        end do
    end function test_isDigit_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_int322str_1() result(assertion)
        use, intrinsic :: iso_fortran_env, only: IK => int32
        implicit none
        logical                     :: assertion
        character(*), parameter     :: this_ref = "12345"
        integer(IK) , parameter     :: this_int = 12345_IK
        character(:), allocatable   :: this
        this = num2str(this_int)
        assertion = this == this_ref .and. len(this) == len(this_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = '", this_ref, "'"
            write(Test%outputUnit,"(*(g0))") "this_str  = " , this_int
            write(Test%outputUnit,"(*(g0))") "this      = '", this, "'"
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_int322str_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_int322str_2() result(assertion)
        use, intrinsic :: iso_fortran_env, only: IK => int32
        implicit none
        logical                     :: assertion
        character(*), parameter     :: this_ref = "0000012345"
        integer(IK) , parameter     :: this_int = 12345_IK
        character(:), allocatable   :: this
        this = num2str(this_int,"(1I10.10)")
        assertion = this == this_ref .and. len(this) == len(this_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = '", this_ref, "'"
            write(Test%outputUnit,"(*(g0))") "this_str  = " , this_int
            write(Test%outputUnit,"(*(g0))") "this      = '", this, "'"
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_int322str_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_int322str_3() result(assertion)
        use, intrinsic :: iso_fortran_env, only: IK => int32
        implicit none
        logical                     :: assertion
        character(*), parameter     :: this_ref = "12345     "
        integer(IK) , parameter     :: this_int = 12345_IK
        character(:), allocatable   :: this
        this = num2str(this_int,minlen=10_IK)
        assertion = this == this_ref .and. len(this) == len(this_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = '", this_ref, "'"
            write(Test%outputUnit,"(*(g0))") "this_str  = " , this_int
            write(Test%outputUnit,"(*(g0))") "this      = '", this, "'"
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_int322str_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_int642str_1() result(assertion)
        use, intrinsic :: iso_fortran_env, only: IK => int64
        implicit none
        logical                     :: assertion
        character(*), parameter     :: this_ref = "98765432112345"
        integer(IK) , parameter     :: this_int = 98765432112345_IK
        character(:), allocatable   :: this
        this = num2str(this_int)
        assertion = this == this_ref .and. len(this) == len(this_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = '", this_ref, "'"
            write(Test%outputUnit,"(*(g0))") "this_str  = " , this_int
            write(Test%outputUnit,"(*(g0))") "this      = '", this, "'"
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_int642str_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_int642str_2() result(assertion)
        use, intrinsic :: iso_fortran_env, only: IK => int64
        implicit none
        logical                     :: assertion
        character(*), parameter     :: this_ref = "00000098765432112345"
        integer(IK) , parameter     :: this_int = 98765432112345_IK
        character(:), allocatable   :: this
        this = num2str(this_int,"(1I20.20)")
        assertion = this == this_ref .and. len(this) == len(this_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = '", this_ref, "'"
            write(Test%outputUnit,"(*(g0))") "this_str  = " , this_int
            write(Test%outputUnit,"(*(g0))") "this      = '", this, "'"
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_int642str_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_int642str_3() result(assertion)
        use, intrinsic :: iso_fortran_env, only: IK => int64, int32
        implicit none
        logical                     :: assertion
        character(*), parameter     :: this_ref = "98765432112345      "
        integer(IK) , parameter     :: this_int = 98765432112345_IK
        character(:), allocatable   :: this
        this = num2str(this_int,minlen=20_int32)
        assertion = this == this_ref .and. len(this) == len(this_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = '", this_ref, "'"
            write(Test%outputUnit,"(*(g0))") "this_str  = " , this_int
            write(Test%outputUnit,"(*(g0))") "this      = '", this, "'"
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_int642str_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_real322str_1() result(assertion)
        use, intrinsic :: iso_fortran_env, only: RK => real32
        implicit none
        logical                     :: assertion
        real(RK)    , parameter     :: this_ref = 1.23456798e-30_RK
        real(RK)                    :: this
        this = str2real32( num2str(this_ref) )
        assertion = this == this_ref
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_real322str_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_real322str_2() result(assertion)
        use, intrinsic :: iso_fortran_env, only: RK => real32
        implicit none
        logical                     :: assertion
        real(RK)    , parameter     :: this_ref = 1.23456798e-30_RK
        real(RK)                    :: this
        this = str2real32( num2str(this_ref,"(g0)") )
        assertion = this == this_ref
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_real322str_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_real322str_3() result(assertion)
        use, intrinsic :: iso_fortran_env, only: RK => real32, IK => int32
        implicit none
        logical                     :: assertion
        real(RK)    , parameter     :: this_ref = 1.23456798e-30_RK
        real(RK)                    :: this
        this = str2real32( num2str(this_ref,"(g0)",20_IK) )
        assertion = this == this_ref
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_real322str_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> minLen can be larger than the length of the constructed string.
    function test_real322str_4() result(assertion)
        use, intrinsic :: iso_fortran_env, only: RK => real32, IK => int32
        implicit none
        logical                     :: assertion
        real(RK)    , parameter     :: this_ref = 1.23456798e-30_RK
        real(RK)                    :: this
        this = str2real32( num2str(this_ref,"(g0)",263_IK) )
        assertion = this == this_ref
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_real322str_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_real642str_1() result(assertion)
        use, intrinsic :: iso_fortran_env, only: RK => real64
        implicit none
        logical                     :: assertion
        real(RK)    , parameter     :: this_ref = 1.23456798e-30_RK
        real(RK)                    :: this
        this = str2real64( num2str(this_ref) )
        assertion = this == this_ref
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_real642str_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_real642str_2() result(assertion)
        use, intrinsic :: iso_fortran_env, only: RK => real64
        implicit none
        logical                     :: assertion
        real(RK)    , parameter     :: this_ref = 1.23456798e-30_RK
        real(RK)                    :: this
        this = str2real64( num2str(this_ref,"(g0)") )
        assertion = this == this_ref
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_real642str_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_real642str_3() result(assertion)
        use, intrinsic :: iso_fortran_env, only: RK => real64, IK => int32
        implicit none
        logical                     :: assertion
        real(RK)    , parameter     :: this_ref = 1.23456798e-30_RK
        real(RK)                    :: this
        this = str2real64( num2str(this_ref,"(g0)",30_IK) )
        assertion = this == this_ref
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_real642str_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> minLen can be larger than the length of the constructed string.
    function test_real642str_4() result(assertion)
        use, intrinsic :: iso_fortran_env, only: RK => real64, IK => int32
        implicit none
        logical                     :: assertion
        real(RK)    , parameter     :: this_ref = 1.23456798e-30_RK
        real(RK)                    :: this
        this = str2real64( num2str(this_ref,"(g0)",263_IK) )
        assertion = this == this_ref
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_real642str_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_real642str_1D_1() result(assertion)
        use, intrinsic :: iso_fortran_env, only: RK => real64
        implicit none
        logical                     :: assertion
        integer     , parameter     :: nthis = 2
        real(RK)    , parameter     :: this_ref(nthis) = [1.23456798e-30_RK, 2.32456798e+30_RK]
        real(RK)                    :: this(nthis)
        character(:), allocatable   :: string
        assertion = .true.
        string = num2str(this_ref)
        read(string,*) this
        assertion = all(this == this_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_real642str_1D_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_real642str_1D_2() result(assertion)
        use, intrinsic :: iso_fortran_env, only: RK => real64
        implicit none
        logical                     :: assertion
        integer     , parameter     :: nthis = 2
        real(RK)    , parameter     :: this_ref(nthis) = [1.23456798e-30_RK, 2.32456798e+30_RK]
        real(RK)                    :: this(nthis)
        character(:), allocatable   :: string
        assertion = .true.
        string = num2str(this_ref,"(*(g0,:,' '))")
        read(string,*) this
        assertion = all(this == this_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_real642str_1D_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_real642str_1D_3() result(assertion)
        use, intrinsic :: iso_fortran_env, only: RK => real64, IK => int32
        implicit none
        logical                     :: assertion
        integer     , parameter     :: nthis = 2
        real(RK)    , parameter     :: this_ref(nthis) = [1.23456798e-30_RK, 2.32456798e+30_RK]
        real(RK)                    :: this(nthis)
        character(:), allocatable   :: string
        assertion = .true.
        string = num2str(this_ref,"(*(g0,:,' '))",63_IK)
        read(string,*) this
        assertion = all(this == this_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_real642str_1D_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> minLen can be larger than the length of the constructed string.
    function test_real642str_1D_4() result(assertion)
        use, intrinsic :: iso_fortran_env, only: RK => real64, IK => int32
        implicit none
        logical                     :: assertion
        integer     , parameter     :: nthis = 2
        real(RK)    , parameter     :: this_ref(nthis) = [1.23456798e-30_RK, 2.32456798e+30_RK]
        real(RK)                    :: this(nthis)
        character(:), allocatable   :: string
        assertion = .true.
        string = num2str(this_ref,"(*(g0,:,' '))",263_IK)
        read(string,*) this
        assertion = all(this == this_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_real642str_1D_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_real642str_2D_1() result(assertion)
        use, intrinsic :: iso_fortran_env, only: RK => real64
        implicit none
        logical                     :: assertion
        integer     , parameter     :: nthis = 2
        real(RK)    , parameter     :: this_ref(nthis,nthis) = reshape([1.23456798e-30_RK, 2.32456798e+30_RK, 1.23456798e-30_RK, 2.32456798e+30_RK], shape=shape(this_ref))
        real(RK)                    :: this(nthis,nthis)
        character(:), allocatable   :: string
        assertion = .true.
        string = num2str(this_ref)
        read(string,*) this
        assertion = all(this == this_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_real642str_2D_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_real642str_2D_2() result(assertion)
        use, intrinsic :: iso_fortran_env, only: RK => real64
        implicit none
        logical                     :: assertion
        integer     , parameter     :: nthis = 2
        real(RK)    , parameter     :: this_ref(nthis,nthis) = reshape([1.23456798e-30_RK, 2.32456798e+30_RK, 1.23456798e-30_RK, 2.32456798e+30_RK], shape=shape(this_ref))
        real(RK)                    :: this(nthis,nthis)
        character(:), allocatable   :: string
        assertion = .true.
        string = num2str(this_ref,"(*(g0,:,' '))")
        read(string,*) this
        assertion = all(this == this_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_real642str_2D_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_real642str_2D_3() result(assertion)
        use, intrinsic :: iso_fortran_env, only: RK => real64, IK => int32
        implicit none
        logical                     :: assertion
        integer     , parameter     :: nthis = 2
        real(RK)    , parameter     :: this_ref(nthis,nthis) = reshape([1.23456798e-30_RK, 2.32456798e+30_RK, 1.23456798e-30_RK, 2.32456798e+30_RK], shape=shape(this_ref))
        real(RK)                    :: this(nthis,nthis)
        character(:), allocatable   :: string
        assertion = .true.
        string = num2str(this_ref,"(*(g0,:,' '))",128_IK)
        read(string,*) this
        assertion = all(this == this_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_real642str_2D_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> minLen can be larger than the length of the constructed string.
    function test_real642str_2D_4() result(assertion)
        use, intrinsic :: iso_fortran_env, only: RK => real64, IK => int32
        implicit none
        logical                     :: assertion
        integer     , parameter     :: nthis = 2
        real(RK)    , parameter     :: this_ref(nthis,nthis) = reshape([1.23456798e-30_RK, 2.32456798e+30_RK, 1.23456798e-30_RK, 2.32456798e+30_RK], shape=shape(this_ref))
        real(RK)                    :: this(nthis,nthis)
        character(:), allocatable   :: string
        assertion = .true.
        string = num2str(this_ref,"(*(g0,:,' '))",256_IK)
        read(string,*) this
        assertion = all(this == this_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "this_ref  = ", this_ref
            write(Test%outputUnit,"(*(g0))") "this      = ", this
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_real642str_2D_4

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_log2str_1() result(assertion)
        implicit none
        logical                     :: assertion
        character(:), allocatable   :: strout
        logical                     :: logic
        logic = .true.
        strout = log2str(logic)
        assertion = strout == "TRUE"
    end function test_log2str_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_log2str_2() result(assertion)
        implicit none
        logical                     :: assertion
        character(:), allocatable   :: strout
        logical                     :: logic
        logic = .false.
        strout = log2str(logic)
        assertion = strout == "FALSE"
    end function test_log2str_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_IntStr_type_1() result(assertion)
        use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
        implicit none
        logical                     :: assertion
        type(IntStr_type)           :: IntStr
        IntStr%str = IntStr%int2str(123_int32)
        assertion = IntStr%str == "123"
    end function test_IntStr_type_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_IntStr_type_2() result(assertion)
        use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
        implicit none
        logical                     :: assertion
        type(IntStr_type)           :: IntStr
        IntStr%str = IntStr%int2str(123_int64)
        assertion = IntStr%str == "123"
    end function test_IntStr_type_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_num2str_1() result(assertion)
        use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
        implicit none
        logical                     :: assertion
        type(RealStr_type)          :: RealStr
        RealStr%str = RealStr%real2str(123._real32,"(F10.4)",15)
        assertion = RealStr%str == "123.0000       "
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "Number=123_real32"
            write(Test%outputUnit,"(*(g0))") "RealStr%str = '", RealStr%str, "'"
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_num2str_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_num2str_2() result(assertion)
        use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
        implicit none
        logical                     :: assertion
        type(RealStr_type)          :: RealStr
        RealStr%str = RealStr%real2str(123._real64,"(F10.4)",15)
        assertion = RealStr%str == "123.0000       "
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "Number=123_real64"
            write(Test%outputUnit,"(*(g0))") "RealStr%str = '", RealStr%str, "'"
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_num2str_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_padString_1() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        integer(IK) , parameter     :: paddedLen = 30_IK
        character(*), parameter     :: string_nonPadded = "ParaMonte"
        character(*), parameter     :: stringPadded_ref = "ParaMonte....................."
        character(*), parameter     :: symbol = "."
        character(:), allocatable   :: stringPadded
        stringPadded = padString(string_nonPadded, symbol, paddedLen)
        assertion = stringPadded == stringPadded_ref .and. len(stringPadded) == len(stringPadded_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "string_nonPadded  : '", string_nonPadded, "'"
            write(Test%outputUnit,"(*(g0))") "stringPadded_ref  : '", stringPadded_ref, "'"
            write(Test%outputUnit,"(*(g0))") "stringPadded      : '", stringPadded, "'"
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_padString_1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_padString_2() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        integer(IK) , parameter     :: paddedLen = 9_IK
        character(*), parameter     :: string_nonPadded = "ParaMonte"
        character(*), parameter     :: stringPadded_ref = "ParaMonte"
        character(*), parameter     :: symbol = "."
        character(:), allocatable   :: stringPadded
        stringPadded = padString(string_nonPadded, symbol, paddedLen)
        assertion = stringPadded == stringPadded_ref .and. len(stringPadded) == len(stringPadded_ref)
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "string_nonPadded  : '", string_nonPadded, "'"
            write(Test%outputUnit,"(*(g0))") "stringPadded_ref  : '", stringPadded_ref, "'"
            write(Test%outputUnit,"(*(g0))") "stringPadded      : '", stringPadded, "'"
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_padString_2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> When `len(string) > paddedLen`, the full string must be returned without any padding.
    function test_padString_3() result(assertion)
        use Constants_mod, only: IK
        implicit none
        logical                     :: assertion
        integer     , parameter     :: paddedLen = 5_IK
        character(*), parameter     :: string_nonPadded = "ParaMonte"
        character(*), parameter     :: stringPadded_ref = "ParaM"
        character(*), parameter     :: symbol = "."
        character(:), allocatable   :: stringPadded
        stringPadded = padString(string_nonPadded, symbol, paddedLen)
        assertion = stringPadded == stringPadded_ref .and. len(stringPadded) == paddedLen .and. stringPadded == stringPadded_ref
        if (Test%isDebugMode .and. .not. assertion) then
        ! LCOV_EXCL_START
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "string_nonPadded  : '", string_nonPadded, "'"
            write(Test%outputUnit,"(*(g0))") "stringPadded_ref  : '", stringPadded_ref, "'"
            write(Test%outputUnit,"(*(g0))") "stringPadded      : '", stringPadded, "'"
            write(Test%outputUnit,"(*(g0))")
        end if
        ! LCOV_EXCL_STOP
    end function test_padString_3

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Test_String_mod ! LCOV_EXCL_LINE