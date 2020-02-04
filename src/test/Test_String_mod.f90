!**********************************************************************************************************************************
!**********************************************************************************************************************************
!
!  ParaMonte: plain powerful parallel Monte Carlo library.
!
!  Copyright (C) 2012-present, The Computational Data Science Lab
!
!  This file is part of ParaMonte library. 
!
!  ParaMonte is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, version 3 of the License.
!
!  ParaMonte is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************

module Test_String_mod

    use String_mod
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_String

    type(Test_type) :: Test

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_String()
        implicit none
        Test = Test_type(moduleName=MODULE_NAME)
        call test_String_type()
        call test_replaceStr()
        call test_getLowerCase()
        call test_getUpperCase()
        call test_log2str()
        call test_num2str()
        call Test%finalize()
#ifdef CAF_ENABLED
        sync all
#endif

    end subroutine test_String

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_String_type()

        implicit none
        type(String_type)   :: String
        integer             :: i

        if (Test%Image%isFirst) call Test%testing("String_type%splitStr()")

        String%value = "  StringString "
        String%Parts  = String%splitStr(string=String%value,delimiter='String')
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "String%value = '", String%value, "'"
            write(Test%outputUnit,"(*(g0))") "String%splitStr(string=String%value,delimiter='String') = ", &
                                               (String%Parts(i)%record,i=1,size(String%Parts))
            write(Test%outputUnit,"(*(g0))")
        end if
        Test%assertion = String%Parts(1)%record == "  " .and. String%Parts(2)%record == " "
        call Test%verify()

        if (Test%Image%isFirst) call Test%testing("String_type%splitStr()")

        String%value = "  Stringstring "
        String%Parts  = String%splitStr(string=String%value,delimiter='str')
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "String%value = '", String%value, "'"
            write(Test%outputUnit,"(*(g0))") "String%splitStr(string=String%value,delimiter='str') = ", &
                                               (String%Parts(i)%record,i=1,size(String%Parts))
            write(Test%outputUnit,"(*(g0))")
        end if
        Test%assertion = String%Parts(1)%record == "  String" .and. String%Parts(2)%record == "ing "
        call Test%verify()

        if (Test%Image%isFirst) call Test%testing("String_type%splitStr()")

        String%value = "   "
        String%Parts  = String%splitStr(string=String%value,delimiter=' ')
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "String%value = '", String%value, "'"
            write(Test%outputUnit,"(*(g0))") "String%splitStr(string=String%value,delimiter=' ') = ", &
                                               ("'"//String%Parts(i)%record//"'",i=1,size(String%Parts))
            write(Test%outputUnit,"(A,*(g0))")
        end if
        Test%assertion = String%Parts(1)%record == ""
        call Test%verify()

        if (Test%Image%isFirst) call Test%testing("String%str2int()")

        Test%assertion = 123456_IK == str2int("123456")
        call Test%verify()

    end subroutine test_String_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_replaceStr()

        implicit none
        type(String_type)   :: String
        character(:), allocatable :: strout

        if (Test%Image%isFirst) call Test%testing("String_type%replaceStr()")

        String%value = "  StringString "
        strout = String%replaceStr(string=String%value, search="Str", substitute="")
        Test%assertion = strout == "  inging "
        call Test%verify()

        String%value = ""
        strout = String%replaceStr(string=String%value, search="", substitute="")
        Test%assertion = strout == ""
        call Test%verify()

    end subroutine test_replaceStr

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getLowerCase()

        implicit none
        type(String_type)   :: String
        character(:), allocatable :: strout

        if (Test%Image%isFirst) call Test%testing("String_type%getLowerCase()")

        String%value = "  StringString !@#$"
        strout = String%getLowerCase(string=String%value)
        Test%assertion = strout == "  stringstring !@#$"
        call Test%verify()

        String%value = ""
        strout = String%getLowerCase(string=String%value)
        Test%assertion = strout == ""
        call Test%verify()

    end subroutine test_getLowerCase

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getUpperCase()

        implicit none
        type(String_type)   :: String
        character(:), allocatable :: strout

        if (Test%Image%isFirst) call Test%testing("String_type%getUpperCase()")

        String%value = "  StringString !@#$"
        strout = String%getUpperCase(string=String%value)
        Test%assertion = strout == "  STRINGSTRING !@#$"
        call Test%verify()

        String%value = ""
        strout = String%getUpperCase(string=String%value)
        Test%assertion = strout == ""
        call Test%verify()

    end subroutine test_getUpperCase

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_log2str()

        implicit none
        character(:), allocatable   :: strout
        logical                     :: logic

        if (Test%Image%isFirst) call Test%testing("log2str()")

        logic = .true.
        strout = log2str(logic)
        Test%assertion = strout == "TRUE"
        call Test%verify()

        logic = .false.
        strout = log2str(logic)
        Test%assertion = strout == "FALSE"
        call Test%verify()

    end subroutine test_log2str

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_IntStr_type()

        use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
        implicit none
        type(IntStr_type)           :: IntStr

        if (Test%Image%isFirst) call Test%testing("IntStr_type%int2str()")

        IntStr%str = IntStr%int2str(123_int32)
        Test%assertion = IntStr%str == "123"
        call Test%verify()

        if (Test%Image%isFirst) call Test%testing("IntStr%int642str()")

        IntStr%str = IntStr%int2str(123_int64)
        Test%assertion = IntStr%str == "123"
        call Test%verify()

    end subroutine test_IntStr_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_num2str()

        use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
        implicit none
        type(RealStr_type)          :: RealStr

        if (Test%Image%isFirst) call Test%testing("RealStr%real322str()")

        RealStr%str = RealStr%real2str(123._real32,"(F10.4)",15)
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "Number=123_real32"
            write(Test%outputUnit,"(*(g0))") "RealStr%str = '", RealStr%str, "'"
            write(Test%outputUnit,"(*(g0))")
        end if
        Test%assertion = RealStr%str == "123.0000       "
        call Test%verify()

        if (Test%Image%isFirst) call Test%testing("RealStr%real642str()")

        RealStr%str = RealStr%real2str(123._real64,"(F10.4)",15)
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "Number=123_real64"
            write(Test%outputUnit,"(*(g0))") "RealStr%str = '", RealStr%str, "'"
            write(Test%outputUnit,"(*(g0))")
        end if
        Test%assertion = RealStr%str == "123.0000       "
        call Test%verify()

    end subroutine test_num2str

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_String_mod
