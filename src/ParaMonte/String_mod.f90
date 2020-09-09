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

module String_mod

    use, intrinsic :: iso_fortran_env, only: int8
    use Constants_mod, only: IK, RK
    use JaggedArray_mod, only: CharVec_type
    implicit none

    public

    character(*), parameter :: MODULE_NAME = "@String_mod"

    integer(int8) :: NUM2STR_MAXLEN = 63_int8

    type :: IntStr_type
        integer(IK)                 :: val
        character(:), allocatable   :: str
    contains
        procedure, nopass           :: int322str, int642str
        generic                     :: int2str => int322str, int642str
    end type IntStr_type

    type :: RealStr_type
        integer(IK)                 :: val
        character(:), allocatable   :: str
    contains
        procedure, nopass           :: real322str, real642str
        generic                     :: real2str => real322str, real642str
    end type RealStr_type

    type :: String_type
        character(:)      , allocatable   :: value
        type(CharVec_type), allocatable   :: Parts(:)
        integer(IK)                       :: nPart = 0_IK
    contains
        procedure, nopass :: replaceStr, splitStr, getLowerCase, getUpperCase, isInteger, isDigit
        procedure, nopass :: str2int, str2real, str2int32, str2int64, str2real32, str2real64
    end type String_type

    interface int2str
        module procedure int322str, int642str
    end interface

    interface real2str
        module procedure real322str, real642str, real642str_1D, real642str_2D
    end interface

    interface num2str
        module procedure int322str, int642str, real642str, real642str_1D, real642str_2D, real322str, log2str
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure recursive function replaceStr(string,search,substitute) result(modifiedString)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: replaceStr
#endif
        implicit none
        character(len=*), intent(in)  :: string, search, substitute
        character(len=:), allocatable :: modifiedString
        integer(IK)                   :: i, stringLen, searchLen
        stringLen = len(string)
        searchLen = len(search)
        if (stringLen==0 .or. searchLen==0) then
            modifiedString = ""
            return
        elseif (stringLen<searchLen) then
            modifiedString = string
            return
        end if
        i = 1
        do
            if (string(i:i+searchLen-1)==search) then
                modifiedString = string(1:i-1) // substitute // replaceStr(string(i+searchLen:stringLen),search,substitute)
                exit
            end if
            if (i+searchLen>stringLen) then
                modifiedString = string
                exit
            end if
            i = i + 1
            cycle
        end do
    end function replaceStr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function splitStr(string,delimiter,nPart) result(Parts)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: splitStr
#endif

        implicit none
        character(len=*)    , intent(in)            :: string, delimiter
        integer(IK)         , intent(out), optional :: nPart
        character(len=:)    , allocatable           :: dummyStr
        type(CharVec_type)  , allocatable           :: Parts(:)
        integer(IK)                                 :: maxNumSplit
        integer(IK)                                 :: stringLen, delimLen, splitCounter, currentPos

        dummyStr  = string
        delimLen  = len(delimiter)
        stringLen = len(dummyStr)

        if (delimLen==0) then
            allocate(Parts(1))
            Parts(1)%record = string
            return
        end if

        maxNumSplit = 1 + stringLen / delimLen
        allocate(Parts(maxNumSplit))
        splitCounter = 1
        loopParseString: do
            if (stringLen<delimLen) then
                Parts(splitCounter)%record = dummyStr
                exit loopParseString
            elseif (stringLen==delimLen) then
                ! Note that in Fortran: 'amir '=='amir' = .true.
                ! https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/275823)
                if (dummyStr==delimiter) then
                    Parts(splitCounter)%record = ""
                else
                    Parts(splitCounter)%record = dummyStr
                end if
                exit loopParseString
            elseif (dummyStr(1:delimLen)==delimiter) then
                dummyStr = dummyStr(delimLen+1:stringLen)
                stringLen = len(dummyStr)
                cycle loopParseString
            else
                currentPos = 2
                loopSearchString: do
                    if (dummyStr(currentPos:currentPos+delimLen-1)==delimiter) then
                        Parts(splitCounter)%record = dummyStr(1:currentPos-1)
                        if (currentPos+delimLen>stringLen) then
                            exit loopParseString
                        else
                            splitCounter = splitCounter + 1
                            dummyStr = dummyStr(currentPos+delimLen:stringLen)
                            stringLen = len(dummyStr)
                            cycle loopParseString
                        end if
                    else
                        currentPos = currentPos + 1
                        if (stringLen<currentPos+delimLen-1) then
                            Parts(splitCounter)%record = dummyStr
                            exit loopParseString
                        end if
                        cycle loopSearchString
                    end if
                end do loopSearchString
            end if
        end do loopParseString
        Parts = Parts(1:splitCounter)
        if (present(nPart)) nPart = splitCounter

    end function splitStr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function isDigit(singleChar) result(stringIsDigit)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: isDigit
#endif
        character(1), intent(in)    :: singleChar
        logical                     :: stringIsDigit
        character(*), parameter     :: Digit(10) = ["0","1","2","3","4","5","6","7","8","9"]
        integer                     :: j
        stringIsDigit = .false.
        loopOverDigit: do j = 1,10
            if (singleChar==Digit(j)) then
                stringIsDigit = .true.
                exit loopOverDigit
            end if
        end do loopOverDigit
    end function isDigit

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function isInteger(string) result(stringIsInteger)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: isInteger
#endif
        character(*), intent(in)    :: string
        logical                     :: stringIsInteger
        character(*), parameter     :: Digit(10) = ["0","1","2","3","4","5","6","7","8","9"]
        integer                     :: i, j
        do i = 1,len(string)
            stringIsInteger = .false.
            loopOverDigit: do j = 1,10
                if (string(i:i)==Digit(j)) then
                    stringIsInteger = .true.
                    exit loopOverDigit
                end if
            end do loopOverDigit
            if (.not.stringIsInteger) return
        end do
    end function isInteger

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLowerCase(string) result(output)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLowerCase
#endif
        character(*), intent(in) :: string
        integer(IK), parameter   :: DUC = ichar('A') - ichar('a')
        character(len(string))   :: output
        character                :: ch
        integer(IK)              :: i
        do i = 1,len(string)
            ch = string(i:i)
            if (ch>='A' .and. ch<='Z') ch = char(ichar(ch)-DUC)
            output(i:i) = ch
        end do
    end function getLowerCase

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getUpperCase(string) result(output)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getUpperCase
#endif
        character(*), intent(in) :: string
        integer(IK), parameter   :: DUC = ichar('A') - ichar('a')
        character(len(string))   :: output
        character                :: ch
        integer(IK)              :: i
        do i = 1,len(string)
            ch = string(i:i)
            if (ch>='a' .and. ch<='z') ch = char(ichar(ch)+DUC)
            output(i:i) = ch
        end do
    end function getUpperCase

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ATTN: legacy code, do not use. The substitute functions are >1 order of magnitude faster. Changes a string to lower case
    pure function getLowerCaseOld(string)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLowerCaseOld
#endif
        implicit None
        character(*), intent(in) :: string
        character(len(string))   :: getLowerCaseOld
        character(26), parameter :: lowerCase = 'abcdefghijklmnopqrstuvwxyz', upperCase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        integer(IK)              :: ic, i
        getLowerCaseOld = string
        do i = 1, len(string)
            ic = INDEX(upperCase, string(i:i))
            if (ic > 0) getLowerCaseOld(i:i) = lowerCase(ic:ic)
        end do
    end function getLowerCaseOld

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ATTN: legacy code, do not use. The substitute functions are >1 order of magnitude faster. Changes a string to upper case
    pure function getUpperCaseOld(string)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getUpperCaseOld
#endif
        implicit None
        character(*), intent(in) :: string
        character(len(string))   :: getUpperCaseOld
        character(26), parameter :: lowerCase = 'abcdefghijklmnopqrstuvwxyz', upperCase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        integer(IK)              :: ic, i
        getUpperCaseOld = string
        do i = 1, len(string)
            ic = INDEX(lowerCase, string(i:i))
            if (ic > 0) getUpperCaseOld(i:i) = upperCase(ic:ic)
        end do
    end function getUpperCaseOld

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function log2str(logicalIn)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: log2str
#endif
        implicit none
        logical     , intent(in)           :: logicalIn
        character(:), allocatable          :: log2str
        allocate(character(NUM2STR_MAXLEN) :: log2str)
        if (logicalIn) then
            log2str = "TRUE"
        else
            log2str = "FALSE"
        end if
    end function log2str

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function int322str(integerIn,formatIn,minLen)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: int322str
#endif
        use, intrinsic :: iso_fortran_env, only: int32
        implicit none
        integer(int32) , intent(in)           :: integerIn
        character(*)   , intent(in), optional :: formatIn
        integer(IK)    , intent(in), optional :: minLen
        character(:)   , allocatable          :: int322str
        allocate(character(NUM2STR_MAXLEN) :: int322str)
        if (present(formatIn)) then
            write(int322str,formatIn) integerIn
        else
            write(int322str,*) integerIn
        end if
        if (present(minLen)) then
            int322str = adjustl(int322str)
            int322str = int322str(1:minLen)
        else
            int322str = trim(adjustl(int322str))
        end if
    end function int322str

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function int642str(integerIn,formatIn,minLen)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: int642str
#endif
        use, intrinsic :: iso_fortran_env, only: int64
        implicit none
        integer(int64), intent(in)           :: integerIn
        character(*)  , intent(in), optional :: formatIn
        integer(IK)   , intent(in), optional :: minLen
        character(:)  , allocatable          :: int642str
        allocate(character(NUM2STR_MAXLEN) :: int642str)
        if (present(formatIn)) then
            write(int642str,formatIn) integerIn
        else
            write(int642str,*) integerIn
        end if
        if (present(minLen)) then
            int642str = adjustl(int642str)
            int642str = int642str(1:minLen)
        else
            int642str = trim(adjustl(int642str))
        end if
    end function int642str

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function real322str(realIn,formatIn,minLen)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: real322str
#endif
        use, intrinsic :: iso_fortran_env, only: real32
        implicit none
        real(real32), intent(in)           :: realIn
        character(*), intent(in), optional :: formatIn
        integer(IK) , intent(in), optional :: minLen
        character(:), allocatable          :: real322str
        allocate(character(NUM2STR_MAXLEN) :: real322str)
        if (present(formatIn)) then
            write(real322str,formatIn) realIn
        else
            write(real322str,*) realIn
        end if
        if (present(minLen)) then
            real322str = adjustl(real322str)
            real322str = real322str(1:minLen)
        else
            real322str = trim(adjustl(real322str))
        end if
    end function real322str

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function real642str(realIn,formatIn,minLen)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: real642str
#endif
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none
        real(real64), intent(in)           :: realIn
        character(*), intent(in), optional :: formatIn
        integer(IK) , intent(in), optional :: minLen
        character(:), allocatable          :: real642str
        allocate(character(NUM2STR_MAXLEN) :: real642str)
        if (present(formatIn)) then
            write(real642str,formatIn) realIn
        else
            write(real642str,*) realIn
        end if
        if (present(minLen)) then
            real642str = adjustl(real642str)
            real642str = real642str(1:minLen)
        else
            real642str = trim(adjustl(real642str))
        end if
    end function real642str

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function real642str_1D(RealIn,formatIn,minLen)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: real642str_1D
#endif
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none
        real(real64), intent(in)           :: RealIn(:)
        character(*), intent(in), optional :: formatIn
        integer(IK) , intent(in), optional :: minLen
        character(:), allocatable          :: real642str_1D
        allocate(character(NUM2STR_MAXLEN*size(RealIn)) :: real642str_1D)
        if (present(formatIn)) then
            write(real642str_1D,formatIn) RealIn
        else
            write(real642str_1D,"(*(g0.15,:,','))") RealIn
        end if
        if (present(minLen)) then
            real642str_1D = adjustl(real642str_1D)
            real642str_1D = real642str_1D(1:minLen)
        else
            real642str_1D = trim(adjustl(real642str_1D))
        end if
    end function real642str_1D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function real642str_2D(RealIn,formatIn,minLen)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: real642str_2D
#endif
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none
        real(real64), intent(in)           :: RealIn(:,:)
        character(*), intent(in), optional :: formatIn
        integer(IK) , intent(in), optional :: minLen
        character(:), allocatable          :: real642str_2D
        allocate(character(NUM2STR_MAXLEN*size(RealIn,1)*size(RealIn,2)) :: real642str_2D)
        if (present(formatIn)) then
            write(real642str_2D,formatIn) RealIn
        else
            write(real642str_2D,"(*(g0.15,:,','))") RealIn
        end if
        if (present(minLen)) then
            real642str_2D = adjustl(real642str_2D)
            real642str_2D = real642str_2D(1:minLen)
        else
            real642str_2D = trim(adjustl(real642str_2D))
        end if
    end function real642str_2D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function str2int(str,iostat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: str2int
#endif
        use Constants_mod, only: IK
        implicit none
        character(*), intent(in)        :: str
        integer, optional, intent(out)  :: iostat
        integer(IK)                     :: str2int
        if (present(iostat)) then
            iostat = 0
            read(str,*,iostat=iostat) str2int
        else
            read(str,*) str2int
        endif
    end function str2int

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function str2int32(str,iostat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: str2int32
#endif
        use, intrinsic :: iso_fortran_env, only: int32
        implicit none
        character(*), intent(in)        :: str
        integer, optional, intent(out)  :: iostat
        integer(int32)                  :: str2int32
        if (present(iostat)) then
            iostat = 0
            read(str,*,iostat=iostat) str2int32
        else
            read(str,*) str2int32
        endif
    end function str2int32

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function str2int64(str,iostat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: str2int64
#endif
        use, intrinsic :: iso_fortran_env, only: int64
        implicit none
        character(*), intent(in)        :: str
        integer, optional, intent(out)  :: iostat
        integer(int64)                  :: str2int64
        if (present(iostat)) then
            iostat = 0
            read(str,*,iostat=iostat) str2int64
        else
            read(str,*) str2int64
        endif
    end function str2int64

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function str2real(str,iostat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: str2real
#endif
        use Constants_mod, only: RK
        implicit none
        character(*), intent(in)        :: str
        integer, optional, intent(out)  :: iostat
        real(RK)                        :: str2real
        if (present(iostat)) then
            iostat = 0
            read(str,*,iostat=iostat) str2real
        else
            read(str,*) str2real
        endif
    end function str2real

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function str2real32(str,iostat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: str2real32
#endif
        use, intrinsic :: iso_fortran_env, only: real32
        implicit none
        character(*), intent(in)        :: str
        integer, optional, intent(out)  :: iostat
        real(real32)                    :: str2real32
        if (present(iostat)) then
            iostat = 0
            read(str,*,iostat=iostat) str2real32
        else
            read(str,*) str2real32
        endif
    end function str2real32

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function str2real64(str,iostat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: str2real64
#endif
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none
        character(*), intent(in)        :: str
        integer, optional, intent(out)  :: iostat
        real(real64)                    :: str2real64
        if (present(iostat)) then
            iostat = 0
            read(str,*,iostat=iostat) str2real64
        else
            read(str,*) str2real64
        endif
    end function str2real64

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module String_mod
