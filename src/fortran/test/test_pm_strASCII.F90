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
!>  This module contains tests of the module [pm_strASCII](@ref pm_strASCII).
!>
!>  \fintest
!>
!>  \author
!>  \AmirShahmoradi

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module test_pm_strASCII

    use pm_strASCII
    use pm_kind, only: LK
    use pm_container, only: strc => css_pdt
    use pm_test, only: test_type, LK
    implicit none

    private
    public :: setTest
    type(test_type) :: test

    type(strc), allocatable :: NonIntSet(:)
    type(strc), allocatable :: IntSet(:)
    type(strc), allocatable :: RealSet(:)
    type(strc), allocatable :: NonRealSet(:)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getLocSpace_1() result(assertion)
        logical(LK) :: assertion
        assertion = .true._LK
        assertion = assertion .and. getLocSpace(SK_"") == 0_IK
        call test%assert(assertion, MODULE_NAME//SK_': The condition `getLocSpace(SK_"") == 0_IK` must hold.')
        assertion = assertion .and. getLocSpace(SK_" ") == 1_IK
        call test%assert(assertion, MODULE_NAME//SK_': The condition `getLocSpace(SK_" ") == 1_IK` must hold.')
        assertion = assertion .and. getLocSpace(SK_" paramonte") == 1_IK
        call test%assert(assertion, MODULE_NAME//SK_': The condition `getLocSpace(SK_" paramonte") == 1_IK` must hold.')
        assertion = assertion .and. getLocSpace(SK_"paramonte ") == 10_IK
        call test%assert(assertion, MODULE_NAME//SK_': The condition `getLocSpace(SK_"paramonte ") == 10_IK` must hold.')
        assertion = assertion .and. getLocSpace(SK_"paramonte") == 0_IK
        call test%assert(assertion, MODULE_NAME//SK_': The condition `getLocSpace(SK_"paramonte") == 0_IK` must hold.')
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isCharDigit_1() result(assertion)
        use pm_val2str, only: getStr
        implicit none
        logical(LK)                 :: assertion
        integer(IK)                 :: i
        assertion = .true._LK
        do i = 1, 127
            if (index(SK_"0123456789", achar(i, kind = SK)) > 0) then
                assertion = assertion .and. isCharDigit(achar(i))
            else
                assertion = assertion .and. .not. isCharDigit(achar(i))
            end if
            call test%assert(assertion, MODULE_NAME//SK_"isCharDigit() must correctly identify "//achar(i)//SK_" as a digit or a non-digit.")
        end do
        assertion = assertion .and. all(isCharDigit([SK_"0", SK_" ", SK_"+", SK_"."]) .eqv. [.true._LK, .false._LK, .false._LK, .false._LK])
        call test%assert(assertion, MODULE_NAME//SK_'isCharDigit() must correctly identify [SK_"0", SK_" ", SK_"+", SK_"."] as digits or non-digits.')
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isStrDigitAll_1() result(assertion)
        use pm_val2str, only: getStr
        implicit none
        logical(LK)                 :: assertion
        integer(IK)                 :: i
        assertion = .true._LK
        do i = 0, 9
            assertion = assertion .and. isStrDigitAll(SK_"0")
            call test%assert(assertion, MODULE_NAME//SK_"isStrDigitAll() must correctly identify "//getStr(i)//SK_" as numeric.")
        end do
        assertion = assertion .and. isStrDigitAll(SK_"0123")
        call test%assert(assertion, MODULE_NAME//SK_"@isStrDigitAll() must correctly identify `0123` as numeric.")
        assertion = assertion .and. .not. isStrDigitAll(SK_"123.")
        call test%assert(assertion, MODULE_NAME//SK_"@isStrDigitAll() must correctly identify `123.` as not numeric.")
        assertion = assertion .and. .not. isStrDigitAll(SK_"123 ")
        call test%assert(assertion, MODULE_NAME//SK_"@isStrDigitAll() must correctly identify `123 ` as not numeric.")
        assertion = assertion .and. .not. isStrDigitAll(SK_" ")
        call test%assert(assertion, MODULE_NAME//SK_"@isStrDigitAll() must correctly identify ` ` as not numeric.")
        assertion = assertion .and. .not. isStrDigitAll(SK_"-123")
        call test%assert(assertion, MODULE_NAME//SK_"@isStrDigitAll() must correctly identify `-123` as not numeric.")
        assertion = assertion .and. .not. isStrDigitAll(SK_"")
        call test%assert(assertion, MODULE_NAME//SK_"@isStrDigitAll() must correctly identify `` as not numeric.")
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isStrReal_1() result(assertion)
        use pm_val2str, only: getStr
        implicit none
        logical(LK)                 :: assertion
        integer(IK)                 :: i
        assertion = .true._LK
        do i = 1, size(IntSet)
            assertion = assertion .and. isStrReal(IntSet(i)%val)
            call test%assert(assertion, MODULE_NAME//SK_"isStrReal() must correctly identify """//IntSet(i)%val//SK_""" as a real number.")
        end do
        do i = 1, size(RealSet)
            assertion = assertion .and. isStrReal(RealSet(i)%val)
            call test%assert(assertion, MODULE_NAME//SK_"isStrReal() must correctly identify """//RealSet(i)%val//SK_""" as a real number.")
        end do
        do i = 1, size(NonRealSet)
            assertion = assertion .and. .not. isStrReal(NonRealSet(i)%val)
            call test%assert(assertion, MODULE_NAME//SK_"isStrReal() must correctly identify """//NonRealSet(i)%val//SK_""" as not a real number.")
        end do
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isStrInteger_1() result(assertion)
        use pm_val2str, only: getStr
        implicit none
        logical(LK)                 :: assertion
        integer(IK)                 :: i
        assertion = .true._LK
        do i = 1, size(IntSet)
            assertion = assertion .and. isStrInteger(IntSet(i)%val)
            call test%assert(assertion, MODULE_NAME//SK_"isStrInteger() must correctly identify """//IntSet(i)%val//SK_""" as an integer number.")
        end do
        do i = 1, size(NonIntSet)
            assertion = assertion .and. .not. isStrInteger(NonIntSet(i)%val)
            call test%assert(assertion, MODULE_NAME//SK_"isStrInteger() must correctly identify """//NonIntSet(i)%val//SK_""" as not an integer number.")
        end do
    end function


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isStrNumber_1() result(assertion)
        use pm_val2str, only: getStr
        implicit none
        logical(LK)                 :: assertion
        integer(IK)                 :: i
        assertion = .true._LK
        do i = 0, 9
            assertion = assertion .and. isStrNumber(SK_"0")
            call test%assert(assertion, MODULE_NAME//SK_"isStrNumber() must correctly identify "//getStr(i)//SK_" as numeric.")
        end do
        assertion = assertion .and. isStrNumber(SK_"0123")
        call test%assert(assertion, MODULE_NAME//SK_"@isStrNumber() must correctly identify `0123` as a number.")
        assertion = assertion .and. isStrNumber(SK_"123.")
        call test%assert(assertion, MODULE_NAME//SK_"@isStrNumber() must correctly identify `123.` as a number.")
        assertion = assertion .and. isStrNumber(SK_"123 ")
        call test%assert(assertion, MODULE_NAME//SK_"@isStrNumber() must correctly identify `123 ` as a number.")
        assertion = assertion .and. isStrNumber(SK_" 13 ")
        call test%assert(assertion, MODULE_NAME//SK_"@isStrNumber() must correctly identify ` 13 ` as a number.")
        assertion = assertion .and. isStrNumber(SK_"-13 ")
        call test%assert(assertion, MODULE_NAME//SK_"@isStrNumber() must correctly identify `-13 ` as a number.")
        assertion = assertion .and. isStrNumber(SK_"-13.")
        call test%assert(assertion, MODULE_NAME//SK_"@isStrNumber() must correctly identify `-13.` as a number.")
        assertion = assertion .and. isStrNumber(SK_"+13.")
        call test%assert(assertion, MODULE_NAME//SK_"@isStrNumber() must correctly identify `+13.` as a number.")
        assertion = assertion .and. isStrNumber(SK_" +3 ")
        call test%assert(assertion, MODULE_NAME//SK_"@isStrNumber() must correctly identify ` +3 ` as a number.")
        assertion = assertion .and. isStrNumber(SK_" +3 ")
        call test%assert(assertion, MODULE_NAME//SK_"@isStrNumber() must correctly identify ` +3 ` as a number.")
        assertion = assertion .and. .not. isStrNumber(SK_" ")
        call test%assert(assertion, MODULE_NAME//SK_"@isStrNumber() must correctly identify ` ` as not numeric.")
        assertion = assertion .and. isStrNumber(SK_"-123")
        call test%assert(assertion, MODULE_NAME//SK_"@isStrNumber() must correctly identify `-123` as not numeric.")
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isCharLower_1() result(assertion)
        logical(LK) :: assertion
        integer(IK) :: i
        assertion = .true._LK
        do i = 1_IK, 127_IK
            assertion = assertion .and. (isCharLower(char(i, kind = SK)) .eqv. any(ALPHA_LOWER_VEC_SK == char(i, kind = SK)))
            call test%assert(assertion, MODULE_NAME//SK_': isCharLower() must properly recognize lower-case characters.')
        end do
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isCharUpper_1() result(assertion)
        logical(LK) :: assertion
        integer(IK) :: i
        assertion = .true._LK
        do i = 1_IK, 127_IK
            assertion = assertion .and. (isCharUpper(char(i, kind = SK)) .eqv. any(ALPHA_UPPER_VEC_SK == char(i, kind = SK)))
            call test%assert(assertion, MODULE_NAME//SK_': isCharUpper() must properly recognize upper-case characters.')
        end do
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isStrLowerAll_1() result(assertion)
        logical(LK) :: assertion
        integer(IK) :: i
        assertion = .true._LK
        assertion = assertion .and. .not. isStrLowerAll(SK_"")
        call test%assert(assertion, MODULE_NAME//SK_': The condition `.not. isStrLowerAll(SK_"")` must hold.')
        do i = 1_IK, 127_IK
            assertion = assertion .and. (isStrLowerAll(char(i, kind = SK)) .eqv. any(ALPHA_LOWER_VEC_SK == char(i, kind = SK)))
            call test%assert(assertion, MODULE_NAME//SK_': isStrLowerAll() must properly recognize upper-case characters.')
        end do
        assertion = assertion .and. isStrLowerAll(SK_"paramonte")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrLowerAll(SK_"paramonte")` must hold.')
        assertion = assertion .and. .not. isStrLowerAll(SK_"PARAMONTE ")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrLowerAll(SK_"PARAMONTE ")` must not hold.')
        assertion = assertion .and. .not. isStrLowerAll(SK_"paramonte ")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrLowerAll(SK_"paramonte ")` must not hold.')
        assertion = assertion .and. .not. isStrLowerAll(SK_"paraMONTE")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrLowerAll(SK_"paraMONTE")` must not hold.')
        assertion = assertion .and. .not. isStrLowerAll(SK_"+paramonte1")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrLowerAll(SK_"+paramonte1")` must not hold.')
        assertion = assertion .and. .not. isStrLowerAll(SK_"paramonte1")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrLowerAll(SK_"paramonte1")` must not hold.')
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isStrUpperAll_1() result(assertion)
        logical(LK) :: assertion
        integer(IK) :: i
        assertion = .true._LK
        assertion = assertion .and. .not. isStrUpperAll(SK_"")
        call test%assert(assertion, MODULE_NAME//SK_': The condition `.not. isStrUpperAll(SK_"")` must hold.')
        do i = 1_IK, 127_IK
            assertion = assertion .and. (isStrUpperAll(char(i, kind = SK)) .eqv. any(ALPHA_UPPER_VEC_SK == char(i, kind = SK)))
            call test%assert(assertion, MODULE_NAME//SK_': isStrUpperAll() must properly recognize upper-case characters.')
        end do
        assertion = assertion .and. isStrUpperAll(SK_"PARAMONTE")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrUpperAll(SK_"PARAMONTE")` must hold.')
        assertion = assertion .and. .not. isStrUpperAll(SK_"PARAMONTE ")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrUpperAll(SK_"PARAMONTE ")` must not hold.')
        assertion = assertion .and. .not. isStrUpperAll(SK_"paramonte ")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrUpperAll(SK_"paramonte ")` must not hold.')
        assertion = assertion .and. .not. isStrUpperAll(SK_"paraMONTE")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrUpperAll(SK_"paraMONTE")` must not hold.')
        assertion = assertion .and. .not. isStrUpperAll(SK_"+PARAMONTE1")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrUpperAll(SK_"+PARAMONTE1")` must not hold.')
        assertion = assertion .and. .not. isStrUpperAll(SK_"PARAMONTE1")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrUpperAll(SK_"PARAMONTE1")` must not hold.')
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isStrAlphaNumAll_1() result(assertion)
        logical(LK) :: assertion
        integer(IK) :: i
        assertion = .true._LK
        assertion = assertion .and. .not. isStrAlphaNumAll(SK_"")
        call test%assert(assertion, MODULE_NAME//SK_': The condition `.not. isStrAlphaNumAll(SK_"")` must hold.')
        do i = 1_IK, 127_IK
            assertion = assertion .and. (isStrAlphaNumAll(char(i, kind = SK)) .eqv. any(ALPHANUM_VEC_SK == char(i, kind = SK)))
            call test%assert(assertion, MODULE_NAME//SK_': isStrAlphaNumAll() must properly recognize upper-case characters.')
        end do
        assertion = assertion .and. isStrAlphaNumAll(SK_"paramonte")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrAlphaNumAll(SK_"paramonte")` must hold.')
        assertion = assertion .and. isStrAlphaNumAll(SK_"paramonte1")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrAlphaNumAll(SK_"paramonte1")` must hold.')
        assertion = assertion .and. isStrAlphaNumAll(SK_"paraMONTE1")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrAlphaNumAll(SK_"paraMONTE1")` must hold.')
        assertion = assertion .and. isStrAlphaNumAll(SK_"PARAMONTE")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrAlphaNumAll(SK_"PARAMONTE")` must hold.')
        assertion = assertion .and. isStrAlphaNumAll(SK_"012345")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrAlphaNumAll(SK_"012345")` must hold.')
        assertion = assertion .and. .not. isStrAlphaNumAll(SK_"PARAMONTE ")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrAlphaNumAll(SK_"PARAMONTE ")` must not hold.')
        assertion = assertion .and. .not. isStrAlphaNumAll(SK_"paramonte+")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrAlphaNumAll(SK_"paramonte+")` must not hold.')
        assertion = assertion .and. .not. isStrAlphaNumAll(SK_" paraMONTE")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrAlphaNumAll(SK_" paraMONTE")` must not hold.')
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isStrAlpha_1() result(assertion)
        logical(LK) :: assertion
        integer(IK) :: i
        assertion = .true._LK
        assertion = assertion .and. .not. isStrAlpha(SK_"")
        call test%assert(assertion, MODULE_NAME//SK_': The condition `.not. isStrAlpha(SK_"")` must hold.')
        do i = 1_IK, 127_IK
            assertion = assertion .and. (isStrAlpha(char(i, kind = SK)) .eqv. any(ALPHA_VEC_SK == char(i, kind = SK)))
            call test%assert(assertion, MODULE_NAME//SK_': isStrAlpha() must properly recognize upper-case characters.')
        end do
        assertion = assertion .and. isStrAlpha(SK_"paramonte")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrAlpha(SK_"paramonte")` must hold.')
        assertion = assertion .and. isStrAlpha(SK_"PARAMONTE")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrAlpha(SK_"PARAMONTE")` must hold.')
        assertion = assertion .and. .not. isStrAlpha(SK_"012345")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrAlpha(SK_"012345")` must not hold.')
        assertion = assertion .and. .not. isStrAlpha(SK_"paramonte1")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrAlpha(SK_"paramonte1")` must not hold.')
        assertion = assertion .and. .not. isStrAlpha(SK_"paraMONTE1")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrAlpha(SK_"paraMONTE1")` must not hold.')
        assertion = assertion .and. .not. isStrAlpha(SK_"PARAMONTE ")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrAlpha(SK_"PARAMONTE ")` must not hold.')
        assertion = assertion .and. .not. isStrAlpha(SK_"paramonte+")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrAlpha(SK_"paramonte+")` must not hold.')
        assertion = assertion .and. .not. isStrAlpha(SK_" paraMONTE")
        call test%assert(assertion, MODULE_NAME//SK_': `isStrAlpha(SK_" paraMONTE")` must not hold.')
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getStrLower_1()  result(assertion)
    
        logical(LK)                     :: assertion
        character(:, SK), allocatable   :: string
        character(:, SK), allocatable   :: strout
        assertion = .true._LK

        string = SK_"  StringString !@#$"
        strout = getStrLower(string)
        assertion = assertion .and. strout == SK_"  stringstring !@#$"
        call test%assert(assertion, MODULE_NAME//SK_': `getStrLower(SK_"'//string//'")` must yield correct answer.')

        string = SK_""
        strout = getStrLower(string)
        assertion = assertion .and. strout == SK_""
        call test%assert(assertion, MODULE_NAME//SK_': `getStrLower(SK_"'//string//'")` must yield correct answer.')

    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_getStrUpper_1()  result(assertion)
    
        logical(LK)                     :: assertion
        character(:, SK), allocatable   :: string
        character(:, SK), allocatable   :: strout
        assertion = .true._LK

        string = SK_"  StringString !@#$"
        strout = getStrUpper(string)
        assertion = assertion .and. strout == SK_"  STRINGSTRING !@#$"
        call test%assert(assertion, MODULE_NAME//SK_': `getStrUpper(SK_"'//string//'")` must yield correct answer.')

        string = SK_""
        strout = getStrUpper(string)
        assertion = assertion .and. strout == SK_""
        call test%assert(assertion, MODULE_NAME//SK_': `getStrUpper(SK_"'//string//'")` must yield correct answer.')

    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_setStrLower_1()  result(assertion)
    
        logical(LK)                     :: assertion
        character(:, SK), allocatable   :: string
        assertion = .true._LK

        string = SK_"  StringString !@#$"
        call setStrLower(string)
        assertion = assertion .and. string == SK_"  stringstring !@#$"
        call test%assert(assertion, MODULE_NAME//SK_': `setStrLower(SK_"'//string//'")` must yield correct answer.')

        string = SK_""
        call setStrLower(string)
        assertion = assertion .and. string == SK_""
        call test%assert(assertion, MODULE_NAME//SK_': `setStrLower(SK_"'//string//'")` must yield correct answer.')

    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_setStrUpper_1()  result(assertion)
    
        logical(LK)                     :: assertion
        character(:, SK), allocatable   :: string
        assertion = .true._LK

        string = SK_"  StringString !@#$"
        call setStrUpper(string)
        assertion = assertion .and. string == SK_"  STRINGSTRING !@#$"
        call test%assert(assertion, MODULE_NAME//SK_': `setStrUpper(SK_"'//string//'")` must yield correct answer.')

        string = SK_""
        call setStrUpper(string)
        assertion = assertion .and. string == SK_""
        call test%assert(assertion, MODULE_NAME//SK_': `setStrUpper(SK_"'//string//'")` must yield correct answer.')

    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    function test_setStrAlphaRand_1()  result(assertion)
!        use pm_distUnif, only: getUnifRand
!        use pm_str, only: getCharVec
!        logical(LK)                     :: assertion
!        character(:, SK), allocatable   :: string
!        integer(IK)                     :: i
!        assertion = .true._LK
!        do i = 1_IK, 127_IK
!            allocate(character(getUnifRand(0_IK,10_IK), SK) :: string)
!            call setStrAlphaRand(string)
!            assertion = assertion .and. (isStrAlpha(string) .or. len(string) == 0)
!            call test%assert(assertion, MODULE_NAME//SK_': `setStrAlphaRand(string)` must yield characters within the alphabetic bounds: "'//string//SK_'"')
!            deallocate(string)
!        end do
!    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    function test_setStrAlphaRandLower_1()  result(assertion)
!        use pm_distUnif, only: getUnifRand
!        use pm_str, only: getCharVec
!        logical(LK)                     :: assertion
!        character(:, SK), allocatable   :: string
!        integer(IK)                     :: i
!        assertion = .true._LK
!        do i = 1_IK, 127_IK
!            allocate(character(getUnifRand(0_IK,10_IK), SK) :: string)
!            call setStrAlphaRandLower(string)
!            assertion = assertion .and. ((all(SK_"a" <= getCharVec(string)) .and. all(getCharVec(string) <= SK_"z")) .or. len(string) == 0)
!            call test%assert(assertion, MODULE_NAME//SK_': `setStrAlphaRandLower(string)` must yield characters within the alphabetic bounds: "'//string//SK_'"')
!            deallocate(string)
!        end do
!    end function
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    function test_setStrAlphaRandUpper_1()  result(assertion)
!        use pm_distUnif, only: getUnifRand
!        use pm_str, only: getCharVec
!        logical(LK)                     :: assertion
!        character(:, SK), allocatable   :: string
!        integer(IK)                     :: i
!        assertion = .true._LK
!        do i = 1_IK, 127_IK
!            allocate(character(getUnifRand(0_IK,10_IK), SK) :: string)
!            call setStrAlphaRandUpper(string)
!            assertion = assertion .and. ((all(SK_"A" <= getCharVec(string)) .and. all(getCharVec(string) <= SK_"Z")) .or. len(string) == 0)
!            call test%assert(assertion, MODULE_NAME//SK_': `setStrAlphaRandUpper(string)` must yield characters within the alphabetic bounds: "'//string//SK_'"')
!            deallocate(string)
!        end do
!    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function test_isStrComplex_1() result(assertion)
        use pm_distUnif, only: setUnifRand
        use pm_val2str, only: getStr
        implicit none
        logical(LK)                     :: assertion
        character(:, SK), allocatable   :: cmplx
        integer(IK)                     :: i
        assertion = .true._LK
        do i = 1_IK, 1000_IK
            cmplx = getComplexStr()
            assertion = assertion .and. isStrComplex(cmplx)
            call test%assert(assertion, MODULE_NAME//SK_"isStrComplex() must correctly identify """//cmplx//SK_""" as a complex number.")
            cmplx = getNonComplexStr()
            assertion = assertion .and. .not. isStrComplex(cmplx)
            call test%assert(assertion, MODULE_NAME//SK_"isStrComplex() must correctly identify """//cmplx//SK_""" as not a complex number.")
        end do
        assertion = assertion .and. .not. isStrComplex(SK_"")
        call test%assert(assertion, MODULE_NAME//SK_"isStrComplex() must correctly identify `""""` as not a complex number.")
        assertion = assertion .and. .not. isStrComplex(SK_"123")
        call test%assert(assertion, MODULE_NAME//SK_"isStrComplex() must correctly identify `""123""` as not a complex number.")
        assertion = assertion .and. .not. isStrComplex(SK_"123.0")
        call test%assert(assertion, MODULE_NAME//SK_"isStrComplex() must correctly identify `""123.0""` as not a complex number.")
        assertion = assertion .and. .not. isStrComplex(SK_"-.01e+5")
        call test%assert(assertion, MODULE_NAME//SK_"isStrComplex() must correctly identify `""-.01e+5""` as not a complex number.")
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a helper function for generating valid complex values in string format, used to verify the functionality of [isStrComplex](@ref pm_strASCII::isStrComplex).
    function getComplexStr(lp, rp, sep) result(cmplx)
        use pm_distUnif, only: getUnifRand
        character(*, SK), intent(in), optional  :: lp, rp, sep
        character(:, SK), allocatable           :: lm, rm, lp_def, rp_def, sep_def, real, imag, cmplx
        lp_def = SK_"("
        rp_def = SK_")"
        sep_def = SK_","
        if (present(lp)) lp_def = lp
        if (present(rp)) rp_def = rp
        if (present(sep)) sep_def = sep
        lm = repeat(SK_" ", getUnifRand(0_IK, 3_IK))
        rm = repeat(SK_" ", getUnifRand(0_IK, 3_IK))
        real = getRealStr()
        imag = getRealStr()
        cmplx = lm//lp_def//real//sep_def//imag//rp_def//rm
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a helper function for generating invalid complex values in string format, used to verify the functionality of [isStrComplex](@ref pm_strASCII::isStrComplex).
    function getNonComplexStr() result(cmplx)
        use pm_distUnif, only: getUnifRand
        use pm_distUnif, only: setUnifRand
        character(:, SK), allocatable           :: cmplx
        character(:, SK), allocatable           :: lm, rm, real, imag
        character(1, SK)                        :: lp, rp, sep
        if (getUnifRand()) then
            real = getRealStr()
            imag = getRealStr()
            call setUnifRand(lp)
            call setUnifRand(rp)
            call setUnifRand(sep)
            cmplx = getComplexStr(lp, rp, sep)
        else
            lm = repeat(SK_" ", getUnifRand(0_IK, 3_IK))
            rm = repeat(SK_" ", getUnifRand(0_IK, 3_IK))
            real = getNonRealStr()
            imag = getNonRealStr()
            sep = SK_","
            lp = SK_"("
            rp = SK_")"
            cmplx = lm//lp//real//sep//imag//rp//rm
        end if
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a helper function for generating valid real values in string format, used to verify the functionality of [isStrReal](@ref pm_strASCII::isStrReal).
    function getRealStr() result(real)
        use pm_distUnif, only: getUnifRand
        character(:, SK), allocatable :: real
        if (getUnifRand()) then
            real = IntSet(getUnifRand(1_IK, size(IntSet, kind = IK)))%val
        else
            real = RealSet(getUnifRand(1_IK, size(RealSet, kind = IK)))%val
        end if
        real = repeat(SK_" ", getUnifRand(0_IK, 3_IK)) // real // repeat(SK_" ", getUnifRand(0_IK, 3_IK))
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a helper function for generating invalid real values in string format, used to verify the functionality of [isStrReal](@ref pm_strASCII::isStrReal).
    function getNonRealStr() result(real)
        use pm_distUnif, only: getUnifRand
        character(:, SK), allocatable :: real
        if (getUnifRand()) then
            do
                real = NonIntSet(getUnifRand(1_IK, size(NonIntSet, kind = IK)))%val
                if (.not. isStrReal(real)) exit
            end do
        else
            real = NonRealSet(getUnifRand(1_IK, size(NonRealSet, kind = IK)))%val
        end if
        real = repeat(SK_" ", getUnifRand(0_IK, 3_IK)) // real // repeat(SK_" ", getUnifRand(0_IK, 3_IK))
    end function

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Define the global testing data.
    subroutine setData()
        IntSet  =   [ strc(SK_"0") & ! LCOV_EXCL_LINE
                    , strc(SK_"1") & ! LCOV_EXCL_LINE
                    , strc(SK_"2") & ! LCOV_EXCL_LINE
                    , strc(SK_"3") & ! LCOV_EXCL_LINE
                    , strc(SK_"4") & ! LCOV_EXCL_LINE
                    , strc(SK_"5") & ! LCOV_EXCL_LINE
                    , strc(SK_"6") & ! LCOV_EXCL_LINE
                    , strc(SK_"7") & ! LCOV_EXCL_LINE
                    , strc(SK_"8") & ! LCOV_EXCL_LINE
                    , strc(SK_"9") & ! LCOV_EXCL_LINE
                    , strc(SK_" 13 ") & ! LCOV_EXCL_LINE
                    , strc(SK_"-13 ") & ! LCOV_EXCL_LINE
                    , strc(SK_" +3 ") & ! LCOV_EXCL_LINE
                    , strc(SK_"123 ") & ! LCOV_EXCL_LINE
                    , strc(SK_"0123") & ! LCOV_EXCL_LINE
                    , strc(SK_"-123") & ! LCOV_EXCL_LINE
                    , strc(SK_" -123 ") & ! LCOV_EXCL_LINE
                    ]
        NonIntSet = [ strc(SK_"0+") & ! LCOV_EXCL_LINE
                    , strc(SK_"") & ! LCOV_EXCL_LINE
                    , strc(SK_" ") & ! LCOV_EXCL_LINE
                    , strc(SK_"+") & ! LCOV_EXCL_LINE
                    , strc(SK_"-") & ! LCOV_EXCL_LINE
                    , strc(SK_".") & ! LCOV_EXCL_LINE
                    , strc(SK_"- 123 ") & ! LCOV_EXCL_LINE
                    , strc(SK_" - 123 ") & ! LCOV_EXCL_LINE
                    , strc(SK_" +1 23 ") & ! LCOV_EXCL_LINE
                    , strc(SK_" +1 23 ") & ! LCOV_EXCL_LINE
                    , strc(SK_" +123. ") & ! LCOV_EXCL_LINE
                    , strc(SK_" +123.1 ") & ! LCOV_EXCL_LINE
                    , strc(SK_"123.1") & ! LCOV_EXCL_LINE
                    , strc(SK_"-423.1") & ! LCOV_EXCL_LINE
                    , strc(SK_"++423") & ! LCOV_EXCL_LINE
                    , strc(SK_"+-3456") & ! LCOV_EXCL_LINE
                    ]
        RealSet =   [ strc(SK_"123.") & ! LCOV_EXCL_LINE
                    , strc(SK_"-13.") & ! LCOV_EXCL_LINE
                    , strc(SK_"+13.") & ! LCOV_EXCL_LINE
                    , strc(SK_" -.1 ") & ! LCOV_EXCL_LINE
                    , strc(SK_" -.1d-12 ") & ! LCOV_EXCL_LINE
                    , strc(SK_" -.1d-123987548756 ") & ! LCOV_EXCL_LINE
                    , strc(SK_" 123.e1") & ! LCOV_EXCL_LINE
                    , strc(SK_" 123.D1 ") & ! LCOV_EXCL_LINE
                    , strc(SK_" 123.D-1") & ! LCOV_EXCL_LINE
                    , strc(SK_" 123.D+1") & ! LCOV_EXCL_LINE
                    , strc(SK_"-123.D+1 ") & ! LCOV_EXCL_LINE
                    , strc(SK_"+123.D+1 ") & ! LCOV_EXCL_LINE
                    , strc(SK_" +123.D+1") & ! LCOV_EXCL_LINE
                    , strc(SK_" +123.e+12456 ") & ! LCOV_EXCL_LINE
                    , strc(SK_" +123.E+12456 ") & ! LCOV_EXCL_LINE
                    , strc(SK_" +123.d+12456 ") & ! LCOV_EXCL_LINE
                    , strc(SK_" +123.D+12456 ") & ! LCOV_EXCL_LINE
                    , strc(SK_" +.1e-12456   ") & ! LCOV_EXCL_LINE
                    , strc(SK_" +.12456 ") & ! LCOV_EXCL_LINE
                    ]
        NonRealSet= [ strc(SK_" -.1.d-123987548756 ") & ! LCOV_EXCL_LINE
                    , strc(SK_" ") & ! LCOV_EXCL_LINE
                    , strc(SK_"") & ! LCOV_EXCL_LINE
                    , strc(SK_".123.") & ! LCOV_EXCL_LINE
                    , strc(SK_"123 .") & ! LCOV_EXCL_LINE
                    , strc(SK_" 123.+1") & ! LCOV_EXCL_LINE
                    , strc(SK_" 123.D 1") & ! LCOV_EXCL_LINE
                    , strc(SK_" 123.D +1") & ! LCOV_EXCL_LINE
                    , strc(SK_" 123.D + 1") & ! LCOV_EXCL_LINE
                    , strc(SK_" 123.D+ 1") & ! LCOV_EXCL_LINE
                    , strc(SK_" 123. D+ 1") & ! LCOV_EXCL_LINE
                    , strc(SK_" + 123.D+1") & ! LCOV_EXCL_LINE
                    , strc(SK_" +123.D+-12456 ") & ! LCOV_EXCL_LINE
                    , strc(SK_" +123.D+.12456 ") & ! LCOV_EXCL_LINE
                    , strc(SK_" +.d+12456") & ! LCOV_EXCL_LINE
                    , strc(SK_" -.1.d-123987548756 ") & ! LCOV_EXCL_LINE
                    , strc(SK_" +.e-12456   ") & ! LCOV_EXCL_LINE
                    , strc(SK_"+ . ") & ! LCOV_EXCL_LINE
                    , strc(SK_"+. ") & ! LCOV_EXCL_LINE
                    , strc(SK_"+ .") & ! LCOV_EXCL_LINE
                    ]
    end subroutine

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setTest()

        implicit none

        call setData()
        test = test_type(MODULE_NAME)
        call test%run(test_isStrReal_1, SK_"test_isStrReal_1")
        call test%run(test_isStrAlpha_1, SK_"test_isStrAlpha_1")
        call test%run(test_isCharDigit_1, SK_"test_isCharDigit_1")
        call test%run(test_isStrNumber_1, SK_"test_isStrNumber_1")
        call test%run(test_isStrDigitAll_1, SK_"test_isStrDigitAll_1")
        call test%run(test_isStrInteger_1, SK_"test_isStrInteger_1")
        call test%run(test_isStrComplex_1, SK_"test_isStrComplex_1")
        call test%run(test_isStrAlphaNumAll_1, SK_"test_isStrAlphaNumAll_1")
        call test%run(test_isCharLower_1, SK_"test_isCharLower_1")
        call test%run(test_isCharUpper_1, SK_"test_isCharUpper_1")
        call test%run(test_getLocSpace_1, SK_"test_getLocSpace_1")
        call test%run(test_isStrLowerAll_1, SK_"test_isStrLowerAll_1")
        call test%run(test_isStrUpperAll_1, SK_"test_isStrUpperAll_1")
        call test%run(test_getStrLower_1, SK_"test_getStrLower_1")
        call test%run(test_getStrUpper_1, SK_"test_getStrUpper_1")
        call test%run(test_setStrLower_1, SK_"test_setStrLower_1")
        call test%run(test_setStrUpper_1, SK_"test_setStrUpper_1")
!        call test%run(test_setStrAlphaRand_1, SK_"test_setStrAlphaRand_1")
!        call test%run(test_setStrAlphaRandLower_1, SK_"test_setStrAlphaRandLower_1")
!        call test%run(test_setStrAlphaRandUpper_1, SK_"test_setStrAlphaRandUpper_1")

        call test%summarize()

    end subroutine setTest

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module test_pm_strASCII ! LCOV_EXCL_LINE