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
!>  This module contains the uncommon and hardly representable ASCII characters as well as
!>  procedures for operating on strings that exclusively contain the 128 ASCII characters.<br>
!>
!>  \details
!>  ASCII, stands for American Standard Code for Information Interchange.<br>
!>  It is a 7-bit character code where every single bit represents a unique character.<br>
!>  The first 32 characters in the ASCII-table are unprintable control codes and are used to control peripherals such as printers.<br>
!>  ASCII printable characters (character code 32-127) Codes 32-127 are common for all the different variations of the ASCII table,
!>  they are called printable characters, represent letters, digits, punctuation marks, and a few miscellaneous symbols.<br>
!>  You will find almost every character on your keyboard. Character 127 represents the command DEL.<br>
!>
!>  \see
!>  [pm_distUnif](@ref pm_distUnif) for uniformly-distributed random character generation.<br>
!>  [https://www.ascii-code.com/](https://www.ascii-code.com/) for information about ASCII characters.<br>
!>
!>  \benchmarks
!>
!>  \benchmark{getStrLower_vs_setStrLower, The runtime performance of [getStrLower](@ref pm_strASCII::getStrLower) vs. [setStrLower](@ref pm_strASCII::setStrLower)}
!>  \include{lineno} benchmark/pm_strASCII/getStrLower_vs_setStrLower/main.F90
!>  \compilefb{getStrLower_vs_setStrLower}
!>  \postprocb{getStrLower_vs_setStrLower}
!>  \include{lineno} benchmark/pm_strASCII/getStrLower_vs_setStrLower/main.py
!>  \visb{getStrLower_vs_setStrLower}
!>  \image html benchmark/pm_strASCII/getStrLower_vs_setStrLower/benchmark.getStrLower_vs_setStrLower.runtime.png width=1000
!>  \image html benchmark/pm_strASCII/getStrLower_vs_setStrLower/benchmark.getStrLower_vs_setStrLower.runtime.ratio.png width=1000
!>  \moralb{getStrLower_vs_setStrLower}
!>      -#  The procedure [getStrLower](@ref pm_strASCII::getStrLower) is a function while
!>          the procedure [setStrLower](@ref pm_strASCII::setStrLower) is a subroutines.<br>
!>          From the benchmark results, it appears that the functional interface performs less efficiently than
!>          the subroutine interface when the input `array` size is small (\f$<300\f$).<br>
!>          Nevertheless, the difference appears to diminish or even change its course for large string sizes.<br>
!>      -#  The results of this benchmark are equally applicable to the [getStrUpper](@ref getStrUpper) vs. [setStrUpper](@ref setStrUpper).<br>
!>
!>  \test
!>  [test_pm_strASCII](@ref test_pm_strASCII)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_strASCII

    use pm_kind, only: SK, IK, LK

    implicit none

    !>  \cond excluded
    integer, private :: i
    !>  \endcond excluded

    character(*, SK), parameter :: MODULE_NAME = SK_"@pm_strASCII"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    integer(IK)     , parameter :: MAX_ASCII_CODE = 127_IK  !<  \public The scalar `integer` constant of default kind \IK representing the maximum standard ASCII code.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ASCII Control characters.

    character(1, SK), parameter :: NUL  = achar(0, SK)      !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Null character.<br>
                                                            !<          \warning    This is a dangerous character that has a special meaning to C-processors and some Fortran compilers like ifort.<br>
                                                            !<                      As such, it should never be used in any string unless the consequences are fully understood.<br>
    character(1, SK), parameter :: SOH  = achar(1, SK)      !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Start of Heading.<br>
    character(1, SK), parameter :: STX  = achar(2, SK)      !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Start of Text.<br>
    character(1, SK), parameter :: ETX  = achar(3, SK)      !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: End of Text.<br>
    character(1, SK), parameter :: EOT  = achar(4, SK)      !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: End of Transmission.<br>
    character(1, SK), parameter :: ENQ  = achar(5, SK)      !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Enquiry.<br>
    character(1, SK), parameter :: ACK  = achar(6, SK)      !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Acknowledgment.<br>
    character(1, SK), parameter :: BEL  = achar(7, SK)      !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Bell (ring character).<br>
    character(1, SK), parameter :: BS   = achar(8, SK)      !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Back Space.<br>
    character(1, SK), parameter :: HT   = achar(9, SK)      !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Horizontal Tab.<br>
    character(1, SK), parameter :: LF   = achar(10, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Line Feed (new line character).<br>
    character(1, SK), parameter :: VT   = achar(11, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Vertical Tab.<br>
    character(1, SK), parameter :: FF   = achar(12, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Form Feed.<br>
    character(1, SK), parameter :: CR   = achar(13, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Carriage Return.<br>
    character(1, SK), parameter :: SO   = achar(14, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Shift Out / X-On.<br>
    character(1, SK), parameter :: SI   = achar(15, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Shift In / X-Off.<br>
    character(1, SK), parameter :: DLE  = achar(16, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Data Line Escape.<br>
    character(1, SK), parameter :: DC1  = achar(17, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Device Control 1 (oft. XON).<br>
    character(1, SK), parameter :: DC2  = achar(18, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Device Control 2.<br>
    character(1, SK), parameter :: DC3  = achar(19, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Device Control 3 (oft. XOFF).<br>
    character(1, SK), parameter :: DC4  = achar(20, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Device Control 4.<br>
    character(1, SK), parameter :: NAK  = achar(21, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Negative Acknowledgement.<br>
    character(1, SK), parameter :: SYN  = achar(22, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Synchronous Idle.<br>
    character(1, SK), parameter :: ETB  = achar(23, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: End of Transmit Block.<br>
    character(1, SK), parameter :: CAN  = achar(24, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Cancel.<br>
    character(1, SK), parameter :: EM   = achar(25, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: End of Medium.<br>
    character(1, SK), parameter :: SUB  = achar(26, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Substitute.<br>
    character(1, SK), parameter :: ESC  = achar(27, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Escape.<br>
    character(1, SK), parameter :: FS   = achar(28, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: File Separator.<br>
    character(1, SK), parameter :: GS   = achar(29, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Group Separator.<br>
    character(1, SK), parameter :: RS   = achar(30, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Record Separator.<br>
    character(1, SK), parameter :: US   = achar(31, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Unit Separator.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Windows line terminator.
    character(2, SK), parameter :: CRLF = CR//LF            !<  \public The scalar `character` constant of default kind \SK of length `2` representing ASCII
                                                            !!          characters Carriage Return and Line Feed that together form the line terminator on Windows systems.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ASCII printable characters.

    character(1, SK), parameter :: SPC  = achar(32, SK)     !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Space.<br>
    character(1, SK), parameter :: DEL  = achar(127, SK)    !<  \public The scalar `character` constant of default kind \SK of length `1` representing ASCII character: Delete.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `character` of default kind \SK containing the ASCII uppercase English letters.<br>
    character(*, SK), parameter :: ALPHA_UPPER_STR_SK = SK_"ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    !>  \brief
    !>  The constant array of type `character` of default kind \SK containing the uppercase English letters.<br>
    character(1, SK), parameter :: ALPHA_UPPER_VEC_SK(*) = [(ALPHA_UPPER_STR_SK(i:i), i = 1, len(ALPHA_UPPER_STR_SK))]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: ALPHA_UPPER_VEC_SK
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `character` of default kind \SK containing the ASCII lowercase English letters.<br>
    character(*, SK), parameter :: ALPHA_LOWER_STR_SK = SK_"abcdefghijklmnopqrstuvwxyz"

    !>  \brief
    !>  The constant array of type `character` of default kind \SK containing the lowercase English letters.<br>
    character(1, SK), parameter :: ALPHA_LOWER_VEC_SK(*) = [(ALPHA_LOWER_STR_SK(i:i), i = 1, len(ALPHA_LOWER_STR_SK))]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: ALPHA_LOWER_VEC_SK
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `character` of default kind \SK containing the ASCII uppercase and lowercase English letters.<br>
    character(*, SK), parameter :: ALPHA_STR_SK = ALPHA_UPPER_STR_SK//ALPHA_LOWER_STR_SK

    !>  \brief
    !>  The constant array of type `character` of default kind \SK containing the ASCII uppercase and lowercase English letters.<br>
    character(1, SK), parameter :: ALPHA_VEC_SK(*) = [ALPHA_UPPER_VEC_SK, ALPHA_LOWER_VEC_SK]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: ALPHA_VEC_SK
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `integer` of default kind \IK representing the distance between `A` and `a` in the **processor** collating sequence.<br>
    integer(IK)     , parameter :: UPPER_MINUS_LOWER_IK = ichar(SK_"A", kind = IK) - ichar(SK_"a", kind = IK)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `character` of default kind \SK containing the English numbers.<br>
    character(*, SK), parameter :: DIGIT_STR_SK = SK_"0123456789"

    !>  \brief
    !>  The constant array of type `character` of default kind \SK containing the English numbers.<br>
    character(1, SK), parameter :: DIGIT_VEC_SK(*) = [(DIGIT_STR_SK(i:i), i = 1, len(DIGIT_STR_SK))]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: DIGIT_VEC_SK
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `character` of default kind \SK containing the English numbers and alphabets.<br>
    character(*, SK), parameter :: ALPHANUM_STR_SK = DIGIT_STR_SK//ALPHA_STR_SK

    !>  \brief
    !>  The constant array of type `character` of default kind \SK containing the English numbers and alphabets.<br>
    character(1, SK), parameter :: ALPHANUM_VEC_SK(*) = [DIGIT_VEC_SK, ALPHA_VEC_SK]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: ALPHANUM_VEC_SK
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `character` of default kind \SK containing characters that might be used for integer number representation.<br>
    character(*, SK), parameter :: INTEGER_STR_SK = DIGIT_STR_SK//SK_"-+"

    !>  \brief
    !>  The constant array of type `character` of default kind \SK containing characters that might be used for integer number representation.<br>
    character(1, SK), parameter :: INTEGER_VEC_SK(*) = [(INTEGER_STR_SK(i:i), i = 1, len(INTEGER_STR_SK))]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: INTEGER_VEC_SK
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `character` of default kind \SK containing characters that might be used for real number representation.<br>
    character(*, SK), parameter :: REAL_STR_SK = DIGIT_STR_SK//SK_".edED"

    !>  \brief
    !>  The constant array of type `character` of default kind \SK containing characters that might be used for real number representation.<br>
    character(1, SK), parameter :: REAL_VEC_SK(*) = [(REAL_STR_SK(i:i), i = 1, len(REAL_STR_SK))]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: REAL_VEC_SK
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  The constant scalar of type `character` of default kind \SK containing characters that might be used for complex number representation.<br>
    character(*, SK), parameter :: COMPLEX_STR_SK = REAL_STR_SK//SK_"(),"

    !>  \brief
    !>  The constant array of type `character` of default kind \SK containing characters that might be used for complex number representation.<br>
    character(1, SK), parameter :: COMPLEX_VEC_SK(*) = [(COMPLEX_STR_SK(i:i), i = 1, len(COMPLEX_STR_SK))]
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: COMPLEX_VEC_SK
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the location of the first occurrence of the whitespace character in the input string, from left to right.<br>
    !>  Otherwise, return `0` if the input string does not contain a whitespace character.<br>
    !>
    !>  \param[in]  str :   The input scalar or array of arbitrary rank of type `character` of default kind \SK.<br>
    !>
    !>  \return
    !>  `locSpace`      :   The output scalar or array of the same shape as the input `str` of type `integer` of default kind \IK
    !>                      whose value is the location of the first occurrence of the whitespace character in the input string
    !>                      from the left to the right.<br>
    !>
    !>  \interface{getLocSpace}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: getLocSpace
    !>      use pm_kind, only: IK
    !>      integer(IK) :: locSpace
    !>
    !>      locSpace = getLocSpace(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  This procedure exists primarily to facilitate coding of other procedures in this module.<br>
    !>  Other algorithms in this library listed below can provide more general solutions to the problem solved by this procedure.<br>
    !>
    !>  \see
    !>  [getLocNonSpace](@ref pm_strASCII::getLocNonSpace)<br>
    !>  [setLoc](@ref pm_arrayFind::setLoc)<br>
    !>  [getLoc](@ref pm_arrayFind::getLoc)<br>
    !>  [getBin](@ref pm_arraySearch::getBin)<br>
    !>
    !>  \example{getLocSpace}
    !>  \include{lineno} example/pm_strASCII/getLocSpace/main.F90
    !>  \compilef{getLocSpace}
    !>  \output{getLocSpace}
    !>  \include{lineno} example/pm_strASCII/getLocSpace/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getLocSpace

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function getLocSpace_SK5(str) result(locSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocSpace_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        integer(IK)                                             :: locSpace
    end function
#endif

#if SK4_ENABLED
    pure elemental module function getLocSpace_SK4(str) result(locSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocSpace_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        integer(IK)                                             :: locSpace
    end function
#endif

#if SK3_ENABLED
    pure elemental module function getLocSpace_SK3(str) result(locSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocSpace_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        integer(IK)                                             :: locSpace
    end function
#endif

#if SK2_ENABLED
    pure elemental module function getLocSpace_SK2(str) result(locSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocSpace_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        integer(IK)                                             :: locSpace
    end function
#endif

#if SK1_ENABLED
    pure elemental module function getLocSpace_SK1(str) result(locSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocSpace_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        integer(IK)                                             :: locSpace
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the location of the first occurrence of a non-whitespace character in the input string, from left to right.<br>
    !>  Otherwise, return `0` if the input string is all whitespace characters.<br>
    !>
    !>  \param[in]  str :   The input scalar or array of arbitrary rank of type `character` of default kind \SK.<br>
    !>
    !>  \return
    !>  `locNonSpace`   :   The output scalar or array of the same shape as the input `str` of type `integer` of default kind \IK
    !>                      whose value is the location of the first occurrence of the whitespace character in the input string
    !>                      from the left to the right.<br>
    !>
    !>  \interface{getLocNonSpace}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: getLocNonSpace
    !>      use pm_kind, only: IK
    !>      integer(IK) :: locNonSpace
    !>
    !>      locNonSpace = getLocNonSpace(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \remark
    !>  This procedure exists primarily to facilitate coding of other procedures in this module.<br>
    !>  Other algorithms in this library listed below can provide more general solutions to the problem solved by this procedure.<br>
    !>
    !>  \see
    !>  [getLocNonSpace](@ref pm_strASCII::getLocNonSpace)<br>
    !>  [setLoc](@ref pm_arrayFind::setLoc)<br>
    !>  [getLoc](@ref pm_arrayFind::getLoc)<br>
    !>  [getBin](@ref pm_arraySearch::getBin)<br>
    !>
    !>  \example{getLocNonSpace}
    !>  \include{lineno} example/pm_strASCII/getLocNonSpace/main.F90
    !>  \compilef{getLocNonSpace}
    !>  \output{getLocNonSpace}
    !>  \include{lineno} example/pm_strASCII/getLocNonSpace/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getLocNonSpace

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function getLocNonSpace_SK5(str) result(locNonSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocNonSpace_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        integer(IK)                                             :: locNonSpace
    end function
#endif

#if SK4_ENABLED
    pure elemental module function getLocNonSpace_SK4(str) result(locNonSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocNonSpace_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        integer(IK)                                             :: locNonSpace
    end function
#endif

#if SK3_ENABLED
    pure elemental module function getLocNonSpace_SK3(str) result(locNonSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocNonSpace_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        integer(IK)                                             :: locNonSpace
    end function
#endif

#if SK2_ENABLED
    pure elemental module function getLocNonSpace_SK2(str) result(locNonSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocNonSpace_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        integer(IK)                                             :: locNonSpace
    end function
#endif

#if SK1_ENABLED
    pure elemental module function getLocNonSpace_SK1(str) result(locNonSpace)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLocNonSpace_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        integer(IK)                                             :: locNonSpace
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input character of length `1` is a digit as in [DIGIT_VEC_SK](@ref pm_strASCII::DIGIT_VEC_SK).<br>
    !>
    !>  \param[in]  chr     :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL of length `1`.<br>
    !>
    !>  \return
    !>  `charIsDigit`       :   The output scalar or array of the same shape as the input `chr` of type `logical` of default kind \LK
    !>                          whose value is `.true.` corresponding to each element of the input `chr` if the element is a digit,
    !>                          otherwise it is `.false.`.<br>
    !>
    !>  \interface{isCharDigit}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: isCharDigit
    !>      use pm_kind, only: LK
    !>      logical(LK) :: charIsDigit
    !>
    !>      charIsDigit = isCharDigit(chr)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>  [isStrAlphaNumAll](@ref pm_strASCII::isStrAlphaNumAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isCharDigit}
    !>  \include{lineno} example/pm_strASCII/isCharDigit/main.F90
    !>  \compilef{isCharDigit}
    !>  \output{isCharDigit}
    !>  \include{lineno} example/pm_strASCII/isCharDigit/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isCharDigit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isCharDigit_SK5(chr) result(charIsDigit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharDigit_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsDigit
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isCharDigit_SK4(chr) result(charIsDigit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharDigit_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsDigit
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isCharDigit_SK3(chr) result(charIsDigit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharDigit_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsDigit
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isCharDigit_SK2(chr) result(charIsDigit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharDigit_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsDigit
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isCharDigit_SK1(chr) result(charIsDigit)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharDigit_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsDigit
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if **all characters** of the input string are numeric (digits as in [DIGIT_VEC_SK](@ref pm_strASCII::DIGIT_VEC_SK)).<br>
    !>
    !>  \param[in]  str     :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `strIsDigitAll`     :   The output scalar or array of the same shape as the input `str` of type `logical` of default kind \LK
    !>                          whose value is `.true.` if **all** characters of the corresponding element of `str` are numeric,
    !>                          otherwise it is `.false.`.<br>
    !>
    !>  \interface{isStrDigitAll}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: isStrDigitAll
    !>      use pm_kind, only: LK
    !>      logical(LK) :: strIsDigit
    !>
    !>      strIsDigit = isStrDigitAll(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The functionality of this interface can be also replicated by [isStrDigit](@ref pm_strASCII::isStrDigit),
    !>  \code{.F90}
    !>
    !>      character(10) :: str
    !>
    !>      isStrDigitAll(str) == all(isStrDigit(str))
    !>
    !>  \endcode
    !>  although this generic interface is potentially faster than [isStrDigit](@ref pm_strASCII::isStrDigit).<br>
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>  [isStrAlphaNumAll](@ref pm_strASCII::isStrAlphaNumAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrDigitAll}
    !>  \include{lineno} example/pm_strASCII/isStrDigitAll/main.F90
    !>  \compilef{isStrDigitAll}
    !>  \output{isStrDigitAll}
    !>  \include{lineno} example/pm_strASCII/isStrDigitAll/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrDigitAll

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isStrDigitAll_SK5(str) result(strIsDigitAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrDigitAll_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsDigitAll
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isStrDigitAll_SK4(str) result(strIsDigitAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrDigitAll_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsDigitAll
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isStrDigitAll_SK3(str) result(strIsDigitAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrDigitAll_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsDigitAll
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isStrDigitAll_SK2(str) result(strIsDigitAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrDigitAll_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsDigitAll
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isStrDigitAll_SK1(str) result(strIsDigitAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrDigitAll_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsDigitAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if **any character** of the input string is numeric (digits as in [DIGIT_VEC_SK](@ref pm_strASCII::DIGIT_VEC_SK)).<br>
    !>
    !>  \param[in]  str     :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `strIsDigitAny`     :   The output scalar or array of the same shape as the input `str` of type `logical` of default kind \LK
    !>                          whose value is `.true.` if **any character** of the corresponding element of `str` is numeric,
    !>                          otherwise it is `.false.`.<br>
    !>
    !>  \interface{isStrDigitAny}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: isStrDigitAny
    !>      use pm_kind, only: LK
    !>      logical(LK) :: strIsDigitAny
    !>
    !>      strIsDigitAny = isStrDigitAny(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The functionality of this interface can be also replicated by [isStrDigit](@ref pm_strASCII::isStrDigit),
    !>  \code{.F90}
    !>
    !>      character(10) :: str
    !>
    !>      isStrDigitAny(str) == any(isStrDigit(str))
    !>
    !>  \endcode
    !>  although this generic interface is potentially faster than [isStrDigit](@ref pm_strASCII::isStrDigit).<br>
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAny](@ref pm_strASCII::isStrDigitAny)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>  [isStrAlphaNumAll](@ref pm_strASCII::isStrAlphaNumAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrDigitAny}
    !>  \include{lineno} example/pm_strASCII/isStrDigitAny/main.F90
    !>  \compilef{isStrDigitAny}
    !>  \output{isStrDigitAny}
    !>  \include{lineno} example/pm_strASCII/isStrDigitAny/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrDigitAny

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isStrDigitAny_SK5(str) result(strIsDigitAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrDigitAny_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsDigitAny
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isStrDigitAny_SK4(str) result(strIsDigitAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrDigitAny_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsDigitAny
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isStrDigitAny_SK3(str) result(strIsDigitAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrDigitAny_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsDigitAny
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isStrDigitAny_SK2(str) result(strIsDigitAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrDigitAny_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsDigitAny
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isStrDigitAny_SK1(str) result(strIsDigitAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrDigitAny_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsDigitAny
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a vector of logical values of the same size as the length of the input scalar string `str`, whose elements
    !>  are `.true.` <b>if and only if</b> the corresponding character in the input `str` is numeric (digits as in [DIGIT_VEC_SK](@ref pm_strASCII::DIGIT_VEC_SK)).<br>
    !>
    !>  \param[in]  str :   The input scalar of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `StrIsNumeric`  :   The output array of rank `1` of the same size as the length of the input scalar `str`, of type `logical` of default kind \LK,
    !>                      whose elements are `.true.` <b>if and only if</b> the corresponding character in the input `str` is numeric,
    !>                      otherwise it is `.false.`.<br>
    !>
    !>  \interface{isStrDigit}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: isStrDigit
    !>      use pm_kind, only: LK
    !>      logical(LK) :: StrIsNumeric(len(str,IK))
    !>
    !>      StrIsNumeric = isStrDigit(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigit](@ref pm_strASCII::isStrDigit)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>  [isStrAlphaNumAll](@ref pm_strASCII::isStrAlphaNumAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrDigit}
    !>  \include{lineno} example/pm_strASCII/isStrDigit/main.F90
    !>  \compilef{isStrDigit}
    !>  \output{isStrDigit}
    !>  \include{lineno} example/pm_strASCII/isStrDigit/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrDigit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function isStrDigit_SK5(str) result(StrIsNumeric)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrDigit_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsNumeric(len(str,IK))
    end function
#endif

#if SK4_ENABLED
    pure module function isStrDigit_SK4(str) result(StrIsNumeric)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrDigit_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsNumeric(len(str,IK))
    end function
#endif

#if SK3_ENABLED
    pure module function isStrDigit_SK3(str) result(StrIsNumeric)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrDigit_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsNumeric(len(str,IK))
    end function
#endif

#if SK2_ENABLED
    pure module function isStrDigit_SK2(str) result(StrIsNumeric)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrDigit_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsNumeric(len(str,IK))
    end function
#endif

#if SK1_ENABLED
    pure module function isStrDigit_SK1(str) result(StrIsNumeric)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrDigit_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsNumeric(len(str,IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if all characters of the input string collectively represent a real number.<br>
    !>
    !>  \details
    !>  The characters used in the representation of a real number are given by [REAL_VEC_SK](@ref pm_strASCII::REAL_VEC_SK)).<br>
    !>  The `real` string must contain only these characters and must collectively represent a meaningful real number.<br>
    !>
    !>  \param[in]  str :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `strIsReal`     :   The output scalar or array of the same shape as the input `str` of type `logical` of default kind \LK
    !>                      whose value is `.true.` corresponding to each element of the input `str` if the element is a real
    !>                      number, otherwise it is `.false.`.<br>
    !>
    !>  \interface{isStrReal}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: isStrReal
    !>      use pm_kind, only: LK
    !>      logical(LK) :: strIsReal
    !>
    !>      strIsReal = isStrReal(str)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  This procedure is sensitive to the presence of white-space characters in the string, in which case the output is `.false.`.<br>
    !>  The must not exist any white-space within the real number or at the beginning or the end of the string.<br>
    !>  If there is the possibility for any trailing or preceding white-space characters in the string,
    !>  pass `trim(adjustl(str))` as the input argument instead.<br>
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>  [isStrAlphaNumAll](@ref pm_strASCII::isStrAlphaNumAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrReal}
    !>  \include{lineno} example/pm_strASCII/isStrReal/main.F90
    !>  \compilef{isStrReal}
    !>  \output{isStrReal}
    !>  \include{lineno} example/pm_strASCII/isStrReal/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrReal

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isStrReal_SK5(str) result(strIsReal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrReal_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsReal
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isStrReal_SK4(str) result(strIsReal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrReal_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsReal
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isStrReal_SK3(str) result(strIsReal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrReal_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsReal
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isStrReal_SK2(str) result(strIsReal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrReal_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsReal
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isStrReal_SK1(str) result(strIsReal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrReal_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsReal
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if all characters of the input string collectively represent a complex number.<br>
    !>
    !>  \details
    !>  The characters used in the representation of a complex number are given by [COMPLEX_VEC_SK](@ref pm_strASCII::COMPLEX_VEC_SK)).<br>
    !>  The `complex` string must contain only these characters and must collectively represent a meaningful complex number.<br>
    !>
    !>  \param[in]  str :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `strIsComplex`  :   The output scalar or array of the same shape as the input `str` of type `logical` of default kind \LK
    !>                      whose value is `.true.` corresponding to each element of the input `str` if the element is a complex
    !>                      number, otherwise it is `.false.`.<br>
    !>
    !>  \interface{isStrComplex}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: isStrComplex
    !>      use pm_kind, only: LK
    !>      logical(LK) :: strIsComplex
    !>
    !>      strIsComplex = isStrComplex(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>  [isStrAlphaNumAll](@ref pm_strASCII::isStrAlphaNumAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrComplex}
    !>  \include{lineno} example/pm_strASCII/isStrComplex/main.F90
    !>  \compilef{isStrComplex}
    !>  \output{isStrComplex}
    !>  \include{lineno} example/pm_strASCII/isStrComplex/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrComplex

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isStrComplex_SK5(str) result(strIsComplex)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrComplex_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsComplex
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isStrComplex_SK4(str) result(strIsComplex)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrComplex_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsComplex
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isStrComplex_SK3(str) result(strIsComplex)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrComplex_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsComplex
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isStrComplex_SK2(str) result(strIsComplex)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrComplex_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsComplex
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isStrComplex_SK1(str) result(strIsComplex)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrComplex_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsComplex
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if all characters of the input string collectively represent an integer.<br>
    !>
    !>  \details
    !>  The characters used in the representation of an integer number are given by [INTEGER_VEC_SK](@ref pm_strASCII::INTEGER_VEC_SK)).<br>
    !>  The `integer` string must contain only these characters and must collectively represent a meaningful integer number.<br>
    !>
    !>  \param[in]  str :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `strIsInteger`  :   The output scalar or array of the same shape as the input `str` of type `logical` of default kind \LK
    !>                      whose value is `.true.` corresponding to each element of the input `str` if the element is an integer
    !>                      number, otherwise it is `.false.`.<br>
    !>
    !>  \interface{isStrInteger}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_strASCII, only: isStrInteger
    !>      logical(LK) :: strIsInteger
    !>
    !>      strIsInteger = isStrInteger(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>  [isStrAlphaNumAll](@ref pm_strASCII::isStrAlphaNumAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrInteger}
    !>  \include{lineno} example/pm_strASCII/isStrInteger/main.F90
    !>  \compilef{isStrInteger}
    !>  \output{isStrInteger}
    !>  \include{lineno} example/pm_strASCII/isStrInteger/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrInteger

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isStrInteger_SK5(str) result(strIsInteger)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrInteger_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsInteger
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isStrInteger_SK4(str) result(strIsInteger)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrInteger_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsInteger
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isStrInteger_SK3(str) result(strIsInteger)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrInteger_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsInteger
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isStrInteger_SK2(str) result(strIsInteger)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrInteger_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsInteger
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isStrInteger_SK1(str) result(strIsInteger)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrInteger_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsInteger
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input string collectively represents an integer, a real, or a complex number.<br>
    !>
    !>  \param[in]  str :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `strIsNumber`   :   The output scalar or array of the same shape as the input `str` of type `logical` of default kind \LK
    !>                      whose value is `.true.` corresponding to each element of the input `str` if the element is an integer,
    !>                      a real, or a complex number, otherwise it is `.false.`.<br>
    !>
    !>  \interface{isStrNumber}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: isStrNumber
    !>      use pm_kind, only: LK
    !>      logical(LK) :: strIsNumber
    !>
    !>      strIsNumber = isStrNumber(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>  [isStrAlphaNumAll](@ref pm_strASCII::isStrAlphaNumAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrNumber}
    !>  \include{lineno} example/pm_strASCII/isStrNumber/main.F90
    !>  \compilef{isStrNumber}
    !>  \output{isStrNumber}
    !>  \include{lineno} example/pm_strASCII/isStrNumber/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrNumber

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isStrNumber_SK5(str) result(strIsNumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrNumber_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsNumber
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isStrNumber_SK4(str) result(strIsNumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrNumber_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsNumber
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isStrNumber_SK3(str) result(strIsNumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrNumber_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsNumber
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isStrNumber_SK2(str) result(strIsNumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrNumber_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsNumber
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isStrNumber_SK1(str) result(strIsNumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrNumber_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsNumber
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input single character contains an uppercase English alphabet [ALPHA_UPPER_VEC_SK](@ref ALPHA_UPPER_VEC_SK).<br>
    !>
    !>  \param[in]  chr     :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL of length `1`.<br>
    !>
    !>  \return
    !>  `charIsUpper`       :   The output scalar or array of the same shape as the input `chr` of type `logical` of default kind \LK
    !>                          whose value is `.true.` corresponding to each element of the input `chr` if the element is an uppercase
    !>                          English letter, otherwise it is `.false.`.<br>
    !>
    !>  \interface{isCharUpper}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_strASCII, only: isCharUpper
    !>      logical(LK) :: charIsUpper
    !>
    !>      charIsUpper = isCharUpper(chr)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>  [isStrAlphaNumAll](@ref pm_strASCII::isStrAlphaNumAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isCharUpper}
    !>  \include{lineno} example/pm_strASCII/isCharUpper/main.F90
    !>  \compilef{isCharUpper}
    !>  \output{isCharUpper}
    !>  \include{lineno} example/pm_strASCII/isCharUpper/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isCharUpper

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isCharUpper_SK5(chr) result(charIsUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharUpper_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsUpper
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isCharUpper_SK4(chr) result(charIsUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharUpper_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsUpper
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isCharUpper_SK3(chr) result(charIsUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharUpper_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsUpper
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isCharUpper_SK2(chr) result(charIsUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharUpper_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsUpper
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isCharUpper_SK1(chr) result(charIsUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharUpper_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsUpper
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input single character contains an lowercase English alphabet [ALPHA_LOWER_VEC_SK](@ref ALPHA_LOWER_VEC_SK).<br>
    !>
    !>  \param[in]  chr :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL of length `1`.<br>
    !>
    !>  \return
    !>  `charIsLower`   :   The output scalar or array of the same shape as the input `chr` of type `logical` of default kind \LK
    !>                      whose value is `.true.` corresponding to each element of the input `chr` if the element is an lowercase
    !>                      English letter, otherwise it is `.false.`.<br>
    !>
    !>  \interface{isCharLower}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_strASCII, only: isCharLower
    !>      logical(LK) :: charIsLower
    !>
    !>      charIsLower = isCharLower(chr)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>  [isStrAlphaNumAll](@ref pm_strASCII::isStrAlphaNumAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isCharLower}
    !>  \include{lineno} example/pm_strASCII/isCharLower/main.F90
    !>  \compilef{isCharLower}
    !>  \output{isCharLower}
    !>  \include{lineno} example/pm_strASCII/isCharLower/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isCharLower

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isCharLower_SK5(chr) result(charIsLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharLower_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsLower
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isCharLower_SK4(chr) result(charIsLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharLower_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsLower
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isCharLower_SK3(chr) result(charIsLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharLower_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsLower
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isCharLower_SK2(chr) result(charIsLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharLower_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsLower
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isCharLower_SK1(chr) result(charIsLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharLower_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsLower
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input string contains **only** uppercase English alphabets [ALPHA_UPPER_VEC_SK](@ref ALPHA_UPPER_VEC_SK).<br>
    !>
    !>  \param[in]  str :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `strIsUpperAll` :   The output scalar or array of the same shape as the input `str` of type `logical` of default kind \LK
    !>                      whose value is `.true.` corresponding to each element of the input `str` if the element is an uppercase
    !>                      English letter, otherwise it is `.false.`.<br>
    !>
    !>  \interface{isStrUpperAll}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_strASCII, only: isStrUpperAll
    !>      logical(LK) :: strIsUpperAll
    !>
    !>      strIsUpperAll = isStrUpperAll(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The functionality of this interface can be also replicated by [isStrUpper](@ref pm_strASCII::isStrUpper),
    !>  \code{.F90}
    !>
    !>      character(10) :: str
    !>
    !>      isStrUpperAll(str) == all(isStrUpper(str))
    !>
    !>  \endcode
    !>  although this generic interface is potentially faster than [isStrUpper](@ref pm_strASCII::isStrUpper).<br>
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>  [isStrAlphaNumAll](@ref pm_strASCII::isStrAlphaNumAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrUpperAll}
    !>  \include{lineno} example/pm_strASCII/isStrUpperAll/main.F90
    !>  \compilef{isStrUpperAll}
    !>  \output{isStrUpperAll}
    !>  \include{lineno} example/pm_strASCII/isStrUpperAll/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrUpperAll

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isStrUpperAll_SK5(str) result(strIsUpperAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrUpperAll_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsUpperAll
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isStrUpperAll_SK4(str) result(strIsUpperAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrUpperAll_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsUpperAll
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isStrUpperAll_SK3(str) result(strIsUpperAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrUpperAll_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsUpperAll
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isStrUpperAll_SK2(str) result(strIsUpperAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrUpperAll_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsUpperAll
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isStrUpperAll_SK1(str) result(strIsUpperAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrUpperAll_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsUpperAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input string contains **at least** one uppercase English alphabet [ALPHA_UPPER_VEC_SK](@ref ALPHA_UPPER_VEC_SK).<br>
    !>
    !>  \param[in]  str :   The input scalar of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `strIsUpperAny` :   The output scalar or array of the same shape as the input `str` of type `logical` of default kind \LK
    !>                      whose value is `.true.` if at least one element of the input `str` is an uppercase
    !>                      English letter, otherwise it is `.false.`.<br>
    !>
    !>  \interface{isStrUpperAny}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_strASCII, only: isStrUpperAny
    !>      logical(LK) :: strIsUpperAny
    !>
    !>      strIsUpperAny = isStrUpperAny(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The functionality of this interface can be also replicated by [isStrUpper](@ref pm_strASCII::isStrUpper),
    !>  \code{.F90}
    !>
    !>      character(10) :: str
    !>
    !>      isStrUpperAny(str) == any(isStrUpper(str))
    !>
    !>  \endcode
    !>  although this generic interface is potentially faster than [isStrUpper](@ref pm_strASCII::isStrUpper).<br>
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAny](@ref pm_strASCII::isStrUpperAny)<br>
    !>  [isStrLowerAny](@ref pm_strASCII::isStrLowerAny)<br>
    !>  [isStrAlphaNumAll](@ref pm_strASCII::isStrAlphaNumAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrUpperAny}
    !>  \include{lineno} example/pm_strASCII/isStrUpperAny/main.F90
    !>  \compilef{isStrUpperAny}
    !>  \output{isStrUpperAny}
    !>  \include{lineno} example/pm_strASCII/isStrUpperAny/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrUpperAny

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isStrUpperAny_SK5(str) result(strIsUpperAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrUpperAny_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsUpperAny
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isStrUpperAny_SK4(str) result(strIsUpperAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrUpperAny_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsUpperAny
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isStrUpperAny_SK3(str) result(strIsUpperAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrUpperAny_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsUpperAny
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isStrUpperAny_SK2(str) result(strIsUpperAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrUpperAny_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsUpperAny
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isStrUpperAny_SK1(str) result(strIsUpperAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrUpperAny_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsUpperAny
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a `logical` vector whose elements are `.true.` <b>if and only if</b> the corresponding
    !>  characters in the input string belong to the uppercase English alphabet [ALPHA_UPPER_VEC_SK](@ref ALPHA_UPPER_VEC_SK).<br>
    !>
    !>  \param[in]  str :   The input scalar of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `StrIsUpper`    :   The output vector of (rank `1`) of type `logical` of default kind \LK, of the same size as the length of the input `str`
    !>                      whose elements are `.true.` <b>if and only if</b> the corresponding characters in the input string belong to the uppercase English alphabet,
    !>                      otherwise they are `.false.`.<br>
    !>
    !>  \interface{isStrUpper}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_strASCII, only: isStrUpper
    !>      logical(LK) :: StrIsUpper(len(str,IK))
    !>
    !>      StrIsUpper = isStrUpper(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAny](@ref pm_strASCII::isStrUpperAny)<br>
    !>  [isStrLowerAny](@ref pm_strASCII::isStrLowerAny)<br>
    !>  [isStrUpper](@ref pm_strASCII::isStrUpper)<br>
    !>  [isStrLower](@ref pm_strASCII::isStrLower)<br>
    !>  [isStrAlphaNumAll](@ref pm_strASCII::isStrAlphaNumAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrUpper}
    !>  \include{lineno} example/pm_strASCII/isStrUpper/main.F90
    !>  \compilef{isStrUpper}
    !>  \output{isStrUpper}
    !>  \include{lineno} example/pm_strASCII/isStrUpper/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrUpper

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function isStrUpper_SK5(str) result(StrIsUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrUpper_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsUpper(len(str,IK))
    end function
#endif

#if SK4_ENABLED
    pure module function isStrUpper_SK4(str) result(StrIsUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrUpper_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsUpper(len(str,IK))
    end function
#endif

#if SK3_ENABLED
    pure module function isStrUpper_SK3(str) result(StrIsUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrUpper_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsUpper(len(str,IK))
    end function
#endif

#if SK2_ENABLED
    pure module function isStrUpper_SK2(str) result(StrIsUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrUpper_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsUpper(len(str,IK))
    end function
#endif

#if SK1_ENABLED
    pure module function isStrUpper_SK1(str) result(StrIsUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrUpper_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsUpper(len(str,IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input string contains **only** lowercase English alphabets [ALPHA_LOWER_VEC_SK](@ref ALPHA_LOWER_VEC_SK).<br>
    !>
    !>  \param[in]  str :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `strIsLowerAny` :   The output scalar or array of the same shape as the input `str` of type `logical` of default kind \LK
    !>                      whose value is `.true.` corresponding to each element of the input `str` if the element is a lowercase
    !>                      English letter, otherwise it is `.false.`.<br>
    !>
    !>  \interface{isStrLowerAll}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_strASCII, only: isStrLowerAll
    !>      logical(LK) :: strIsLowerAny
    !>
    !>      strIsLowerAny = isStrLowerAll(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The functionality of this interface can be also replicated by [isStrLower](@ref pm_strASCII::isStrLower),
    !>  \code{.F90}
    !>
    !>      character(10) :: str
    !>
    !>      isStrLowerAll(str) == all(isStrLower(str))
    !>
    !>  \endcode
    !>  although this generic interface is potentially faster than [isStrLower](@ref pm_strASCII::isStrLower).<br>
    !>
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>  [isStrAlphaNumAll](@ref pm_strASCII::isStrAlphaNumAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrLowerAll}
    !>  \include{lineno} example/pm_strASCII/isStrLowerAll/main.F90
    !>  \compilef{isStrLowerAll}
    !>  \output{isStrLowerAll}
    !>  \include{lineno} example/pm_strASCII/isStrLowerAll/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrLowerAll

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isStrLowerAll_SK5(str) result(strIsLowerAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrLowerAll_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsLowerAll
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isStrLowerAll_SK4(str) result(strIsLowerAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrLowerAll_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsLowerAll
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isStrLowerAll_SK3(str) result(strIsLowerAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrLowerAll_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsLowerAll
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isStrLowerAll_SK2(str) result(strIsLowerAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrLowerAll_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsLowerAll
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isStrLowerAll_SK1(str) result(strIsLowerAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrLowerAll_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsLowerAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input string contains **at least** one lowercase English alphabet [ALPHA_LOWER_VEC_SK](@ref ALPHA_LOWER_VEC_SK).<br>
    !>
    !>  \param[in]  str :   The input scalar of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `strIsLowerAny` :   The output scalar or array of the same shape as the input `str` of type `logical` of default kind \LK
    !>                      whose value is `.true.` if at least one element of the input `str` is a lowercase
    !>                      English letter, otherwise it is `.false.`.<br>
    !>
    !>  \interface{isStrLowerAny}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_strASCII, only: isStrLowerAny
    !>      logical(LK) :: strIsLowerAny
    !>
    !>      strIsLowerAny = isStrLowerAny(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The functionality of this interface can be also replicated by [isStrLower](@ref pm_strASCII::isStrLower),
    !>  \code{.F90}
    !>
    !>      character(10) :: str
    !>
    !>      isStrLowerAny(str) == any(isStrLower(str))
    !>
    !>  \endcode
    !>  although this generic interface is potentially faster than [isStrLower](@ref pm_strASCII::isStrLower).<br>
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAny](@ref pm_strASCII::isStrUpperAny)<br>
    !>  [isStrLowerAny](@ref pm_strASCII::isStrLowerAny)<br>
    !>  [isStrAlphaNumAll](@ref pm_strASCII::isStrAlphaNumAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrLowerAny}
    !>  \include{lineno} example/pm_strASCII/isStrLowerAny/main.F90
    !>  \compilef{isStrLowerAny}
    !>  \output{isStrLowerAny}
    !>  \include{lineno} example/pm_strASCII/isStrLowerAny/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrLowerAny

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isStrLowerAny_SK5(str) result(strIsLowerAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrLowerAny_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsLowerAny
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isStrLowerAny_SK4(str) result(strIsLowerAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrLowerAny_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsLowerAny
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isStrLowerAny_SK3(str) result(strIsLowerAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrLowerAny_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsLowerAny
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isStrLowerAny_SK2(str) result(strIsLowerAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrLowerAny_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsLowerAny
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isStrLowerAny_SK1(str) result(strIsLowerAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrLowerAny_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsLowerAny
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a `logical` vector whose elements are `.true.` <b>if and only if</b> the corresponding
    !>  characters in the input string belong to the uppercase English alphabet [ALPHA_UPPER_VEC_SK](@ref ALPHA_UPPER_VEC_SK).<br>
    !>
    !>  \param[in]  str :   The input scalar of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `StrIsLower`    :   The output vector of (rank `1`) of type `logical` of default kind \LK, of the same size as the length of the input `str`
    !>                      whose elements are `.true.` <b>if and only if</b> the corresponding characters in the input string belong to the lowercase English alphabet,
    !>                      otherwise they are `.false.`.<br>
    !>
    !>  \interface{isStrLower}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_strASCII, only: isStrLower
    !>      logical(LK) :: CaseIsLower(len(str,IK))
    !>
    !>      StrIsLower = isStrLower(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpper](@ref pm_strASCII::isStrUpper)<br>
    !>  [isStrLower](@ref pm_strASCII::isStrLower)<br>
    !>  [isStrAlphaNumAll](@ref pm_strASCII::isStrAlphaNumAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrLower}
    !>  \include{lineno} example/pm_strASCII/isStrLower/main.F90
    !>  \compilef{isStrLower}
    !>  \output{isStrLower}
    !>  \include{lineno} example/pm_strASCII/isStrLower/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrLower

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function isStrLower_SK5(str) result(StrIsLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrLower_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsLower(len(str,IK))
    end function
#endif

#if SK4_ENABLED
    pure module function isStrLower_SK4(str) result(StrIsLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrLower_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsLower(len(str,IK))
    end function
#endif

#if SK3_ENABLED
    pure module function isStrLower_SK3(str) result(StrIsLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrLower_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsLower(len(str,IK))
    end function
#endif

#if SK2_ENABLED
    pure module function isStrLower_SK2(str) result(StrIsLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrLower_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsLower(len(str,IK))
    end function
#endif

#if SK1_ENABLED
    pure module function isStrLower_SK1(str) result(StrIsLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrLower_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsLower(len(str,IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input (single) character is a digit or
    !>  an uppercase or lowercase English alphabet [ALPHA_VEC_SK](@ref ALPHA_VEC_SK).<br>
    !>
    !>  \param[in]  chr     :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL of `len = 1`.<br>
    !>
    !>  \return
    !>  `charIsAlphaNum`    :   The output scalar or array of the same shape as the input `chr` of type `logical` of default kind \LK
    !>                          whose value is `.true.` corresponding to each element of the input `chr` if the element belongs to
    !>                          the set of English digits or English alphabets (lowercase or uppercase), otherwise it is `.false.`.<br>
    !>
    !>  \interface{isCharAlphaNum}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_strASCII, only: isCharAlphaNum
    !>      logical(LK) :: charIsAlphaNum
    !>
    !>      charIsAlphaNum = isCharAlphaNum(chr)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>  [isCharAlphaNum](@ref pm_strASCII::isCharAlphaNum)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isCharAlphaNum}
    !>  \include{lineno} example/pm_strASCII/isCharAlphaNum/main.F90
    !>  \compilef{isCharAlphaNum}
    !>  \output{isCharAlphaNum}
    !>  \include{lineno} example/pm_strASCII/isCharAlphaNum/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isCharAlphaNum

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isCharAlphaNum_SK5(chr) result(charIsAlphaNum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharAlphaNum_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsAlphaNum
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isCharAlphaNum_SK4(chr) result(charIsAlphaNum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharAlphaNum_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsAlphaNum
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isCharAlphaNum_SK3(chr) result(charIsAlphaNum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharAlphaNum_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsAlphaNum
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isCharAlphaNum_SK2(chr) result(charIsAlphaNum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharAlphaNum_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsAlphaNum
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isCharAlphaNum_SK1(chr) result(charIsAlphaNum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharAlphaNum_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsAlphaNum
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input string is **all** alphanumeric,
    !>  containing only digits or the English alphabet [ALPHANUM_VEC_SK](@ref ALPHANUM_VEC_SK).<br>
    !>
    !>  \param[in]  str     :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `strIsAlphaNumAll`  :   The output scalar or array of the same shape as the input `str` of type `logical` of default kind \LK
    !>                          whose value is `.true.` corresponding to each element of the input `str` if the scalar string is all digits,
    !>                          or English alphabets, otherwise it is `.false.`.<br>
    !>
    !>  \interface{isStrAlphaNumAll}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_strASCII, only: isStrAlphaNumAll
    !>      logical(LK) :: strIsAlphaNumAll
    !>
    !>      strIsAlphaNumAll = isStrAlphaNumAll(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The functionality of this interface can be also replicated by [isStrAlphaNum](@ref pm_strASCII::isStrAlphaNum),
    !>  \code{.F90}
    !>
    !>      character(10) :: str
    !>
    !>      isStrAlphaNumAll(str) == all(isStrAlphaNum(str))
    !>
    !>  \endcode
    !>  although this generic interface is potentially faster than [isStrAlphaNum](@ref pm_strASCII::isStrAlphaNum).<br>
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>  [isStrAlphaNumAll](@ref pm_strASCII::isStrAlphaNumAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrAlphaNumAll}
    !>  \include{lineno} example/pm_strASCII/isStrAlphaNumAll/main.F90
    !>  \compilef{isStrAlphaNumAll}
    !>  \output{isStrAlphaNumAll}
    !>  \include{lineno} example/pm_strASCII/isStrAlphaNumAll/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrAlphaNumAll

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isStrAlphaNumAll_SK5(str) result(strIsAlphaNumAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaNumAll_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaNumAll
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isStrAlphaNumAll_SK4(str) result(strIsAlphaNumAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaNumAll_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaNumAll
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isStrAlphaNumAll_SK3(str) result(strIsAlphaNumAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaNumAll_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaNumAll
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isStrAlphaNumAll_SK2(str) result(strIsAlphaNumAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaNumAll_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaNumAll
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isStrAlphaNumAll_SK1(str) result(strIsAlphaNumAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaNumAll_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaNumAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if **any** characters of the input string are alphanumeric,
    !>  that is, only digits or the English alphabet [ALPHANUM_VEC_SK](@ref ALPHANUM_VEC_SK).<br>
    !>
    !>  \param[in]  str     :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `strIsAlphaNumAny`  :   The output scalar or array of the same shape as the input `str` of type `logical` of default kind \LK.<br>
    !>                          If it `.true.` <b>if and only if</b> any character of the input string `str` is digits or English alphabet.<br>
    !>
    !>  \interface{isStrAlphaNumAny}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_strASCII, only: isStrAlphaNumAny
    !>      logical(LK) :: strIsAlphaNumAny
    !>
    !>      strIsAlphaNumAny = isStrAlphaNumAny(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The functionality of this interface can be also replicated by [isStrAlphaNum](@ref pm_strASCII::isStrAlphaNum),
    !>  \code{.F90}
    !>
    !>      character(10) :: str
    !>
    !>      isStrAlphaNumAny(str) == any(isStrAlphaNum(str))
    !>
    !>  \endcode
    !>  although this generic interface is potentially faster than [isStrAlphaNum](@ref pm_strASCII::isStrAlphaNum).<br>
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>  [isStrAlphaNumAny](@ref pm_strASCII::isStrAlphaNumAny)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrAlphaNumAny}
    !>  \include{lineno} example/pm_strASCII/isStrAlphaNumAny/main.F90
    !>  \compilef{isStrAlphaNumAny}
    !>  \output{isStrAlphaNumAny}
    !>  \include{lineno} example/pm_strASCII/isStrAlphaNumAny/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrAlphaNumAny

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isStrAlphaNumAny_SK5(str) result(strIsAlphaNumAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaNumAny_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaNumAny
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isStrAlphaNumAny_SK4(str) result(strIsAlphaNumAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaNumAny_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaNumAny
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isStrAlphaNumAny_SK3(str) result(strIsAlphaNumAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaNumAny_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaNumAny
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isStrAlphaNumAny_SK2(str) result(strIsAlphaNumAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaNumAny_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaNumAny
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isStrAlphaNumAny_SK1(str) result(strIsAlphaNumAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaNumAny_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaNumAny
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a `logical` vector whose elements are `.true.` <b>if and only if</b> the corresponding
    !>  characters in the input string belong to the uppercase English alphabet [ALPHA_UPPER_VEC_SK](@ref ALPHA_UPPER_VEC_SK).<br>
    !>
    !>  \param[in]  str :   The input scalar of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `StrIsLower`    :   The output vector of (rank `1`) of type `logical` of default kind \LK, of the same size as the length of the input `str`
    !>                      whose elements are `.true.` <b>if and only if</b> the corresponding characters in the input string belong to the lowercase English alphabet,
    !>                      otherwise they are `.false.`.<br>
    !>
    !>  \interface{isStrLower}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_strASCII, only: isStrLower
    !>      logical(LK) :: CaseIsLower(len(str,IK))
    !>
    !>      StrIsLower = isStrLower(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpper](@ref pm_strASCII::isStrUpper)<br>
    !>  [isStrLower](@ref pm_strASCII::isStrLower)<br>
    !>  [isStrAlphaNumAll](@ref pm_strASCII::isStrAlphaNumAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrLower}
    !>  \include{lineno} example/pm_strASCII/isStrLower/main.F90
    !>  \compilef{isStrLower}
    !>  \output{isStrLower}
    !>  \include{lineno} example/pm_strASCII/isStrLower/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrAlphaNum

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function isStrAlphaNum_SK5(str) result(StrIsAlphaNum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaNum_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsAlphaNum(len(str,IK))
    end function
#endif

#if SK4_ENABLED
    pure module function isStrAlphaNum_SK4(str) result(StrIsAlphaNum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaNum_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsAlphaNum(len(str,IK))
    end function
#endif

#if SK3_ENABLED
    pure module function isStrAlphaNum_SK3(str) result(StrIsAlphaNum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaNum_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsAlphaNum(len(str,IK))
    end function
#endif

#if SK2_ENABLED
    pure module function isStrAlphaNum_SK2(str) result(StrIsAlphaNum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaNum_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsAlphaNum(len(str,IK))
    end function
#endif

#if SK1_ENABLED
    pure module function isStrAlphaNum_SK1(str) result(StrIsAlphaNum)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaNum_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsAlphaNum(len(str,IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input (single) character is an
    !>  uppercase or lowercase English alphabet [ALPHA_VEC_SK](@ref ALPHA_VEC_SK).<br>
    !>
    !>  \param[in]  chr     :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL of `len = 1`.<br>
    !>
    !>  \return
    !>  `charIsAlpha`       :   The output scalar or array of the same shape as the input `chr` of type `logical` of default kind \LK
    !>                          whose value is `.true.` corresponding to each element of the input `chr` if the element belongs to
    !>                          the English alphabets (lowercase or uppercase), otherwise it is `.false.`.<br>
    !>
    !>  \interface{isCharAlpha}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_strASCII, only: isCharAlpha
    !>      logical(LK) :: charIsAlpha
    !>
    !>      charIsAlpha = isCharAlpha(chr)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>  [isCharAlpha](@ref pm_strASCII::isCharAlpha)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isCharAlpha}
    !>  \include{lineno} example/pm_strASCII/isCharAlpha/main.F90
    !>  \compilef{isCharAlpha}
    !>  \output{isCharAlpha}
    !>  \include{lineno} example/pm_strASCII/isCharAlpha/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isCharAlpha

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isCharAlpha_SK5(chr) result(charIsAlpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharAlpha_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsAlpha
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isCharAlpha_SK4(chr) result(charIsAlpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharAlpha_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsAlpha
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isCharAlpha_SK3(chr) result(charIsAlpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharAlpha_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsAlpha
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isCharAlpha_SK2(chr) result(charIsAlpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharAlpha_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsAlpha
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isCharAlpha_SK1(chr) result(charIsAlpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isCharAlpha_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(1,SKG)            , intent(in)                :: chr
        logical(LK)                                             :: charIsAlpha
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input string is **all** alphabetic,
    !>  containing only the English alphabet [ALPHA_VEC_SK](@ref ALPHA_VEC_SK).<br>
    !>
    !>  \param[in]  str :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `strIsAlphaAll` :   The output scalar or array of the same shape as the input `str` of type `logical` of default kind \LK
    !>                      whose value is `.true.` corresponding to each element of the input `str` if the scalar string is all
    !>                      English alphabets, otherwise it is `.false.`.<br>
    !>
    !>  \interface{isStrAlphaAll}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_strASCII, only: isStrAlphaAll
    !>      logical(LK) :: strIsAlphaAll
    !>
    !>      strIsAlphaAll = isStrAlphaAll(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The functionality of this interface can be also replicated by [isStrAlpha](@ref pm_strASCII::isStrAlpha),
    !>  \code{.F90}
    !>
    !>      character(10) :: str
    !>
    !>      isStrAlphaAll(str) == all(isStrAlpha(str))
    !>
    !>  \endcode
    !>  although this generic interface is potentially faster than [isStrAlpha](@ref pm_strASCII::isStrAlpha).<br>
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>  [isStrAlphaAll](@ref pm_strASCII::isStrAlphaAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrAlphaAll}
    !>  \include{lineno} example/pm_strASCII/isStrAlphaAll/main.F90
    !>  \compilef{isStrAlphaAll}
    !>  \output{isStrAlphaAll}
    !>  \include{lineno} example/pm_strASCII/isStrAlphaAll/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrAlphaAll

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isStrAlphaAll_SK5(str) result(strIsAlphaAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaAll_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaAll
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isStrAlphaAll_SK4(str) result(strIsAlphaAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaAll_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaAll
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isStrAlphaAll_SK3(str) result(strIsAlphaAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaAll_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaAll
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isStrAlphaAll_SK2(str) result(strIsAlphaAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaAll_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaAll
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isStrAlphaAll_SK1(str) result(strIsAlphaAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaAll_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if **any** characters of the input string are alphabetic,
    !>  that is, only the English alphabet [ALPHA_VEC_SK](@ref ALPHA_VEC_SK).<br>
    !>
    !>  \param[in]  str     :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `strIsAlphaAny`  :   The output scalar or array of the same shape as the input `str` of type `logical` of default kind \LK.<br>
    !>                          If it `.true.` <b>if and only if</b> any character of the input string `str` is English alphabet.<br>
    !>
    !>  \interface{isStrAlphaAny}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_strASCII, only: isStrAlphaAny
    !>      logical(LK) :: strIsAlphaAny
    !>
    !>      strIsAlphaAny = isStrAlphaAny(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  The functionality of this interface can be also replicated by [isStrAlpha](@ref pm_strASCII::isStrAlpha),
    !>  \code{.F90}
    !>
    !>      character(10) :: str
    !>
    !>      isStrAlphaAny(str) == any(isStrAlpha(str))
    !>
    !>  \endcode
    !>  although this generic interface is potentially faster than [isStrAlpha](@ref pm_strASCII::isStrAlpha).<br>
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>  [isStrAlphaAny](@ref pm_strASCII::isStrAlphaAny)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrAlphaAny}
    !>  \include{lineno} example/pm_strASCII/isStrAlphaAny/main.F90
    !>  \compilef{isStrAlphaAny}
    !>  \output{isStrAlphaAny}
    !>  \include{lineno} example/pm_strASCII/isStrAlphaAny/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrAlphaAny

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function isStrAlphaAny_SK5(str) result(strIsAlphaAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaAny_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaAny
    end function
#endif

#if SK4_ENABLED
    pure elemental module function isStrAlphaAny_SK4(str) result(strIsAlphaAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaAny_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaAny
    end function
#endif

#if SK3_ENABLED
    pure elemental module function isStrAlphaAny_SK3(str) result(strIsAlphaAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaAny_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaAny
    end function
#endif

#if SK2_ENABLED
    pure elemental module function isStrAlphaAny_SK2(str) result(strIsAlphaAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaAny_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaAny
    end function
#endif

#if SK1_ENABLED
    pure elemental module function isStrAlphaAny_SK1(str) result(strIsAlphaAny)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlphaAny_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: strIsAlphaAny
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a `logical` vector whose elements are `.true.` <b>if and only if</b> the corresponding
    !>  characters in the input string belong to the uppercase English alphabet [ALPHA_UPPER_VEC_SK](@ref ALPHA_UPPER_VEC_SK).<br>
    !>
    !>  \param[in]  str :   The input scalar of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `StrIsLower`    :   The output vector of (rank `1`) of type `logical` of default kind \LK, of the same size as the length of the input `str`
    !>                      whose elements are `.true.` <b>if and only if</b> the corresponding characters in the input string belong to the lowercase English alphabet,
    !>                      otherwise they are `.false.`.<br>
    !>
    !>  \interface{isStrLower}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_strASCII, only: isStrLower
    !>      logical(LK) :: CaseIsLower(len(str,IK))
    !>
    !>      StrIsLower = isStrLower(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \see
    !>  [isStrReal](@ref pm_strASCII::isStrReal)<br>
    !>  [isCharDigit](@ref pm_strASCII::isCharDigit)<br>
    !>  [isStrNumber](@ref pm_strASCII::isStrNumber)<br>
    !>  [isStrDigitAll](@ref pm_strASCII::isStrDigitAll)<br>
    !>  [isStrInteger](@ref pm_strASCII::isStrInteger)<br>
    !>  [isStrComplex](@ref pm_strASCII::isStrComplex)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpper](@ref pm_strASCII::isStrUpper)<br>
    !>  [isStrLower](@ref pm_strASCII::isStrLower)<br>
    !>  [isStrAlphaAll](@ref pm_strASCII::isStrAlphaAll)<br>
    !>  [isStrAlpha](@ref pm_strASCII::isStrAlpha)<br>
    !>
    !>  \example{isStrLower}
    !>  \include{lineno} example/pm_strASCII/isStrLower/main.F90
    !>  \compilef{isStrLower}
    !>  \output{isStrLower}
    !>  \include{lineno} example/pm_strASCII/isStrLower/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isStrAlpha

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function isStrAlpha_SK5(str) result(StrIsAlpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlpha_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsAlpha(len(str,IK))
    end function
#endif

#if SK4_ENABLED
    pure module function isStrAlpha_SK4(str) result(StrIsAlpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlpha_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsAlpha(len(str,IK))
    end function
#endif

#if SK3_ENABLED
    pure module function isStrAlpha_SK3(str) result(StrIsAlpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlpha_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsAlpha(len(str,IK))
    end function
#endif

#if SK2_ENABLED
    pure module function isStrAlpha_SK2(str) result(StrIsAlpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlpha_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsAlpha(len(str,IK))
    end function
#endif

#if SK1_ENABLED
    pure module function isStrAlpha_SK1(str) result(StrIsAlpha)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isStrAlpha_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        logical(LK)                                             :: StrIsAlpha(len(str,IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the input character where the lowercase English alphabet is converted to uppercase letter.<br>
    !>
    !>  \param[in]  chr :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL of length `1`.<br>
    !>
    !>  \return
    !>  `chrUpper`      :   The output scalar or array of the same type, kind, rank, and `len` type parameter as the input `chr`
    !>                      containing the input character where the lowercase English alphabet is converted to uppercase letter.<br>
    !>
    !>  \interface{getCharUpper}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: getCharUpper
    !>      character(len(chr),kind(chr)):: chrUpper(size(chr))
    !>
    !>      chrUpper = getCharUpper(chr)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getCharUpper](@ref pm_strASCII::getCharUpper)<br>
    !>  [setCharUpper](@ref pm_strASCII::setCharUpper)<br>
    !>  [getCharLower](@ref pm_strASCII::getCharLower)<br>
    !>  [setCharLower](@ref pm_strASCII::setCharLower)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>
    !>  \example{getCharUpper}
    !>  \include{lineno} example/pm_strASCII/getCharUpper/main.F90
    !>  \compilef{getCharUpper}
    !>  \output{getCharUpper}
    !>  \include{lineno} example/pm_strASCII/getCharUpper/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getCharUpper

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function getCharUpper_SK5(chr) result(chrUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharUpper_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(1,SKG)            , intent(in)                :: chr
        character(1,SKG)                                        :: chrUpper
    end function
#endif

#if SK4_ENABLED
    pure elemental module function getCharUpper_SK4(chr) result(chrUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharUpper_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(1,SKG)            , intent(in)                :: chr
        character(1,SKG)                                        :: chrUpper
    end function
#endif

#if SK3_ENABLED
    pure elemental module function getCharUpper_SK3(chr) result(chrUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharUpper_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(1,SKG)            , intent(in)                :: chr
        character(1,SKG)                                        :: chrUpper
    end function
#endif

#if SK2_ENABLED
    pure elemental module function getCharUpper_SK2(chr) result(chrUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharUpper_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(1,SKG)            , intent(in)                :: chr
        character(1,SKG)                                        :: chrUpper
    end function
#endif

#if SK1_ENABLED
    pure elemental module function getCharUpper_SK1(chr) result(chrUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharUpper_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(1,SKG)            , intent(in)                :: chr
        character(1,SKG)                                        :: chrUpper
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Replace any lowercase English alphabet in the input character with the corresponding uppercase English letter.<br>
    !>
    !>  \param[inout]   chr :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL of length `1`.<br>
    !>                          &nbsp; On output, any lowercase English letter will be replaced with the corresponding uppercase letter.<br>
    !>
    !>  \interface{setCharUpper}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: setCharUpper
    !>
    !>      call setCharUpper(chr)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setCharUpper](@ref pm_strASCII::setCharUpper)<br>
    !>  [setCharUpper](@ref pm_strASCII::setCharUpper)<br>
    !>  [getCharLower](@ref pm_strASCII::getCharLower)<br>
    !>  [setCharLower](@ref pm_strASCII::setCharLower)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>
    !>  \example{setCharUpper}
    !>  \include{lineno} example/pm_strASCII/setCharUpper/main.F90
    !>  \compilef{setCharUpper}
    !>  \output{setCharUpper}
    !>  \include{lineno} example/pm_strASCII/setCharUpper/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setCharUpper

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module subroutine setCharUpper_SK5(chr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCharUpper_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(1,SKG)            , intent(inout)             :: chr
    end subroutine
#endif

#if SK4_ENABLED
    pure elemental module subroutine setCharUpper_SK4(chr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCharUpper_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(1,SKG)            , intent(inout)             :: chr
    end subroutine
#endif

#if SK3_ENABLED
    pure elemental module subroutine setCharUpper_SK3(chr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCharUpper_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(1,SKG)            , intent(inout)             :: chr
    end subroutine
#endif

#if SK2_ENABLED
    pure elemental module subroutine setCharUpper_SK2(chr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCharUpper_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(1,SKG)            , intent(inout)             :: chr
    end subroutine
#endif

#if SK1_ENABLED
    pure elemental module subroutine setCharUpper_SK1(chr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCharUpper_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(1,SKG)            , intent(inout)             :: chr
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the input character where the uppercase English alphabet is converted to lowercase letter.<br>
    !>
    !>  \param[in]  chr :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL of length `1`.<br>
    !>
    !>  \return
    !>  `chrUpper`      :   The output scalar or array of the same type, kind, rank, and `len` type parameter as the input `chr`
    !>                      containing the input character where the uppercase English alphabet is converted to uppercase letter.<br>
    !>
    !>  \interface{getCharLower}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: getCharLower
    !>      character(len(chr),kind(chr)):: chrLower(size(chr))
    !>
    !>      chrLower = getCharLower(chr)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getCharLower](@ref pm_strASCII::getCharLower)<br>
    !>  [setCharLower](@ref pm_strASCII::setCharLower)<br>
    !>  [getCharLower](@ref pm_strASCII::getCharLower)<br>
    !>  [setCharLower](@ref pm_strASCII::setCharLower)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>
    !>  \example{getCharLower}
    !>  \include{lineno} example/pm_strASCII/getCharLower/main.F90
    !>  \compilef{getCharLower}
    !>  \output{getCharLower}
    !>  \include{lineno} example/pm_strASCII/getCharLower/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getCharLower

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function getCharLower_SK5(chr) result(chrLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharLower_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(1,SKG)            , intent(in)                :: chr
        character(1,SKG)                                        :: chrLower
    end function
#endif

#if SK4_ENABLED
    pure elemental module function getCharLower_SK4(chr) result(chrLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharLower_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(1,SKG)            , intent(in)                :: chr
        character(1,SKG)                                        :: chrLower
    end function
#endif

#if SK3_ENABLED
    pure elemental module function getCharLower_SK3(chr) result(chrLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharLower_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(1,SKG)            , intent(in)                :: chr
        character(1,SKG)                                        :: chrLower
    end function
#endif

#if SK2_ENABLED
    pure elemental module function getCharLower_SK2(chr) result(chrLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharLower_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(1,SKG)            , intent(in)                :: chr
        character(1,SKG)                                        :: chrLower
    end function
#endif

#if SK1_ENABLED
    pure elemental module function getCharLower_SK1(chr) result(chrLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCharLower_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(1,SKG)            , intent(in)                :: chr
        character(1,SKG)                                        :: chrLower
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Replace any uppercase English alphabet in the input character with the corresponding lowercase English letter.<br>
    !>
    !>  \param[inout]   chr :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL of length `1`.<br>
    !>                          &nbsp; On output, any uppercase English letter will be replaced with the corresponding lowercase letter.<br>
    !>
    !>  \interface{setCharLower}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: setCharLower
    !>
    !>      call setCharLower(chr)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setCharLower](@ref pm_strASCII::setCharLower)<br>
    !>  [setCharLower](@ref pm_strASCII::setCharLower)<br>
    !>  [getCharLower](@ref pm_strASCII::getCharLower)<br>
    !>  [setCharLower](@ref pm_strASCII::setCharLower)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>
    !>  \example{setCharLower}
    !>  \include{lineno} example/pm_strASCII/setCharLower/main.F90
    !>  \compilef{setCharLower}
    !>  \output{setCharLower}
    !>  \include{lineno} example/pm_strASCII/setCharLower/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setCharLower

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module subroutine setCharLower_SK5(chr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCharLower_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(1,SKG)            , intent(inout)             :: chr
    end subroutine
#endif

#if SK4_ENABLED
    pure elemental module subroutine setCharLower_SK4(chr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCharLower_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(1,SKG)            , intent(inout)             :: chr
    end subroutine
#endif

#if SK3_ENABLED
    pure elemental module subroutine setCharLower_SK3(chr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCharLower_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(1,SKG)            , intent(inout)             :: chr
    end subroutine
#endif

#if SK2_ENABLED
    pure elemental module subroutine setCharLower_SK2(chr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCharLower_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(1,SKG)            , intent(inout)             :: chr
    end subroutine
#endif

#if SK1_ENABLED
    pure elemental module subroutine setCharLower_SK1(chr)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCharLower_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(1,SKG)            , intent(inout)             :: chr
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the input string where the lowercase English alphabets are all converted to uppercase letters.<br>
    !>
    !>  \param[in]  str :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `strUpper`      :   The output scalar or array of the same type, kind, and rank as the input `str` containing
    !>                      the input string where the lowercase English alphabets are all converted to uppercase letters.<br>
    !>
    !>  \interface{getStrUpper}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: getStrUpper
    !>      character(len(str,IK),kind(str)) :: strUpper(size(str))
    !>
    !>      strUpper = getStrUpper(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getStrUpper](@ref pm_strASCII::getStrUpper)<br>
    !>  [setStrUpper](@ref pm_strASCII::setStrUpper)<br>
    !>  [getStrLower](@ref pm_strASCII::getStrLower)<br>
    !>  [setStrLower](@ref pm_strASCII::setStrLower)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>
    !>  \example{getStrUpper}
    !>  \include{lineno} example/pm_strASCII/getStrUpper/main.F90
    !>  \compilef{getStrUpper}
    !>  \output{getStrUpper}
    !>  \include{lineno} example/pm_strASCII/getStrUpper/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getStrUpper

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function getStrUpper_SK5(str) result(strUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrUpper_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        character(len(str,IK),SKG)                              :: strUpper
    end function
#endif

#if SK4_ENABLED
    pure elemental module function getStrUpper_SK4(str) result(strUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrUpper_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        character(len(str,IK),SKG)                              :: strUpper
    end function
#endif

#if SK3_ENABLED
    pure elemental module function getStrUpper_SK3(str) result(strUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrUpper_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        character(len(str,IK),SKG)                              :: strUpper
    end function
#endif

#if SK2_ENABLED
    pure elemental module function getStrUpper_SK2(str) result(strUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrUpper_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        character(len(str,IK),SKG)                              :: strUpper
    end function
#endif

#if SK1_ENABLED
    pure elemental module function getStrUpper_SK1(str) result(strUpper)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrUpper_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        character(len(str,IK),SKG)                              :: strUpper
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Replace all lowercase English alphabets in the input string with the corresponding uppercase English letters.<br>
    !>
    !>  \param[inout]   str :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL.<br>
    !>                          On output, all instances of lowercase English letters will be replaced with the corresponding uppercase letters.<br>
    !>
    !>  \interface{setStrUpper}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: setStrUpper
    !>
    !>      call setStrUpper(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setStrUpper](@ref pm_strASCII::setStrUpper)<br>
    !>  [setStrUpper](@ref pm_strASCII::setStrUpper)<br>
    !>  [getStrLower](@ref pm_strASCII::getStrLower)<br>
    !>  [setStrLower](@ref pm_strASCII::setStrLower)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>
    !>  \example{setStrUpper}
    !>  \include{lineno} example/pm_strASCII/setStrUpper/main.F90
    !>  \compilef{setStrUpper}
    !>  \output{setStrUpper}
    !>  \include{lineno} example/pm_strASCII/setStrUpper/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setStrUpper

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module subroutine setStrUpper_SK5(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStrUpper_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout)             :: str
    end subroutine
#endif

#if SK4_ENABLED
    pure elemental module subroutine setStrUpper_SK4(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStrUpper_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout)             :: str
    end subroutine
#endif

#if SK3_ENABLED
    pure elemental module subroutine setStrUpper_SK3(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStrUpper_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout)             :: str
    end subroutine
#endif

#if SK2_ENABLED
    pure elemental module subroutine setStrUpper_SK2(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStrUpper_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout)             :: str
    end subroutine
#endif

#if SK1_ENABLED
    pure elemental module subroutine setStrUpper_SK1(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStrUpper_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout)             :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the input string where the uppercase English alphabets are all converted to lowercase letters.<br>
    !>
    !>  \param[in]  str :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL.<br>
    !>
    !>  \return
    !>  `strLower`      :   The output scalar or array of the same type, kind, and rank as the input `str` containing
    !>                      the input string where the uppercase English alphabets are all converted to lowercase letters.<br>
    !>
    !>  \interface{getStrLower}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: getStrLower
    !>      character(len(str,IK),kind(str)) :: strLower(size(str))
    !>
    !>      strLower = getStrLower(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getStrLower](@ref pm_strASCII::getStrLower)<br>
    !>  [setStrLower](@ref pm_strASCII::setStrLower)<br>
    !>  [getStrLower](@ref pm_strASCII::getStrLower)<br>
    !>  [setStrLower](@ref pm_strASCII::setStrLower)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>
    !>  \example{getStrLower}
    !>  \include{lineno} example/pm_strASCII/getStrLower/main.F90
    !>  \compilef{getStrLower}
    !>  \output{getStrLower}
    !>  \include{lineno} example/pm_strASCII/getStrLower/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getStrLower

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module function getStrLower_SK5(str) result(strLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrLower_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        character(len(str,IK),SKG)                              :: strLower
    end function
#endif

#if SK4_ENABLED
    pure elemental module function getStrLower_SK4(str) result(strLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrLower_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        character(len(str,IK),SKG)                              :: strLower
    end function
#endif

#if SK3_ENABLED
    pure elemental module function getStrLower_SK3(str) result(strLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrLower_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        character(len(str,IK),SKG)                              :: strLower
    end function
#endif

#if SK2_ENABLED
    pure elemental module function getStrLower_SK2(str) result(strLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrLower_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        character(len(str,IK),SKG)                              :: strLower
    end function
#endif

#if SK1_ENABLED
    pure elemental module function getStrLower_SK1(str) result(strLower)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrLower_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        character(len(str,IK),SKG)                              :: strLower
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Replace all uppercase English alphabets in the input string with the corresponding lowercase English letters.<br>
    !>
    !>  \param[inout]   str :   The input scalar or array of arbitrary rank of type `character` of kind \SKALL.<br>
    !>                          On output, all instances of uppercase English letters will be replaced with the corresponding lowercase letters.<br>
    !>
    !>  \interface{setStrLower}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: setStrLower
    !>
    !>      call setStrLower(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [setStrLower](@ref pm_strASCII::setStrLower)<br>
    !>  [setStrLower](@ref pm_strASCII::setStrLower)<br>
    !>  [getStrLower](@ref pm_strASCII::getStrLower)<br>
    !>  [setStrLower](@ref pm_strASCII::setStrLower)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>
    !>  \example{setStrLower}
    !>  \include{lineno} example/pm_strASCII/setStrLower/main.F90
    !>  \compilef{setStrLower}
    !>  \output{setStrLower}
    !>  \include{lineno} example/pm_strASCII/setStrLower/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setStrLower

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure elemental module subroutine setStrLower_SK5(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStrLower_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout)             :: str
    end subroutine
#endif

#if SK4_ENABLED
    pure elemental module subroutine setStrLower_SK4(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStrLower_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout)             :: str
    end subroutine
#endif

#if SK3_ENABLED
    pure elemental module subroutine setStrLower_SK3(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStrLower_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout)             :: str
    end subroutine
#endif

#if SK2_ENABLED
    pure elemental module subroutine setStrLower_SK2(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStrLower_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout)             :: str
    end subroutine
#endif

#if SK1_ENABLED
    pure elemental module subroutine setStrLower_SK1(str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStrLower_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout)             :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the input string quoted with double-quotation marks where all
    !>  instances of double-quotations within the input string are escaped in the Fortran-style.<br>
    !>
    !>  \details
    !>  The Fortran-style quotation escaping is done by duplicating each quotation mark within the string.<br>
    !>
    !>  \param[in]  str :   The input scalar of type `character` of kind \SKALL containing the string to be quoted.<br>
    !>
    !>  \return
    !>  `strQuoted`     :   The output `allocatable` scalar or array of the same type, kind, and rank as the input `str` containing
    !>                      the input string quoted with double-quotation marks where all instances of double-quotes are escaped.<br>
    !>
    !>  \interface{getStrQuoted}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: getStrQuoted
    !>      character(:,kind(str)), allocatable :: strQuoted(size(str))
    !>
    !>      strQuoted = getStrQuoted(str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \see
    !>  [getStrLower](@ref pm_strASCII::getStrLower)<br>
    !>  [setStrLower](@ref pm_strASCII::setStrLower)<br>
    !>  [getStrLower](@ref pm_strASCII::getStrLower)<br>
    !>  [setStrLower](@ref pm_strASCII::setStrLower)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>
    !>  \example{getStrQuoted}
    !>  \include{lineno} example/pm_strASCII/getStrQuoted/main.F90
    !>  \compilef{getStrQuoted}
    !>  \output{getStrQuoted}
    !>  \include{lineno} example/pm_strASCII/getStrQuoted/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)<br>
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to include an optional argument `delim` of type
    !>  `character(1,SKG)` that contains the user-specified quotation mark (other than double-quote).<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getStrQuoted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function getStrQuoted_SK5(str) result(strQuoted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrQuoted_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        character(:,SKG)            , allocatable               :: strQuoted
    end function
#endif

#if SK4_ENABLED
    pure module function getStrQuoted_SK4(str) result(strQuoted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrQuoted_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        character(:,SKG)            , allocatable               :: strQuoted
    end function
#endif

#if SK3_ENABLED
    pure module function getStrQuoted_SK3(str) result(strQuoted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrQuoted_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        character(:,SKG)            , allocatable               :: strQuoted
    end function
#endif

#if SK2_ENABLED
    pure module function getStrQuoted_SK2(str) result(strQuoted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrQuoted_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        character(:,SKG)            , allocatable               :: strQuoted
    end function
#endif

#if SK1_ENABLED
    pure module function getStrQuoted_SK1(str) result(strQuoted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getStrQuoted_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        character(:,SKG)            , allocatable               :: strQuoted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the input string quoted with double-quotation marks where all
    !>  instances of double-quotations within the input string are escaped in the Fortran-style.<br>
    !>
    !>  \details
    !>  The Fortran-style quotation escaping is done by duplicating each quotation mark within the string.<br>
    !>
    !>  \param[out] strQuoted   :   The output `allocatable` scalar `character` of kind \SKALL containing the input string
    !>                              quoted with double-quotation marks where all instances of double-quotes are escaped.<br>
    !>  \param[in]  str         :   The input scalar of the same type and kind as `strQuoted`, of arbitrary length type parameter containing the string to be quoted.<br>
    !>
    !>  \interface{setStrQuoted}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: setStrQuoted
    !>
    !>      call setStrQuoted(strQuoted, str)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \see
    !>  [setStrLower](@ref pm_strASCII::setStrLower)<br>
    !>  [setStrLower](@ref pm_strASCII::setStrLower)<br>
    !>  [getStrLower](@ref pm_strASCII::getStrLower)<br>
    !>  [setStrLower](@ref pm_strASCII::setStrLower)<br>
    !>  [isCharUpper](@ref pm_strASCII::isCharUpper)<br>
    !>  [isCharLower](@ref pm_strASCII::isCharLower)<br>
    !>  [isStrUpperAll](@ref pm_strASCII::isStrUpperAll)<br>
    !>  [isStrLowerAll](@ref pm_strASCII::isStrLowerAll)<br>
    !>
    !>  \example{setStrQuoted}
    !>  \include{lineno} example/pm_strASCII/setStrQuoted/main.F90
    !>  \compilef{setStrQuoted}
    !>  \output{setStrQuoted}
    !>  \include{lineno} example/pm_strASCII/setStrQuoted/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \todo
    !>  \pmed This generic interface can be extended to include an optional argument `delim` of type
    !>  `character(1,SKG)` that contains the user-specified quotation mark (other than double-quote).<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setStrQuoted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setStrQuoted_SK5(strQuoted, str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStrQuoted_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(:,SKG)    , intent(out)   , allocatable       :: strQuoted
        character(*,SKG)    , intent(in)                        :: str
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setStrQuoted_SK4(strQuoted, str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStrQuoted_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(:,SKG)    , intent(out)   , allocatable       :: strQuoted
        character(*,SKG)    , intent(in)                        :: str
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setStrQuoted_SK3(strQuoted, str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStrQuoted_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(:,SKG)    , intent(out)   , allocatable       :: strQuoted
        character(*,SKG)    , intent(in)                        :: str
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setStrQuoted_SK2(strQuoted, str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStrQuoted_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(:,SKG)    , intent(out)   , allocatable       :: strQuoted
        character(*,SKG)    , intent(in)                        :: str
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setStrQuoted_SK1(strQuoted, str)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStrQuoted_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(:,SKG)    , intent(out)   , allocatable       :: strQuoted
        character(*,SKG)    , intent(in)                        :: str
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the input C-style escaped string where all instances of C-style escape sequences
    !>  are converted to the corresponding ASCII characters or left intact if they are not convertible.<br>
    !>
    !>  \details
    !>  See the documentation of [setAsciiFromEscaped](@ref pm_strASCII::setAsciiFromEscaped) for further details.<br>
    !>
    !>  \param[in]  str :   The input scalar `character` of kind \SKALL containing the C-style escaped string to be converted to ASCII.<br>
    !>
    !>  \return
    !>  `ascii`         :   The output `allocatable` scalar of the same type and kind as the input `str` containing
    !>                      the input string where all instances of escape seuquences with an ASCII representation are
    !>                      replaced with their corresponding ASCII character (, otherwise, left as is).<br>
    !>
    !>  \interface{getAsciiFromEscaped}
    !>  \code{.F90}
    !>
    !>      use pm_strASCII, only: getAsciiFromEscaped
    !>      character(:,kind(str)), allocatable :: ascii
    !>
    !>      ascii = getAsciiFromEscaped(str)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>  The impurity of the procedures is caused by the dependence on other conditionally `impure` procedure.<br>
    !>
    !>  \see
    !>  [setAsciiFromEscaped](@ref pm_strASCII::setAsciiFromEscaped)<br>
    !>  [isFailedList](@ref pm_sysPath::isFailedList)<br>
    !>
    !>  \example{getAsciiFromEscaped}
    !>  \include{lineno} example/pm_strASCII/getAsciiFromEscaped/main.F90
    !>  \compilef{getAsciiFromEscaped}
    !>  \output{getAsciiFromEscaped}
    !>  \include{lineno} example/pm_strASCII/getAsciiFromEscaped/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:02 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getAsciiFromEscaped

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getAsciiFromEscapedNew_SK5(str) result(ascii)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getAsciiFromEscapedNew_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        character(:,SKG)            , allocatable               :: ascii
    end function
#endif

#if SK4_ENABLED
    PURE module function getAsciiFromEscapedNew_SK4(str) result(ascii)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getAsciiFromEscapedNew_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        character(:,SKG)            , allocatable               :: ascii
    end function
#endif

#if SK3_ENABLED
    PURE module function getAsciiFromEscapedNew_SK3(str) result(ascii)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getAsciiFromEscapedNew_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        character(:,SKG)            , allocatable               :: ascii
    end function
#endif

#if SK2_ENABLED
    PURE module function getAsciiFromEscapedNew_SK2(str) result(ascii)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getAsciiFromEscapedNew_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        character(:,SKG)            , allocatable               :: ascii
    end function
#endif

#if SK1_ENABLED
    PURE module function getAsciiFromEscapedNew_SK1(str) result(ascii)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getAsciiFromEscapedNew_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        character(:,SKG)            , allocatable               :: ascii
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the input C-style escaped string where all instances of C-style escape sequences
    !>  are converted to the corresponding ASCII characters or left intact if they are not convertible.<br>
    !>
    !>  \details
    !>  Escape sequences are used in the programming languages C and C++ whose conventions are also followed by many other languages.<br>
    !>  An **escape sequence** is a sequence of characters that does not represent itself when used inside a character or string literal.<br>
    !>  Instead, it is translated into another character or a sequence of characters that may be difficult or impossible to represent directly.<br>
    !>  In C, all escape sequences consist of two or more characters, the first of which is the backslash, `\` (called the **escape character**).<br>
    !>  The remaining characters determine the interpretation of the escape sequence. For example, `\n` is an escape sequence that denotes a newline character.<br>
    !>  The following escape sequences are defined in standard C which commonly (but not always) represent an ASCII (frequently nongraphical) character:<br>
    !>
    !>  Escape sequence         |   ASCII Octal |   ASCII Decimal   |   ASCII Hex   |   ASCII Character Representation
    !>  ------------------------|---------------|-------------------|---------------|------------------------------------
    !>  \f$\ms{\a}\f$           |   07          |   07              |   07          |   Alert (Beep, Bell) (added in C89)
    !>  \f$\ms{\b}\f$           |   10          |   08              |   08          |   Backspace
    !>  \f$\ms{\f}\f$           |   14          |   12              |   0C          |   Formfeed Page Break
    !>  \f$\ms{\n}\f$           |   12          |   10              |   0A          |   Newline (Line Feed)
    !>  \f$\ms{\r}\f$           |   15          |   13              |   0D          |   Carriage Return
    !>  \f$\ms{\t}\f$           |   11          |   09              |   09          |   Horizontal Tab
    !>  \f$\ms{\v}\f$           |   13          |   11              |   0B          |   Vertical Tab
    !>  \f$\ms{\\}\f$           |   134         |   92              |   5C          |   Backslash
    !>  \f$\ms{\'}\f$           |   47          |   39              |   27          |   Apostrophe or single quotation mark
    !>  \f$\ms{\"}\f$           |   42          |   34              |   22          |   Double quotation mark
    !>  \f$\ms{\?}\f$           |   77          |   63              |   3F          |   Question mark (used to avoid trigraphs)
    !>  \f$\ms{\nnn}\f$         |   any         |   any             |   any         |   The byte whose numerical value is given by `nnn` interpreted as an **octal number**
    !>  \f$\ms{\xhh…}\f$        |   any         |   any             |   any         |   The byte whose numerical value is given by `hh…` interpreted as a **hexadecimal number**
    !>  \f$\ms{\uhhhh}\f$       |   none        |   none            |   none        |   Unicode code point below **10000 hexadecimal** (added in C99)
    !>  \f$\ms{\Uhhhhhhhh}\f$   |   none        |   none            |   none        |   Unicode code point where `h` is a **hexadecimal digit**
    !>
    !>  The following remarks are in order:
    !>
    !>  Octal Escape Sequences
    !>  ----------------------
    !>
    !>  An **octal escape sequence** consists of `\` followed by **one, two, or three** octal digits.<br>
    !>  <ol>
    !>      <li>    The octal escape sequence ends when it either contains three octal digits already, or the next character is not an octal digit.<br>
    !>              For example,
    !>              <ol>
    !>                  <li>    `\11` is a single octal escape sequence denoting a byte with numerical value `9` (`11` in octal), rather than the escape sequence `\1` followed by the digit `1`.<br>
    !>                  <li>    However, `\1111` is the octal escape sequence `\111` followed by the digit `1`.<br>
    !>              </ol>
    !>      <li>    In order to denote the byte with numerical value `1`, followed by the digit `1`, one could use `"\1""1"`, since C automatically concatenates adjacent string literals.<br>
    !>      <li>    Some three-digit octal escape sequences may be too large to fit in a single byte.<br>
    !>              This results in an implementation-defined value for the byte actually produced.<br>
    !>      <li>    The escape sequence `\0` is a commonly used octal escape sequence, which denotes the **null character**, with value zero.<br>
    !>      <li>    The procedures of this generic interface leave any octal sequence that is non-convertible to ASCII character intact (as is) in the output.<br>
    !>  </ol>
    !>
    !>  Hex Escape Sequences
    !>  --------------------
    !>
    !>  A **hex escape sequence** must have **at least one hex digit** following `\x`, with no upper bound.<br>
    !>  <ol>
    !>      <li>    It continues for as many hex digits as there are.<br>
    !>              For example, `\xABCDEFG` denotes the byte with the numerical value `ABCDEF16`, followed by the letter `G`, which is not a hex digit.<br>
    !>      <li>    However, if the resulting integer value is too large to fit in a single byte, the actual numerical value assigned is implementation-defined.<br>
    !>              Most platforms have `8`-bit character types, which limits a useful hex escape sequence to **two hex digits**.<br>
    !>              However, hex escape sequences longer than two hex digits might be useful inside a wide character or wide string.<br>
    !>      <li>    The procedures of this generic interface leave any hex sequence that is non-convertible to ASCII character intact (as is) in the output.<br>
    !>      <li>    The Hex alphabetical digits, if any are present, must be upper-case letters.<br>
    !>  </ol>
    !>
    !>  Universal Character Names
    !>  -------------------------
    !>
    !>  The C99 standard also supports escape sequences that denote Unicode code points in string literals.<br>
    !>  Such escape sequences are called **universal character names** and have the form `\uhhhh` or `\Uhhhhhhhh`, where `h` stands for a hex digit.<br>
    !>  <ol>
    !>      <li>    **Unlike** other escape sequences considered, a universal character name may expand into more than one code unit.<br>
    !>      <li>    The sequence `\uhhhh` denotes the code point `hhhh`, interpreted as a **hexadecimal number**.<br>
    !>      <li>    The sequence `\Uhhhhhhhh` denotes the code point `hhhhhhhh`, interpreted as a **hexadecimal number**.<br>
    !>      <li>    The code points located at `U+10000` or higher must be denoted with the `\U` syntax, whereas lower code points may use `\u` or `\U`.<br>
    !>      <li>    The code point is converted into a sequence of code units in the encoding of the destination type on the target system.<br>
    !>      <li>    The procedures of this generic interface leave any UCN that is non-convertible to ASCII character intact (as is) in the output.<br>
    !>      <li>    The Hex alphabetical digits in the UCN, if any are present, must be upper-case letters.<br>
    !>  </ol>
    !>
    !>  This functionality of this generic interface is highly useful for handling and maniulating C-style strings,
    !>  for example, in processing the output of runtime shell commands, or processing user-specified strings that
    !>  may or should contain nongraphical ASCII characters in a portable way.<br>
    !>
    !>  \param[inout]   str     :   The scalar `character` of kind \SKALL containing the C-style escaped string to be converted to ASCII.<br>
    !>                              -#  If the output argument `ascii` is present, then `str` has `intent(in)` attribute.<br>
    !>                              -#  If the output argument `ascii` is missing, then `str` has `intent(inout)` attribute.<br>
    !>                                  On output, the contents of `str(1:endloc)` will be overwritten with the ASCII equivalent.<br>
    !>  \param[out]     ascii   :   The output scalar of the same type and kind as the input `str` of length type parameter equal
    !>                              to or larger than the length of `str`, containing the input string where all instances of escape
    !>                              seuquences with an ASCII representation are replaced with their corresponding ASCII character.<br>
    !>                              (**optional**. If missing, the output ASCII will be written to `str`.)
    !>  \param[out]     endloc  :   The output scalar of type `integer` of default kind \IK, containing the position of
    !>                              the last character of the resulting ASCII-converted string (in either `str` or `ascii`).<br>
    !>
    !>  \interface{setAsciiFromEscaped}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: IK
    !>      use pm_strASCII, only: setAsciiFromEscaped
    !>      character(len(str,IK), kind(str)) :: ascii
    !>      integer(IK) :: endloc
    !>
    !>      call setAsciiFromEscaped(str, endloc)
    !>      call setAsciiFromEscaped(str, ascii, endloc)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>  The impurity of the procedures is caused by the dependence on other conditionally `impure` procedure.<br>
    !>
    !>  \note
    !>  Returning the output inside the input `str` should generally lead to faster runtime performance.<br>
    !>
    !>  \see
    !>  [getAsciiFromEscaped](@ref pm_strASCII::getAsciiFromEscaped)<br>
    !>  [isFailedList](@ref pm_sysPath::isFailedList)<br>
    !>
    !>  \example{setAsciiFromEscaped}
    !>  \include{lineno} example/pm_strASCII/setAsciiFromEscaped/main.F90
    !>  \compilef{setAsciiFromEscaped}
    !>  \output{setAsciiFromEscaped}
    !>  \include{lineno} example/pm_strASCII/setAsciiFromEscaped/main.out.F90
    !>
    !>  \test
    !>  [test_pm_strASCII](@ref test_pm_strASCII)
    !>
    !>  \todo
    !>  \phigh
    !>  A performance benchmarking of the different interfaces of this generic interface should be added in the future.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:02 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setAsciiFromEscaped

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setAsciiFromEscapedRep_SK5(str, endloc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAsciiFromEscapedRep_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout)             :: str
        integer(IK)                 , intent(out)               :: endloc
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setAsciiFromEscapedRep_SK4(str, endloc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAsciiFromEscapedRep_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout)             :: str
        integer(IK)                 , intent(out)               :: endloc
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setAsciiFromEscapedRep_SK3(str, endloc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAsciiFromEscapedRep_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout)             :: str
        integer(IK)                 , intent(out)               :: endloc
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setAsciiFromEscapedRep_SK2(str, endloc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAsciiFromEscapedRep_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout)             :: str
        integer(IK)                 , intent(out)               :: endloc
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setAsciiFromEscapedRep_SK1(str, endloc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAsciiFromEscapedRep_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout)             :: str
        integer(IK)                 , intent(out)               :: endloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setAsciiFromEscapedNew_SK5(str, ascii, endloc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAsciiFromEscapedNew_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                :: str
        character(*,SKG)            , intent(out)               :: ascii
        integer(IK)                 , intent(out)               :: endloc
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setAsciiFromEscapedNew_SK4(str, ascii, endloc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAsciiFromEscapedNew_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                :: str
        character(*,SKG)            , intent(out)               :: ascii
        integer(IK)                 , intent(out)               :: endloc
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setAsciiFromEscapedNew_SK3(str, ascii, endloc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAsciiFromEscapedNew_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                :: str
        character(*,SKG)            , intent(out)               :: ascii
        integer(IK)                 , intent(out)               :: endloc
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setAsciiFromEscapedNew_SK2(str, ascii, endloc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAsciiFromEscapedNew_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                :: str
        character(*,SKG)            , intent(out)               :: ascii
        integer(IK)                 , intent(out)               :: endloc
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setAsciiFromEscapedNew_SK1(str, ascii, endloc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setAsciiFromEscapedNew_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                :: str
        character(*,SKG)            , intent(out)               :: ascii
        integer(IK)                 , intent(out)               :: endloc
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \cond excluded
    impure elemental subroutine setStrAlphaRand(strAlphaRand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setStrAlphaRand
#endif
        use pm_kind, only: SK, IK
        character(*, SK), intent(inout) :: strAlphaRand
        integer(IK)     , parameter     :: BOUND_LOWER_CASE_L = ichar("a", kind = IK)
        integer(IK)     , parameter     :: BOUND_UPPER_CASE_L = ichar("z", kind = IK)
        integer(IK)     , parameter     :: BOUND_LOWER_CASE_U = ichar("A", kind = IK)
        integer(IK)     , parameter     :: BOUND_UPPER_CASE_U = ichar("Z", kind = IK)
        real                            :: randUnif, randUnifUL
        integer(IK)                     :: i
        do i = 1, len(strAlphaRand, IK)
            call random_number(randUnifUL)
            call random_number(randUnif)
            if (randUnifUL < 0.5) then
                strAlphaRand(i:i) = char( BOUND_LOWER_CASE_L + floor(randUnif * (BOUND_UPPER_CASE_L - BOUND_LOWER_CASE_L + 1_IK)) )
            else
                strAlphaRand(i:i) = char( BOUND_LOWER_CASE_U + floor(randUnif * (BOUND_UPPER_CASE_U - BOUND_LOWER_CASE_U + 1_IK)) )
            end if
        end do
    end subroutine setStrAlphaRand
    !>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_strASCII ! LCOV_EXCL_LINE