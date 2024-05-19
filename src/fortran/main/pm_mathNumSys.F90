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
!>  This module contains procedures and generic interfaces for converting numbers to different bases in different numeral systems.
!>
!>  \see
!>  [pm_arraySpace](@ref pm_arraySpace)<br>
!>
!>  \test
!>  [test_pm_mathNumSys](@ref test_pm_mathNumSys)
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Tuesday 02:49 AM, August 10, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathNumSys

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathNumSys"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate the decimal value corresponding to the input numeral string in the requested base.
    !>
    !>  \param[in]  numeral :   The input scalar, or array of the same shape as other array-like arguments,
    !>                          of type `character` of kind \SKALL, containing the input value in an arbitrary numeral system.<br>
    !>                          The input `numeral` must be non-empty and contain only ASCII **digital** and the first `max(base-10,0)` **upper-case** English characters.
    !>  \param[in]  base    :   The input scalar, or array of the same shape as other array-like arguments, of type `integer`
    !>                          of kind \IKALL containing the base of the numeral system of the input `numeral`.
    !>
    !>  \return
    !>  `decimal`           :   The output scalar, or array of the same shape as array-like input arguments, of the same type and kind as
    !>                          the input `base` containing the decimal value corresponding to the input `numeral` in the specified `base`.
    !>
    !>  \interface{getDecimal}
    !>  \code{.F90}
    !>
    !>      use pm_mathNumSys, only: getDecimal
    !>
    !>      decimal = getDecimal(numeral, base)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \warning
    !>  The input `base` must be a positive integer larger than `1`.<br>
    !>  The input `numeral` must be non-empty and contain only ASCII digital and the first `max(base-10,0)` **upper-case** English characters.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  While `base` can be arbitrarily large, the specific integer kind used for `base` may be insufficient for the returned `decimal` value.<br>
    !>  This can lead to integer overflows that may remain unnoticed in non-debug library builds if `numeral` or `base` is too large for the specified kind of `decimal`.<br>
    !>  Additionally, keep in mind that there are only `36` symbols for numeral system digits that can fully cover all numbers up to base `36`.<br>
    !>  It is the responsibility of the user to ensure the specified kind for the input `decimal` is sufficient for the output `decimal` value.
    !>
    !>  \see
    !>  [getNumeral](@ref pm_mathNumSys::getNumeral)<br>
    !>
    !>  \example{getDecimal}
    !>  \include{lineno} example/pm_mathNumSys/getDecimal/main.F90
    !>  \compilef{getDecimal}
    !>  \output{getDecimal}
    !>  \include{lineno} example/pm_mathNumSys/getDecimal/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathNumSys](@ref test_pm_mathNumSys)
    !>
    !>  \final{getDecimal}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 02:49 AM, August 10, 2021, Dallas, TX
    interface getDecimal

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED && IK5_ENABLED
    PURE elemental module function getDecimal_SK5_IK5(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK5_IK5
#endif
        use pm_kind, only: SKG => SK5, IKG => IK5
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK5_ENABLED && IK4_ENABLED
    PURE elemental module function getDecimal_SK5_IK4(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK5_IK4
#endif
        use pm_kind, only: SKG => SK5, IKG => IK4
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK5_ENABLED && IK3_ENABLED
    PURE elemental module function getDecimal_SK5_IK3(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK5_IK3
#endif
        use pm_kind, only: SKG => SK5, IKG => IK3
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK5_ENABLED && IK2_ENABLED
    PURE elemental module function getDecimal_SK5_IK2(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK5_IK2
#endif
        use pm_kind, only: SKG => SK5, IKG => IK2
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK5_ENABLED && IK1_ENABLED
    PURE elemental module function getDecimal_SK5_IK1(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK5_IK1
#endif
        use pm_kind, only: SKG => SK5, IKG => IK1
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK4_ENABLED && IK5_ENABLED
    PURE elemental module function getDecimal_SK4_IK5(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK4_IK5
#endif
        use pm_kind, only: SKG => SK4, IKG => IK5
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK4_ENABLED && IK4_ENABLED
    PURE elemental module function getDecimal_SK4_IK4(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK4_IK4
#endif
        use pm_kind, only: SKG => SK4, IKG => IK4
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK4_ENABLED && IK3_ENABLED
    PURE elemental module function getDecimal_SK4_IK3(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK4_IK3
#endif
        use pm_kind, only: SKG => SK4, IKG => IK3
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK4_ENABLED && IK2_ENABLED
    PURE elemental module function getDecimal_SK4_IK2(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK4_IK2
#endif
        use pm_kind, only: SKG => SK4, IKG => IK2
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK4_ENABLED && IK1_ENABLED
    PURE elemental module function getDecimal_SK4_IK1(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK4_IK1
#endif
        use pm_kind, only: SKG => SK4, IKG => IK1
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK3_ENABLED && IK5_ENABLED
    PURE elemental module function getDecimal_SK3_IK5(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK3_IK5
#endif
        use pm_kind, only: SKG => SK3, IKG => IK5
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK3_ENABLED && IK4_ENABLED
    PURE elemental module function getDecimal_SK3_IK4(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK3_IK4
#endif
        use pm_kind, only: SKG => SK3, IKG => IK4
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK3_ENABLED && IK3_ENABLED
    PURE elemental module function getDecimal_SK3_IK3(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK3_IK3
#endif
        use pm_kind, only: SKG => SK3, IKG => IK3
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK3_ENABLED && IK2_ENABLED
    PURE elemental module function getDecimal_SK3_IK2(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK3_IK2
#endif
        use pm_kind, only: SKG => SK3, IKG => IK2
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK3_ENABLED && IK1_ENABLED
    PURE elemental module function getDecimal_SK3_IK1(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK3_IK1
#endif
        use pm_kind, only: SKG => SK3, IKG => IK1
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK2_ENABLED && IK5_ENABLED
    PURE elemental module function getDecimal_SK2_IK5(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK2_IK5
#endif
        use pm_kind, only: SKG => SK2, IKG => IK5
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK2_ENABLED && IK4_ENABLED
    PURE elemental module function getDecimal_SK2_IK4(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK2_IK4
#endif
        use pm_kind, only: SKG => SK2, IKG => IK4
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK2_ENABLED && IK3_ENABLED
    PURE elemental module function getDecimal_SK2_IK3(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK2_IK3
#endif
        use pm_kind, only: SKG => SK2, IKG => IK3
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK2_ENABLED && IK2_ENABLED
    PURE elemental module function getDecimal_SK2_IK2(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK2_IK2
#endif
        use pm_kind, only: SKG => SK2, IKG => IK2
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK2_ENABLED && IK1_ENABLED
    PURE elemental module function getDecimal_SK2_IK1(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK2_IK1
#endif
        use pm_kind, only: SKG => SK2, IKG => IK1
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK1_ENABLED && IK5_ENABLED
    PURE elemental module function getDecimal_SK1_IK5(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK1_IK5
#endif
        use pm_kind, only: SKG => SK1, IKG => IK5
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK1_ENABLED && IK4_ENABLED
    PURE elemental module function getDecimal_SK1_IK4(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK1_IK4
#endif
        use pm_kind, only: SKG => SK1, IKG => IK4
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK1_ENABLED && IK3_ENABLED
    PURE elemental module function getDecimal_SK1_IK3(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK1_IK3
#endif
        use pm_kind, only: SKG => SK1, IKG => IK3
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK1_ENABLED && IK2_ENABLED
    PURE elemental module function getDecimal_SK1_IK2(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK1_IK2
#endif
        use pm_kind, only: SKG => SK1, IKG => IK2
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

#if SK1_ENABLED && IK1_ENABLED
    PURE elemental module function getDecimal_SK1_IK1(numeral, base) result(decimal)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDecimal_SK1_IK1
#endif
        use pm_kind, only: SKG => SK1, IKG => IK1
        character(*,SKG)        , intent(in)            :: numeral
        integer(IKG)            , intent(in)            :: base
        integer(IKG)                                    :: decimal
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate the `numeral` in the specified `base` corresponding to the input `decimal` number.
    !>
    !>  \param[in]  decimal :   The input scalar of type `integer` of kind \IKALL,
    !>                          containing the input non-negative decimal number to be converted to the requested `base`.<br>
    !>                          Note that `decimal` argument must be passed as `value`. This is automatically done by the compiler for a Fortran user.
    !>  \param[in]  base    :   The input scalar of type `integer` of the same kind as the `integer` kind of `decimal`,
    !>                          containing the base of the numeral system of the input `numeral`.<br>
    !>                          The input `base` must be larger than `1` and less than `36`.
    !>
    !>  \return
    !>  `numeral`           :   The `allocatable` output scalar `character` of default kind \SK containing the
    !>                          numeral in the specified `base`, corresponding to the input` decimal` number,
    !>                          comprised of only ASCII **digital** character and, if `base > 10`, also
    !>                          possibly the first `max(base-10,0)` **upper-case** English characters.
    !>
    !>  \interface{getNumeral}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: SK
    !>      use pm_mathNumSys, only: getNumeral
    !>      character(:, SK), allocatable :: numeral
    !>
    !>      numeral = getNumeral(decimal, base)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The input `base` must be a positive integer larger than `1`.<br>
    !>  The input `value` must be non-empty and contain only ASCII digital and the first `max(base-10,0)` **upper-case** English characters.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  While `base` can be arbitrarily large, the specific integer kind used for `base` may be insufficient for the returned `decimal` value.<br>
    !>  This can lead to integer overflows that may remain unnoticed in non-debug library builds if `numeral` or `base` is too large for the specified kind of `decimal`.<br>
    !>  Additionally, keep in mind that there are only `36` symbols for numeral system digits that can fully cover all numbers up to base `36`.<br>
    !>  It is the responsibility of the user to ensure the specified kind for the input `decimal` is sufficient for the output `decimal` value.
    !>
    !>  \see
    !>  [getDecimal](@ref pm_mathNumSys::getDecimal)<br>
    !>
    !>  \example{getNumeral}
    !>  \include{lineno} example/pm_mathNumSys/getNumeral/main.F90
    !>  \compilef{getNumeral}
    !>  \output{getNumeral}
    !>  \include{lineno} example/pm_mathNumSys/getNumeral/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathNumSys](@ref test_pm_mathNumSys)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.8.0 20221119}
    !>  \desc
    !>  The \ifort{2021.8.0 20221119} cannot run the following function call in the examples of [getDecimal](@ref pm_mathNumSys::getDecimal),
    !>  \code{.F90}
    !>      getNumeral(huge(1_IK), base = 36_IK)
    !>  \endcode
    !>  yielding the following error,<br>
    !>  \code{.sh}
    !>      forrtl: severe (174): SIGSEGV, segmentation fault occurred
    !>      image              PC                Routine            Line        Source
    !>      libc.so.6          000014F31F6D5520  Unknown               Unknown  Unknown
    !>      libparamonte_fort  000014F324F06B7D  pm_mathnumsys_MP_         118  pm_mathNumSys@routines.inc.F90
    !>      main.exe           000000000040B21A  MAIN__                     39  main.F90
    !>      main.exe           000000000040833D  Unknown               Unknown  Unknown
    !>      libc.so.6          000014F31F6BCD90  Unknown               Unknown  Unknown
    !>      libc.so.6          000014F31F6BCE40  __libc_start_main     Unknown  Unknown
    !>      main.exe           0000000000408255  Unknown               Unknown  Unknown
    !>  \endcode
    !>  This segmentation fault error resembles the compiler bug observed in [getFormat](@ref pm_io::getFormat).<br>
    !>  It seems to be related to the `value` attribute of the input argument `decimal`.<br>
    !>  \remedy
    !>  For now, the `value` attribute of the input argument `decimal` of
    !>  [getNumeral](@ref pm_mathNumSys::getNumeral) is converted to `intent(in)` to resolve the compiler bug.<br>
    !>
    !>  \todo
    !>  \pmed
    !>  There is a dependency on the kind type parameter of the `integer` input arguments requiring `range(0_IKG) < 1024`.<br>
    !>  The algorithm must be modified to become kind-agnostic.<br>
    !>
    !>  \final{getNumeral}
    !>
    !>  \author
    !>  \AmirShahmoradi, Tuesday 03:59 AM, August 10, 2021, Dallas, TX
    interface getNumeral

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getNumeral_IK5(decimal, base) result(numeral)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNumeral_IK5
#endif
        use pm_kind, only: SKG => SK, IKG => IK5
        integer(IKG)            , intent(in)            :: decimal
        integer(IKG)            , intent(in)            :: base
        character(:,SKG)        , allocatable           :: numeral
    end function
#endif

#if IK4_ENABLED
    PURE module function getNumeral_IK4(decimal, base) result(numeral)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNumeral_IK4
#endif
        use pm_kind, only: SKG => SK, IKG => IK4
        integer(IKG)            , intent(in)            :: decimal
        integer(IKG)            , intent(in)            :: base
        character(:,SKG)        , allocatable           :: numeral
    end function
#endif

#if IK3_ENABLED
    PURE module function getNumeral_IK3(decimal, base) result(numeral)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNumeral_IK3
#endif
        use pm_kind, only: SKG => SK, IKG => IK3
        integer(IKG)            , intent(in)            :: decimal
        integer(IKG)            , intent(in)            :: base
        character(:,SKG)        , allocatable           :: numeral
    end function
#endif

#if IK2_ENABLED
    PURE module function getNumeral_IK2(decimal, base) result(numeral)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNumeral_IK2
#endif
        use pm_kind, only: SKG => SK, IKG => IK2
        integer(IKG)            , intent(in)            :: decimal
        integer(IKG)            , intent(in)            :: base
        character(:,SKG)        , allocatable           :: numeral
    end function
#endif

#if IK1_ENABLED
    PURE module function getNumeral_IK1(decimal, base) result(numeral)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNumeral_IK1
#endif
        use pm_kind, only: SKG => SK, IKG => IK1
        integer(IKG)            , intent(in)            :: decimal
        integer(IKG)            , intent(in)            :: base
        character(:,SKG)        , allocatable           :: numeral
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the number of digits in the input integer.
    !>
    !>  \details
    !>  The number of digits of an integer is independent of its sign.<br>
    !>  For example, both `100` and `-100` have three digits in base `10`.<br>
    !>
    !>  \param[in]  val :   The input scalar or array of arbitrary rank of type `integer` of kind \IKALL
    !>                      containing the integer whose number of digits is to be computed on return.<br>
    !>                      Note that `val` has the `value` attribute in the procedure interfaces.
    !>
    !>  \return
    !>  `count`         :   The output scalar or array of the same shape as the input `val` of type `integer` od default kind \IK, containing the number of digits of the input integer `val`.
    !>
    !>  \interface{getCountDigit}
    !>  \code{.F90}
    !>
    !>      use pm_mathNumSys, only: getCountDigit
    !>      use pm_kind, only: IK
    !>      integer(IK) :: count
    !>
    !>      count = getCountDigit(val)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getCountDigit](@ref pm_mathNumSys::getCountDigit)<br>
    !>
    !>  \example{getCountDigit}
    !>  \include{lineno} example/pm_mathNumSys/getCountDigit/main.F90
    !>  \compilef{getCountDigit}
    !>  \output{getCountDigit}
    !>  \include{lineno} example/pm_mathNumSys/getCountDigit/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathNumSys](@ref test_pm_mathNumSys)
    !>
    !>  \todo
    !>  \pmed
    !>  The performance of this algorithm can be improved by invoking the proper approach to computing the number of digits.
    !>
    !>  \final{getCountDigit}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface getCountDigit

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure elemental module function getCountDigit_IK5(val) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountDigit_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG), intent(in)    :: val
        integer(IK)                 :: count
    end function
#endif

#if IK4_ENABLED
    pure elemental module function getCountDigit_IK4(val) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountDigit_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG), intent(in)    :: val
        integer(IK)                 :: count
    end function
#endif

#if IK3_ENABLED
    pure elemental module function getCountDigit_IK3(val) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountDigit_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG), intent(in)    :: val
        integer(IK)                 :: count
    end function
#endif

#if IK2_ENABLED
    pure elemental module function getCountDigit_IK2(val) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountDigit_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG), intent(in)    :: val
        integer(IK)                 :: count
    end function
#endif

#if IK1_ENABLED
    pure elemental module function getCountDigit_IK1(val) result(count)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCountDigit_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG), intent(in)    :: val
        integer(IK)                 :: count
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathNumSys ! LCOV_EXCL_LINE