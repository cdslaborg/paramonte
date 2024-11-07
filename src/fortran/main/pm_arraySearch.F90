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
!>  This module contains procedures and generic interfaces for finding the specific array index whose element
!>  has the largest value smaller than the input value in arrays of various types.<br>
!>
!>  \details
!>  <b>The input search array must be in ascending order.</b>
!>
!>  \benchmarks
!>
!>  \benchmark{getBin-isLessArg_vs_default, The runtime performance of [getBin](@ref pm_arraySearch::getBin) with and without `isLess` external optional input argument.}
!>  \include{lineno} benchmark/pm_arraySearch/getBin-isLessArg_vs_default/main.F90
!>  \compilefb{getBin-isLessArg_vs_default}
!>  \postprocb{getBin-isLessArg_vs_default}
!>  \include{lineno} benchmark/pm_arraySearch/getBin-isLessArg_vs_default/main.py
!>  \visb{getBin-isLessArg_vs_default}
!>  \image html benchmark/pm_arraySearch/getBin-isLessArg_vs_default/benchmark.getBin-isLessArg_vs_default.runtime.png width=1000
!>  \image html benchmark/pm_arraySearch/getBin-isLessArg_vs_default/benchmark.getBin-isLessArg_vs_default.runtime.ratio.png width=1000
!>  \moralb{getBin-isLessArg_vs_default}
!>      -#  The procedures under the generic interface [getBin](@ref pm_arraySearch::getBin) take an optional `isLess()` external-function input argument
!>          which can perform custom user-defined comparisons between the input `value` and `array` elements.<br>
!>          By default, the comparison is the `<` operator.<br>
!>      -#  The presence of an external `isLess` argument without inlining does appear to deteriorate the runtime performance of [getBin](@ref pm_arraySearch::getBin).<br>
!>          However, the exact amount of this penalty appears to significantly depend on the computing platform and the ability of the compiler to inline the external function.<br>
!>
!>  \test
!>  [test_pm_arraySearch](@ref test_pm_arraySearch)
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX
module pm_arraySearch

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_arraySearch"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `bin`, the index of the element of the input **ascending-ordered** `array`,
    !>  such that `array(bin) <= value < array(bin+1)` holds for the input `value`.<br>
    !>
    !>  \details
    !>  The following conditions hold.<br>
    !>  <ul>
    !>      <li>    If `value < array(1)`, then `bin = 0_IK`.<br>
    !>      <li>    If `array(size(array)) <= value`, then `bin = size(array, kind = IK)`.<br>
    !>      <li>    The <b>less-than</b> `<` comparison operator can be also customized by a user-defined comparison
    !>              supplied as the external function `isLess` input argument for arrays that are not ascending-ordered.<br>
    !>              In such a case, the output `bin` satisfies both `.not.(value < array(bin))` and `value < array(bin+1)`
    !>              where the `<` binary operator is defined by the `isLess()` user-supplied function.<br>
    !>  </ul>
    !>
    !>  \param[in]  array       :   The input `contiguous` array of shape `(:)` of either <br>
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter, or <br>
    !>                                  <li>    type `integer` of kind \IKALL, or <br>
    !>                                  <li>    type `complex` of kind \CKALL, or <br>
    !>                                  <li>    type `real` of kind \RKALL,<br>
    !>                              </ul>
    !>                              or,
    !>                              <ul>
    !>                                  <li>    a scalar assumed-length `character` of kind \SKALL, <br>
    !>                              </ul>
    !>                              whose elements will have to be searched for the largest value smaller than the input `value`.<br>
    !>                              The input `array` must be sorted in **ascending-order** unless an appropriate `isLess` external
    !>                              comparison function is supplied according to which the `array` effectively behaves as if it is ascending-ordered.<br>
    !>                              If the array is of type `complex`, then **only its real component will be compared** unless the comparison
    !>                              is defined by the `external` user-specified input comparison function `isLess`.
    !>  \param[in]  value       :   The input scalar of the same type and kind as the input array whose value is to be searched in `array`.
    !>  \param      isLess      :   The `external` user-specified function that takes two input scalar arguments of the same type and kind as the input `array`.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if the two input arguments meet the user-defined comparison criterion.<br>
    !>                              Otherwise, it is `.false.`. If `array` is a scalar `character`, the two input arguments to `isLess()` will be scalar `character` values of
    !>                              the same length as the input `value` (which could be larger than 1).<br>
    !>                              The following illustrates the generic interface of `isLess()`,
    !>                              \code{.F90}
    !>                                  function isLess(value, segment) result(isLess)
    !>                                      use pm_kind, only: LK
    !>                                      TYPE(KIND)  , intent(in)    :: value, segment
    !>                                      logical(LK)                 :: isLess
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `array`, which can be one of the following,
    !>                              \code{.F90}
    !>                                  integer(IK)             , intent(in) :: value, segment
    !>                                  complex(CK)             , intent(in) :: value, segment
    !>                                  real(RK)                , intent(in) :: value, segment
    !>                                  character(len(value),SK), intent(in) :: value, segment
    !>                              \endcode
    !>                              This `external` function is extremely useful where a user-defined comparison check other than `<` is desired,
    !>                              for example, when the array segments should match the input `value` only within a given threshold or,
    !>                              when the case-sensitivity in character comparisons do not matter, or when the input `array` can be considered
    !>                              as an ascending-ordered sequence under certain conditions, for example, when only the magnitudes of the array
    !>                              elements are considered or when the array elements are multiplied by `-1` or inverted
    !>                              (that is, when the array is in descending-order). See below for example use cases.<br>
    !>                              (**optional**, the default comparison operator is `<`.)
    !>
    !>  \return
    !>  `bin`                   :   The output scalar `integer` of default kind \IK representing the index of the element of `array`
    !>                              for which `array(bin) <= value < array(bin+1)` holds or if `isLess()` is specified, the following
    !>                              conditions hold,<br>
    !>                              \code{.F90}
    !>                                  isLess(value, array(bin+1)) .and. .not. isLess(value, array(bin)) ! evaluates to .true._LK
    !>                              \endcode
    !>
    !>  \interface{getBin}
    !>  \code{.F90}
    !>
    !>      use pm_arraySearch, only: getBin
    !>
    !>      bin = getBin(array, value)
    !>      bin = getBin(array, value, isLess)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The procedures under this generic interface are `impure` when the user-specified `external` procedure `isLess` is specified as input argument.<br>
    !>
    !>  \warning
    !>  Note that in Fortran, trailing blanks are ignored in character comparison, that is, `"Fortran" == "Fortran "` yields `.true.`.<br>
    !>
    !>  \warning
    !>  The input `array` must be a **non-empty** sequence of values that is sorted according to the criterion specified by the `external` function `isLess()`
    !>  or if it is missing, then the values are sorted in **ascending order**.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  Be mindful of scenario where there are duplicate values in the input `array`.<br>
    !>  In such cases, the returned `bin` is not necessarily the index of the first or the last occurrence of such value.<br>
    !>  See below for illustrative examples.<br>
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getLoc](@ref pm_arrayFind::getLoc)<br>
    !>  [setLoc](@ref pm_arrayFind::setLoc)<br>
    !>  [setReplaced](@ref pm_arrayReplace::setReplaced)<br>
    !>  [getReplaced](@ref pm_arrayReplace::getReplaced)<br>
    !>  [setInserted](@ref pm_arrayInsert::setInserted)<br>
    !>  [setSplit](@ref pm_arraySplit::setSplit)<br>
    !>
    !>  \example{getBin}
    !>  \include{lineno} example/pm_arraySearch/getBin/main.F90
    !>  \compilef{getBin}
    !>  \output{getBin}
    !>  \include{lineno} example/pm_arraySearch/getBin/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arraySearch](@ref test_pm_arraySearch)
    !>
    !>  \final{getBin}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getBin

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getBinDefCom_D0_D0_SK5(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if SK4_ENABLED
    PURE module function getBinDefCom_D0_D0_SK4(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if SK3_ENABLED
    PURE module function getBinDefCom_D0_D0_SK3(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if SK2_ENABLED
    PURE module function getBinDefCom_D0_D0_SK2(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if SK1_ENABLED
    PURE module function getBinDefCom_D0_D0_SK1(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getBinDefCom_D1_D0_SK5(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if SK4_ENABLED
    PURE module function getBinDefCom_D1_D0_SK4(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if SK3_ENABLED
    PURE module function getBinDefCom_D1_D0_SK3(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if SK2_ENABLED
    PURE module function getBinDefCom_D1_D0_SK2(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if SK1_ENABLED
    PURE module function getBinDefCom_D1_D0_SK1(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getBinDefCom_D1_D0_IK5(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if IK4_ENABLED
    PURE module function getBinDefCom_D1_D0_IK4(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if IK3_ENABLED
    PURE module function getBinDefCom_D1_D0_IK3(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if IK2_ENABLED
    PURE module function getBinDefCom_D1_D0_IK2(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if IK1_ENABLED
    PURE module function getBinDefCom_D1_D0_IK1(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getBinDefCom_D1_D0_CK5(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if CK4_ENABLED
    PURE module function getBinDefCom_D1_D0_CK4(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if CK3_ENABLED
    PURE module function getBinDefCom_D1_D0_CK3(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if CK2_ENABLED
    PURE module function getBinDefCom_D1_D0_CK2(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if CK1_ENABLED
    PURE module function getBinDefCom_D1_D0_CK1(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getBinDefCom_D1_D0_RK5(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if RK4_ENABLED
    PURE module function getBinDefCom_D1_D0_RK4(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if RK3_ENABLED
    PURE module function getBinDefCom_D1_D0_RK3(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if RK2_ENABLED
    PURE module function getBinDefCom_D1_D0_RK2(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

#if RK1_ENABLED
    PURE module function getBinDefCom_D1_D0_RK1(array, value) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinDefCom_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: value
        integer(IK)                                         :: bin
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getBinCusCom_D0_D0_SK5(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if SK4_ENABLED
    module function getBinCusCom_D0_D0_SK4(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if SK3_ENABLED
    module function getBinCusCom_D0_D0_SK3(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if SK2_ENABLED
    module function getBinCusCom_D0_D0_SK2(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if SK1_ENABLED
    module function getBinCusCom_D0_D0_SK1(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                :: array
        character(*,SKG)        , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getBinCusCom_D1_D0_SK5(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if SK4_ENABLED
    module function getBinCusCom_D1_D0_SK4(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if SK3_ENABLED
    module function getBinCusCom_D1_D0_SK3(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if SK2_ENABLED
    module function getBinCusCom_D1_D0_SK2(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if SK1_ENABLED
    module function getBinCusCom_D1_D0_SK1(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in), contiguous    :: array(:)
        character(*,SKG)        , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getBinCusCom_D1_D0_IK5(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if IK4_ENABLED
    module function getBinCusCom_D1_D0_IK4(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if IK3_ENABLED
    module function getBinCusCom_D1_D0_IK3(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if IK2_ENABLED
    module function getBinCusCom_D1_D0_IK2(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if IK1_ENABLED
    module function getBinCusCom_D1_D0_IK1(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in), contiguous    :: array(:)
        integer(IKG)            , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getBinCusCom_D1_D0_CK5(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if CK4_ENABLED
    module function getBinCusCom_D1_D0_CK4(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if CK3_ENABLED
    module function getBinCusCom_D1_D0_CK3(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if CK2_ENABLED
    module function getBinCusCom_D1_D0_CK2(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if CK1_ENABLED
    module function getBinCusCom_D1_D0_CK1(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in), contiguous    :: array(:)
        complex(CKG)            , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getBinCusCom_D1_D0_RK5(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if RK4_ENABLED
    module function getBinCusCom_D1_D0_RK4(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if RK3_ENABLED
    module function getBinCusCom_D1_D0_RK3(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if RK2_ENABLED
    module function getBinCusCom_D1_D0_RK2(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

#if RK1_ENABLED
    module function getBinCusCom_D1_D0_RK1(array, value, isLess) result(bin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBinCusCom_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in), contiguous    :: array(:)
        real(RKG)               , intent(in)                :: value
        procedure(logical(LK))                              :: isLess
        integer(IK)                                         :: bin
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arraySearch ! LCOV_EXCL_LINE