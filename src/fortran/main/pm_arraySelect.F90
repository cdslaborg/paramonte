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
!>  This module contains procedures and generic interfaces for selecting the \f$k\f$th smallest element in unsorted arrays of various types.
!>
!>  \benchmarks
!>
!>  \benchmark{getSelected_vs_setSelected, The runtime performance of [getSelected](@ref pm_arraySelect::getSelected) vs. [setSelected](@ref pm_arraySelect::setSelected)}
!>  \include{lineno} benchmark/pm_arraySelect/getSelected_vs_setSelected/main.F90
!>  \compilefb{getSelected_vs_setSelected}
!>  \postprocb{getSelected_vs_setSelected}
!>  \include{lineno} benchmark/pm_arraySelect/getSelected_vs_setSelected/main.py
!>  \visb{getSelected_vs_setSelected}
!>  \image html benchmark/pm_arraySelect/getSelected_vs_setSelected/benchmark.getSelected_vs_setSelected.runtime.png width=1000
!>  \image html benchmark/pm_arraySelect/getSelected_vs_setSelected/benchmark.getSelected_vs_setSelected.runtime.ratio.png width=1000
!>  \moralb{getSelected_vs_setSelected}
!>      -#  The procedures under the generic interface [getSelected](@ref pm_arraySelect::getSelected) are functions while
!>          the procedures under the generic interface [setSelected](@ref pm_arraySelect::setSelected) are subroutines.<br>
!>          From the benchmark results, it appears that the functional interface performs slightly less efficiently than the subroutine interface.<br>
!>          This is entirely due to the fact that the functional interface makes a full copy of the input array to keep the input array intact.<br>
!>      -#  Note that in this benchmark `rank = 1` was used to minimize the time spent on sorting, such that only the array copy time becomes the dominant factor.<br>
!>          Nevertheless, the benchmark results point to the relatively minor effect of array copying on the runtime performance of functional vs. subroutine procedures.<br>
!>      -#  There is, however, a \f$30\%-300\%\f$ performance loss with the use of the functional interface when `size(array) < 10`.<br>
!>
!>  \test
!>  [test_pm_arraySelect](@ref test_pm_arraySelect)
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

module pm_arraySelect

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_arraySelect"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    interface setSelectedIndex
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    PURE module subroutine setSelectedIndex_D0_SK(index, array, rank, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndex_D0_SK
!#endif
!        use pm_kind, only: LK, IK, SK
!        integer(IK)     , intent(out)                   :: index
!        character(*, SK), intent(in)                    :: array
!        integer(IK)     , intent(in)                    :: rank
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    PURE module subroutine setSelectedIndex_D1_SK(index, array, rank, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndex_D1_SK
!#endif
!        use pm_kind, only: LK, IK, SK
!        integer(IK)     , intent(out)                   :: index
!        character(*, SK), intent(in)    , contiguous    :: array(:)
!        integer(IK)     , intent(in)                    :: rank
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK4_ENABLED
!
!    PURE module subroutine setSelectedIndex_D1_IK4(index, array, rank, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndex_D1_IK4
!#endif
!        use pm_kind, only: LK, IK, IKG => IK4
!        integer(IK)     , intent(out)                   :: index
!        integer(IKG)    , intent(in)    , contiguous    :: array(:)
!        integer(IK)     , intent(in)                    :: rank
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK3_ENABLED
!
!    PURE module subroutine setSelectedIndex_D1_IK3(index, array, rank, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndex_D1_IK3
!#endif
!        use pm_kind, only: LK, IK, IKG => IK3
!        integer(IK)     , intent(out)                   :: index
!        integer(IKG)    , intent(in)    , contiguous    :: array(:)
!        integer(IK)     , intent(in)                    :: rank
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK2_ENABLED
!
!    PURE module subroutine setSelectedIndex_D1_IK2(index, array, rank, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndex_D1_IK2
!#endif
!        use pm_kind, only: LK, IK, IKG => IK2
!        integer(IK)     , intent(out)                   :: index
!        integer(IKG)    , intent(in)    , contiguous    :: array(:)
!        integer(IK)     , intent(in)                    :: rank
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK1_ENABLED
!
!    PURE module subroutine setSelectedIndex_D1_IK1(index, array, rank, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndex_D1_IK1
!#endif
!        use pm_kind, only: LK, IK, IKG => IK1
!        integer(IK)     , intent(out)                   :: index
!        integer(IKG)    , intent(in)    , contiguous    :: array(:)
!        integer(IK)     , intent(in)                    :: rank
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK3_ENABLED
!
!    PURE module subroutine setSelectedIndex_D1_RK3(index, array, rank, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndex_D1_RK3
!#endif
!        use pm_kind, only: LK, IK, RK => RK3
!        integer(IK)     , intent(out)                   :: index
!        real(RK)        , intent(in)    , contiguous    :: array(:)
!        integer(IK)     , intent(in)                    :: rank
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK2_ENABLED
!
!    PURE module subroutine setSelectedIndex_D1_RK2(index, array, rank, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndex_D1_RK2
!#endif
!        use pm_kind, only: LK, IK, RK => RK2
!        integer(IK)     , intent(out)                   :: index
!        real(RK)        , intent(in)    , contiguous    :: array(:)
!        integer(IK)     , intent(in)                    :: rank
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK1_ENABLED
!
!    PURE module subroutine setSelectedIndex_D1_RK1(index, array, rank, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndex_D1_RK1
!#endif
!        use pm_kind, only: LK, IK, RK => RK1
!        integer(IK)     , intent(out)                   :: index
!        real(RK)        , intent(in)    , contiguous    :: array(:)
!        integer(IK)     , intent(in)                    :: rank
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    PURE module subroutine setSelectedIndex_D1_SSK(index, array, rank, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndex_D1_SSK
!#endif
!        use pm_kind, only: IK, LK
!        use pm_container, only: Container => css_pdt
!        integer(IK)     , intent(out)                   :: index
!        type(Container) , intent(in)                    :: array(:)
!        integer(IK)     , intent(in)                    :: rank
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    module subroutine setSelectedIndexCustom_D0_SK(index, array, rank, isSorted, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndexCustom_D0_SK
!#endif
!        use pm_kind, only: LK, IK, SK
!        integer(IK)     , intent(out)                   :: index
!        character(*, SK), intent(in)                    :: array
!        integer(IK)     , intent(in)                    :: rank
!        logical(LK)     , external                      :: isSorted
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    module subroutine setSelectedIndexCustom_D1_SK(index, array, rank, isSorted, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndexCustom_D1_SK
!#endif
!        use pm_kind, only: LK, IK, SK
!        integer(IK)     , intent(out)                   :: index
!        character(*, SK), intent(in)    , contiguous    :: array(:)
!        integer(IK)     , intent(in)                    :: rank
!        logical(LK)     , external                      :: isSorted
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK4_ENABLED
!
!    module subroutine setSelectedIndexCustom_D1_IK4(index, array, rank, isSorted, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndexCustom_D1_IK4
!#endif
!        use pm_kind, only: LK, IK, IKG => IK4
!        integer(IK)     , intent(out)                   :: index
!        integer(IKG)    , intent(in)    , contiguous    :: array(:)
!        integer(IK)     , intent(in)                    :: rank
!        logical(LK)     , external                      :: isSorted
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK3_ENABLED
!
!    module subroutine setSelectedIndexCustom_D1_IK3(index, array, rank, isSorted, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndexCustom_D1_IK3
!#endif
!        use pm_kind, only: LK, IK, IKG => IK3
!        integer(IK)     , intent(out)                   :: index
!        integer(IKG)    , intent(in)    , contiguous    :: array(:)
!        integer(IK)     , intent(in)                    :: rank
!        logical(LK)     , external                      :: isSorted
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK2_ENABLED
!
!    module subroutine setSelectedIndexCustom_D1_IK2(index, array, rank, isSorted, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndexCustom_D1_IK2
!#endif
!        use pm_kind, only: LK, IK, IKG => IK2
!        integer(IK)     , intent(out)                   :: index
!        integer(IKG)    , intent(in)    , contiguous    :: array(:)
!        integer(IK)     , intent(in)                    :: rank
!        logical(LK)     , external                      :: isSorted
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if IK1_ENABLED
!
!    module subroutine setSelectedIndexCustom_D1_IK1(index, array, rank, isSorted, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndexCustom_D1_IK1
!#endif
!        use pm_kind, only: LK, IK, IKG => IK1
!        integer(IK)     , intent(out)                   :: index
!        integer(IKG)    , intent(in)    , contiguous    :: array(:)
!        integer(IK)     , intent(in)                    :: rank
!        logical(LK)     , external                      :: isSorted
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK3_ENABLED
!
!    module subroutine setSelectedIndexCustom_D1_RK3(index, array, rank, isSorted, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndexCustom_D1_RK3
!#endif
!        use pm_kind, only: LK, IK, RK => RK3
!        integer(IK)     , intent(out)                   :: index
!        real(RK)        , intent(in)    , contiguous    :: array(:)
!        integer(IK)     , intent(in)                    :: rank
!        logical(LK)     , external                      :: isSorted
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK2_ENABLED
!
!    module subroutine setSelectedIndexCustom_D1_RK2(index, array, rank, isSorted, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndexCustom_D1_RK2
!#endif
!        use pm_kind, only: LK, IK, RK => RK2
!        integer(IK)     , intent(out)                   :: index
!        real(RK)        , intent(in)    , contiguous    :: array(:)
!        integer(IK)     , intent(in)                    :: rank
!        logical(LK)     , external                      :: isSorted
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK1_ENABLED
!
!    module subroutine setSelectedIndexCustom_D1_RK1(index, array, rank, isSorted, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndexCustom_D1_RK1
!#endif
!        use pm_kind, only: LK, IK, RK => RK1
!        integer(IK)     , intent(out)                   :: index
!        real(RK)        , intent(in)    , contiguous    :: array(:)
!        integer(IK)     , intent(in)                    :: rank
!        logical(LK)     , external                      :: isSorted
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    module subroutine setSelectedIndexCustom_D1_SSK(index, array, rank, isSorted, lb, ub)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedIndexCustom_D1_SSK
!#endif
!        use pm_kind, only: IK, LK
!        use pm_container, only: Container => css_pdt
!        integer(IK)     , intent(out)                   :: index
!        type(Container) , intent(in)                    :: array(:)
!        integer(IK)     , intent(in)                    :: rank
!        logical(LK)     , external                      :: isSorted
!        integer(IK)     , intent(in)    , optional      :: lb, ub
!    end subroutine
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface setSelectedIndex

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the `rank`th smallest value in the input array by first sorting its elements in ascending order
    !>  (optionally only between the specified indices `[lb, ub]`).<br>
    !>
    !>  \details
    !>  By default, the array will be sorted in ascending order,
    !>  unless the user-specific `isSorted()` external function is supplied to define a custom sorting criterion.
    !>
    !>  \param[in]      array       :   The input `contiguous` array of shape `(:)` of either <br>
    !>                                  <ul>
    !>                                      <li>    type [css_type](@ref pm_container::css_type) (scalar string container of default kind \SK) or<br>
    !>                                      <li>    type [css_pdt](@ref pm_container::css_pdt) (scalar string container of kind \SKALL) or<br>
    !>                                      <li>    type `character` of kind \SKALL, or <br>
    !>                                      <li>    type `logical` of kind \LKALL, or <br>
    !>                                      <li>    type `integer` of kind \IKALL, or <br>
    !>                                      <li>    type `complex` of kind \CKALL, or <br>
    !>                                      <li>    type `real` of kind \RKALL, <br>
    !>                                  </ul>
    !>                                  or,
    !>                                  <ul>
    !>                                      <li>    a scalar assumed-length `character` of kind \SKALL, <br>
    !>                                  </ul>
    !>                                  containing the potentially unsorted sequence of values that is to be sorted
    !>                                  and whose `rank`th smallest value is to be reported as the output `selection`.
    !>  \param[in]      rank        :   The input `integer` of default kind \IK, representing the index of the `rank`th smallest value in the input `array`
    !>                                  or the `rank`th ordered value according to the user-supplied external `logical` function `isSorted()`.
    !>  \param          isSorted    :   The `external` user-specified function that takes two input **scalar** arguments of the same type
    !>                                  and kind as the input `array`.<br>
    !>                                  It returns a scalar `logical` of default kind \LK that is `.true.` if the first
    !>                                  input scalar argument is sorted with respect to the second input argument according to the user-defined condition
    !>                                  within `isSorted`, otherwise, it is `.false.`.<br>
    !>                                  If `array` is a Fortran string (i.e., an assumed-length scalar `character`),
    !>                                  then both input arguments to `isSorted()` are `character(1,SK)` of default kind \SK.<br>
    !>                                  The following illustrates the generic interface of `isSorted()`,
    !>                                  \code{.F90}
    !>                                      function isSorted(a,b) result (sorted)
    !>                                          use pm_kind, only: IK, LK
    !>                                          TYPE(KIND)  , intent(in)    :: a, b
    !>                                          logical(LK)                 :: sorted
    !>                                      end function
    !>                                  \endcode
    !>                                  where `TYPE(KIND)` is the same as the type and kind of the input argument `array`, which can be one of the following.
    !>                                  \code{.F90}
    !>                                      use pm_container, only: css_pdt, css_type
    !>                                      character(*, SK), intent(in) :: a, b
    !>                                      character(1, SK), intent(in) :: a, b
    !>                                      type(css_type)  , intent(in) :: a, b
    !>                                      type(css_pdt)   , intent(in) :: a, b
    !>                                      integer(IK)     , intent(in) :: a, b
    !>                                      real(RK)        , intent(in) :: a, b
    !>                                  \endcode
    !>                                  where the specified kind type parameters (`SK`, `IK`, `LK`, `CK`, `RK`) can refer to any of the supported kinds by the processor.<br>
    !>                                  This user-defined equivalence check is extremely useful where a user-defined sorting criterion other than simple ascending order
    !>                                  is needed, for example, when the case-sensitivity of an input string or array of strings  is irrelevant or when sorting of
    !>                                  the absolute values matters excluding the signs of the numbers, or when descending order is desired.<br>
    !>                                  In such cases, user can define a custom sorting condition within the user-defined external function `isSorted` to achieve the goal.<br>
    !>                                  (**optional**, the default sorting condition is ascending order, that is `a < b`.)
    !>  \param[in]          lb      :   The `optional` input scalar `integer` of default kind \IK representing the lower bound of the segment of array below which
    !>                                  the array elements are assumed to be in ascending order, such that only the elements starting with and beyond `lb` require sorting.<br>
    !>                                  If the initial segment of the input array is ordered, specifying this lower bound will lead to better performance of the algorithm.<br>
    !>                                  (**optional**, default = `1_IK`)
    !>  \param[in]          ub      :   The `optional` input scalar `integer` of default kind \IK representing the upper bound of the segment of array above which
    !>                                  the array elements are assumed to be in ascending order, such that only the elements ending with and below `ub` require sorting.<br>
    !>                                  If the final segment of the input array is ordered, specifying this upper bound will lead to better performance of the algorithm.<br>
    !>                                  (**optional**, default = `size(array, kind = IK)`)
    !>
    !>  \return
    !>  `selection`                 :   The output scalar of the same type and kind as `array` containing the `rank`th smallest (or ordered) value in the input `array`.<br>
    !>
    !>  \interface{getSelected}
    !>  \code{.F90}
    !>
    !>      use pm_arraySelect, only: getSelected
    !>
    !>      selection = getSelected(array, rank, lb = lb, ub = ub)
    !>      selection = getSelected(array, rank, isSorted, lb = lb, ub = ub)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  If `lb` is present as an input argument, its value must be larger than `0` and less than `ub` or its default value.<br>
    !>  Furthermore, the specified `rank` must be larger than or equal to the specified value for `lb`.<br>
    !>  If `ub` is present as an input argument, its value must be less than or equal to the length of `array` and larger than `ub` or its default value.<br>
    !>  Furthermore, the specified `rank` must be smaller than or equal to the specified value for `ub`.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  The functions under this generic interface are `impure` when the `external` input argument `isSorted()` is `present`.
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [setSelected](@ref pm_arraySelect::setSelected)<br>
    !>  [setRankDense](@ref pm_arrayRank::setRankDense)<br>
    !>  [setRankOrdinal](@ref pm_arrayRank::setRankOrdinal)<br>
    !>  [setRankFractional](@ref pm_arrayRank::setRankFractional)<br>
    !>  [setRankStandard](@ref pm_arrayRank::setRankStandard)<br>
    !>  [setRankModified](@ref pm_arrayRank::setRankModified)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [getRemoved](@ref pm_arrayRemove::getRemoved)<br>
    !>  [setReversed](@ref pm_arrayReverse::setReversed)<br>
    !>  [getReplaced](@ref pm_arrayReplace::getReplaced)<br>
    !>  [setReplaced](@ref pm_arrayReplace::setReplaced)<br>
    !>  [setSplit](@ref pm_arraySplit::setSplit)<br>
    !>
    !>  \example{getSelected}
    !>  \include{lineno} example/pm_arraySelect/getSelected/main.F90
    !>  \compilef{getSelected}
    !>  \output{getSelected}
    !>  \include{lineno} example/pm_arraySelect/getSelected/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arraySelect](@ref test_pm_arraySelect)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to higher-rank input objects.<br>
    !>
    !>  \final{getSelected}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! DefCom

    interface getSelected

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getSelectedDefCom_D0_SK5(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)                                            :: selection
    end function
#endif

#if SK4_ENABLED
    PURE module function getSelectedDefCom_D0_SK4(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)                                            :: selection
    end function
#endif

#if SK3_ENABLED
    PURE module function getSelectedDefCom_D0_SK3(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)                                            :: selection
    end function
#endif

#if SK2_ENABLED
    PURE module function getSelectedDefCom_D0_SK2(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)                                            :: selection
    end function
#endif

#if SK1_ENABLED
    PURE module function getSelectedDefCom_D0_SK1(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)                                            :: selection
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getSelectedDefCom_D1_SK5(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG)                                :: selection
    end function
#endif

#if SK4_ENABLED
    PURE module function getSelectedDefCom_D1_SK4(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG)                                :: selection
    end function
#endif

#if SK3_ENABLED
    PURE module function getSelectedDefCom_D1_SK3(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG)                                :: selection
    end function
#endif

#if SK2_ENABLED
    PURE module function getSelectedDefCom_D1_SK2(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG)                                :: selection
    end function
#endif

#if SK1_ENABLED
    PURE module function getSelectedDefCom_D1_SK1(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG)                                :: selection
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getSelectedDefCom_D1_IK5(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                                                :: selection
    end function
#endif

#if IK4_ENABLED
    PURE module function getSelectedDefCom_D1_IK4(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                                                :: selection
    end function
#endif

#if IK3_ENABLED
    PURE module function getSelectedDefCom_D1_IK3(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                                                :: selection
    end function
#endif

#if IK2_ENABLED
    PURE module function getSelectedDefCom_D1_IK2(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                                                :: selection
    end function
#endif

#if IK1_ENABLED
    PURE module function getSelectedDefCom_D1_IK1(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                                                :: selection
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getSelectedDefCom_D1_LK5(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                                                :: selection
    end function
#endif

#if LK4_ENABLED
    PURE module function getSelectedDefCom_D1_LK4(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                                                :: selection
    end function
#endif

#if LK3_ENABLED
    PURE module function getSelectedDefCom_D1_LK3(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                                                :: selection
    end function
#endif

#if LK2_ENABLED
    PURE module function getSelectedDefCom_D1_LK2(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                                                :: selection
    end function
#endif

#if LK1_ENABLED
    PURE module function getSelectedDefCom_D1_LK1(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                                                :: selection
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getSelectedDefCom_D1_CK5(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                                                :: selection
    end function
#endif

#if CK4_ENABLED
    PURE module function getSelectedDefCom_D1_CK4(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                                                :: selection
    end function
#endif

#if CK3_ENABLED
    PURE module function getSelectedDefCom_D1_CK3(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                                                :: selection
    end function
#endif

#if CK2_ENABLED
    PURE module function getSelectedDefCom_D1_CK2(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                                                :: selection
    end function
#endif

#if CK1_ENABLED
    PURE module function getSelectedDefCom_D1_CK1(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                                                :: selection
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getSelectedDefCom_D1_RK5(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                                                   :: selection
    end function
#endif

#if RK4_ENABLED
    PURE module function getSelectedDefCom_D1_RK4(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                                                   :: selection
    end function
#endif

#if RK3_ENABLED
    PURE module function getSelectedDefCom_D1_RK3(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                                                   :: selection
    end function
#endif

#if RK2_ENABLED
    PURE module function getSelectedDefCom_D1_RK2(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                                                   :: selection
    end function
#endif

#if RK1_ENABLED
    PURE module function getSelectedDefCom_D1_RK1(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                                                   :: selection
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getSelectedDefCom_D1_PSSK5(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))                                          :: selection
    end function
#endif

#if SK4_ENABLED
    module function getSelectedDefCom_D1_PSSK4(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))                                          :: selection
    end function
#endif

#if SK3_ENABLED
    module function getSelectedDefCom_D1_PSSK3(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))                                          :: selection
    end function
#endif

#if SK2_ENABLED
    module function getSelectedDefCom_D1_PSSK2(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))                                          :: selection
    end function
#endif

#if SK1_ENABLED
    module function getSelectedDefCom_D1_PSSK1(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))                                          :: selection
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getSelectedDefCom_D1_BSSK(array, rank, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedDefCom_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_type)                                              :: selection
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! CusCom

    interface getSelected

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getSelectedCusCom_D0_SK5(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)                                            :: selection
    end function
#endif

#if SK4_ENABLED
    module function getSelectedCusCom_D0_SK4(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)                                            :: selection
    end function
#endif

#if SK3_ENABLED
    module function getSelectedCusCom_D0_SK3(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)                                            :: selection
    end function
#endif

#if SK2_ENABLED
    module function getSelectedCusCom_D0_SK2(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)                                            :: selection
    end function
#endif

#if SK1_ENABLED
    module function getSelectedCusCom_D0_SK1(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)                                            :: selection
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getSelectedCusCom_D1_SK5(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG)                                :: selection
    end function
#endif

#if SK4_ENABLED
    module function getSelectedCusCom_D1_SK4(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG)                                :: selection
    end function
#endif

#if SK3_ENABLED
    module function getSelectedCusCom_D1_SK3(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG)                                :: selection
    end function
#endif

#if SK2_ENABLED
    module function getSelectedCusCom_D1_SK2(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG)                                :: selection
    end function
#endif

#if SK1_ENABLED
    module function getSelectedCusCom_D1_SK1(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG)                                :: selection
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getSelectedCusCom_D1_IK5(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                                                :: selection
    end function
#endif

#if IK4_ENABLED
    module function getSelectedCusCom_D1_IK4(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                                                :: selection
    end function
#endif

#if IK3_ENABLED
    module function getSelectedCusCom_D1_IK3(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                                                :: selection
    end function
#endif

#if IK2_ENABLED
    module function getSelectedCusCom_D1_IK2(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                                                :: selection
    end function
#endif

#if IK1_ENABLED
    module function getSelectedCusCom_D1_IK1(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                                                :: selection
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getSelectedCusCom_D1_LK5(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                                                :: selection
    end function
#endif

#if LK4_ENABLED
    module function getSelectedCusCom_D1_LK4(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                                                :: selection
    end function
#endif

#if LK3_ENABLED
    module function getSelectedCusCom_D1_LK3(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                                                :: selection
    end function
#endif

#if LK2_ENABLED
    module function getSelectedCusCom_D1_LK2(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                                                :: selection
    end function
#endif

#if LK1_ENABLED
    module function getSelectedCusCom_D1_LK1(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                                                :: selection
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getSelectedCusCom_D1_CK5(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                                                :: selection
    end function
#endif

#if CK4_ENABLED
    module function getSelectedCusCom_D1_CK4(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                                                :: selection
    end function
#endif

#if CK3_ENABLED
    module function getSelectedCusCom_D1_CK3(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                                                :: selection
    end function
#endif

#if CK2_ENABLED
    module function getSelectedCusCom_D1_CK2(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                                                :: selection
    end function
#endif

#if CK1_ENABLED
    module function getSelectedCusCom_D1_CK1(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                                                :: selection
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getSelectedCusCom_D1_RK5(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                                                   :: selection
    end function
#endif

#if RK4_ENABLED
    module function getSelectedCusCom_D1_RK4(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                                                   :: selection
    end function
#endif

#if RK3_ENABLED
    module function getSelectedCusCom_D1_RK3(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                                                   :: selection
    end function
#endif

#if RK2_ENABLED
    module function getSelectedCusCom_D1_RK2(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                                                   :: selection
    end function
#endif

#if RK1_ENABLED
    module function getSelectedCusCom_D1_RK1(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                                                   :: selection
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getSelectedCusCom_D1_PSSK5(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))                                          :: selection
    end function
#endif

#if SK4_ENABLED
    module function getSelectedCusCom_D1_PSSK4(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))                                          :: selection
    end function
#endif

#if SK3_ENABLED
    module function getSelectedCusCom_D1_PSSK3(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))                                          :: selection
    end function
#endif

#if SK2_ENABLED
    module function getSelectedCusCom_D1_PSSK2(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))                                          :: selection
    end function
#endif

#if SK1_ENABLED
    module function getSelectedCusCom_D1_PSSK1(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))                                          :: selection
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getSelectedCusCom_D1_BSSK(array, rank, isSorted, lb, ub) result(selection)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSelectedCusCom_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_type)                                              :: selection
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the `rank`th smallest (or ordered) value in the input array by first sorting its elements in ascending order
    !>  (optionally only between the specified indices `[lb, ub]`).<br>
    !>
    !>  \details
    !>  By default, the array will be sorted in ascending order,
    !>  unless the user-specific `isSorted()` external function is supplied to define a custom sorting criterion.
    !>
    !>  \param[out]     selection   :   The output scalar of either of the same type and kind as the input `array` argument,
    !>                                  containing the `rank`th smallest value in the input `array`.<br>
    !>  \param[inout]   array       :   The input/output `contiguous` array of shape `(:)` of either
    !>                                  <ul>
    !>                                      <li>    type [css_type](@ref pm_container::css_type) (scalar string container of default kind \SK) or<br>
    !>                                      <li>    type [css_pdt](@ref pm_container::css_pdt) (scalar string container of kind \SKALL) or<br>
    !>                                      <li>    type `character` of kind \SKALL, or <br>
    !>                                      <li>    type `logical` of kind \LKALL, or <br>
    !>                                      <li>    type `integer` of kind \IKALL, or <br>
    !>                                      <li>    type `complex` of kind \CKALL, or <br>
    !>                                      <li>    type `real` of kind \RKALL, <br>
    !>                                  </ul>
    !>                                  or,
    !>                                  <ul>
    !>                                      <li>    a scalar assumed-length `character` of kind \SKALL, <br>
    !>                                  </ul>
    !>                                  containing the potentially unsorted sequence of values that is to be
    !>                                  sorted and whose `rank`th smallest value is to be reported as the output `selection`.<br>
    !>                                  On return, the first `rank`th elements of `array` will be sorted in ascending order or in
    !>                                  the user-specified order as determined by the user-supplied external `isSorted()` argument.<br>
    !>  \param[in]      rank        :   The input `integer` of default kind \IK, representing the index of the `rank`th smallest value in the input `array`
    !>                                  or the `rank`th ordered value according to the user-supplied external `logical` function `isSorted()`.<br>
    !>  \param          isSorted    :   The `external` user-specified function that takes two input **scalar** arguments of the same type and kind as the input `array`.<br>
    !>                                  It returns a scalar `logical` of default kind \LK that is `.true.` if the first input scalar argument is sorted with respect to
    !>                                  the second input argument according to the user-defined condition within `isSorted`, otherwise, it is `.false.`.<br>
    !>                                  If `array` is a Fortran string (i.e., an assumed-length scalar `character`), then both
    !>                                  input arguments to `isSorted()` are `character(1,SK)` of default kind \SK.<br>
    !>                                  The following illustrates the generic interface of `isSorted()`,
    !>                                  \code{.F90}
    !>                                      function isSorted(a,b) result (sorted)
    !>                                          use pm_kind, only: IK, LK
    !>                                          TYPE(KIND)  , intent(in)    :: a, b
    !>                                          logical(LK)                 :: sorted
    !>                                      end function
    !>                                  \endcode
    !>                                  where `TYPE(KIND)` is the same as the type and kind of the input argument `array`, which can be one of the following.
    !>                                  \code{.F90}
    !>                                      use pm_container, only: css_pdt, css_type
    !>                                      character(*, SK), intent(in) :: a, b
    !>                                      character(1, SK), intent(in) :: a, b
    !>                                      type(css_type)  , intent(in) :: a, b
    !>                                      type(css_pdt)   , intent(in) :: a, b
    !>                                      integer(IK)     , intent(in) :: a, b
    !>                                      real(RK)        , intent(in) :: a, b
    !>                                  \endcode
    !>                                  where the specified kind type parameters (`SK`, `IK`, `LK`, `CK`, `RK`) can refer to any of the supported kinds by the processor.<br>
    !>                                  This user-defined equivalence check is extremely useful where a user-defined sorting criterion other than simple ascending order
    !>                                  is needed, for example, when the case-sensitivity of an input string or array of strings  is irrelevant or when sorting of
    !>                                  the absolute values matters excluding the signs of the numbers, or when descending order is desired.<br>
    !>                                  In such cases, user can define a custom sorting condition within the user-defined external function `isSorted` to achieve the goal.<br>
    !>                                  (**optional**, the default sorting condition is ascending order, that is `a < b`.)
    !>  \param[in]          lb      :   The `optional` input scalar `integer` of default kind \IK representing the lower bound of the segment of array below which
    !>                                  the array elements are assumed to be in ascending order, such that only the elements starting with and beyond `lb` require sorting.<br>
    !>                                  If the initial segment of the input array is ordered, specifying this lower bound will lead to better performance of the algorithm.<br>
    !>                                  (**optional**, default = `1_IK`)
    !>  \param[in]          ub      :   The `optional` input scalar `integer` of default kind \IK representing the upper bound of the segment of array above which
    !>                                  the array elements are assumed to be in ascending order, such that only the elements ending with and below `ub`
    !>                                  require sorting.<br>
    !>                                  If the final segment of the input array is ordered, specifying this upper bound will lead to better performance of the algorithm.<br>
    !>                                  (**optional**, default = `size(array, kind = IK)`)
    !>
    !>  \interface{setSelected}
    !>  \code{.F90}
    !>
    !>      use pm_arraySelect, only: setSelected
    !>
    !>      call setSelected(selection, array, rank, lb = lb, ub = ub)
    !>      call setSelected(selection, array, rank, isSorted, lb = lb, ub = ub)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  If `lb` is present as an input argument, its value must be larger than `0` and less than `ub` or its default value.<br>
    !>  Furthermore, the specified `rank` must be larger than or equal to the specified value for `lb`.<br>
    !>  If `ub` is present as an input argument, its value must be less than or equal to the length of `array` and larger than `ub` or its default value.<br>
    !>  Furthermore, the specified `rank` must be smaller than or equal to the specified value for `ub`.<br>
    !>  \vericons
    !>
    !>  \warning
    !>  The functions under this generic interface are `impure` when the `external` input argument `isSorted()` is `present`.
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getSelected](@ref pm_arraySelect::getSelected)<br>
    !>  [setRankDense](@ref pm_arrayRank::setRankDense)<br>
    !>  [setRankOrdinal](@ref pm_arrayRank::setRankOrdinal)<br>
    !>  [setRankFractional](@ref pm_arrayRank::setRankFractional)<br>
    !>  [setRankStandard](@ref pm_arrayRank::setRankStandard)<br>
    !>  [setRankModified](@ref pm_arrayRank::setRankModified)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [getRemoved](@ref pm_arrayRemove::getRemoved)<br>
    !>  [setReversed](@ref pm_arrayReverse::setReversed)<br>
    !>  [getReplaced](@ref pm_arrayReplace::getReplaced)<br>
    !>  [setReplaced](@ref pm_arrayReplace::setReplaced)<br>
    !>  [setSplit](@ref pm_arraySplit::setSplit)<br>
    !>
    !>  \example{setSelected}
    !>  \include{lineno} example/pm_arraySelect/setSelected/main.F90
    !>  \compilef{setSelected}
    !>  \output{setSelected}
    !>  \include{lineno} example/pm_arraySelect/setSelected/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arraySelect](@ref test_pm_arraySelect)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to higher-rank input objects.<br>
    !>
    !>  \final{setSelected}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! DefCom

    interface setSelected

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSelectedDefCom_D0_SK5(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)            , intent(out)                   :: selection
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSelectedDefCom_D0_SK4(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)            , intent(out)                   :: selection
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSelectedDefCom_D0_SK3(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)            , intent(out)                   :: selection
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSelectedDefCom_D0_SK2(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)            , intent(out)                   :: selection
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSelectedDefCom_D0_SK1(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)            , intent(out)                   :: selection
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSelectedDefCom_D1_SK5(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG), intent(out)                   :: selection
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSelectedDefCom_D1_SK4(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG), intent(out)                   :: selection
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSelectedDefCom_D1_SK3(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG), intent(out)                   :: selection
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSelectedDefCom_D1_SK2(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG), intent(out)                   :: selection
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSelectedDefCom_D1_SK1(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG), intent(out)                   :: selection
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setSelectedDefCom_D1_IK5(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setSelectedDefCom_D1_IK4(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setSelectedDefCom_D1_IK3(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setSelectedDefCom_D1_IK2(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setSelectedDefCom_D1_IK1(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                , intent(out)                   :: selection
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setSelectedDefCom_D1_LK5(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setSelectedDefCom_D1_LK4(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setSelectedDefCom_D1_LK3(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setSelectedDefCom_D1_LK2(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setSelectedDefCom_D1_LK1(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                , intent(out)                   :: selection
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setSelectedDefCom_D1_CK5(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setSelectedDefCom_D1_CK4(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setSelectedDefCom_D1_CK3(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setSelectedDefCom_D1_CK2(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setSelectedDefCom_D1_CK1(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                , intent(out)                   :: selection
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setSelectedDefCom_D1_RK5(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                   , intent(out)                   :: selection
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setSelectedDefCom_D1_RK4(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                   , intent(out)                   :: selection
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setSelectedDefCom_D1_RK3(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                   , intent(out)                   :: selection
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setSelectedDefCom_D1_RK2(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                   , intent(out)                   :: selection
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setSelectedDefCom_D1_RK1(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                   , intent(out)                   :: selection
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setSelectedDefCom_D1_PSSK5(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))          , intent(out)                   :: selection
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSelectedDefCom_D1_PSSK4(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))          , intent(out)                   :: selection
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSelectedDefCom_D1_PSSK3(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))          , intent(out)                   :: selection
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSelectedDefCom_D1_PSSK2(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))          , intent(out)                   :: selection
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSelectedDefCom_D1_PSSK1(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))          , intent(out)                   :: selection
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSelectedDefCom_D1_BSSK(selection, array, rank, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedDefCom_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_type)              , intent(out)                   :: selection
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! CusCom

    interface setSelected

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSelectedCusCom_D0_SK5(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)            , intent(out)                   :: selection
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSelectedCusCom_D0_SK4(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)            , intent(out)                   :: selection
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSelectedCusCom_D0_SK3(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)            , intent(out)                   :: selection
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSelectedCusCom_D0_SK2(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)            , intent(out)                   :: selection
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSelectedCusCom_D0_SK1(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout)                 :: array
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(1,SKG)            , intent(out)                   :: selection
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSelectedCusCom_D1_SK5(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG), intent(out)                   :: selection
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSelectedCusCom_D1_SK4(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG), intent(out)                   :: selection
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSelectedCusCom_D1_SK3(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG), intent(out)                   :: selection
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSelectedCusCom_D1_SK2(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG), intent(out)                   :: selection
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSelectedCusCom_D1_SK1(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        character(len(array,IK),SKG), intent(out)                   :: selection
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSelectedCusCom_D1_IK5(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSelectedCusCom_D1_IK4(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSelectedCusCom_D1_IK3(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSelectedCusCom_D1_IK2(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSelectedCusCom_D1_IK1(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        integer(IKG)                , intent(out)                   :: selection
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSelectedCusCom_D1_LK5(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSelectedCusCom_D1_LK4(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSelectedCusCom_D1_LK3(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSelectedCusCom_D1_LK2(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSelectedCusCom_D1_LK1(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        logical(LKG)                , intent(out)                   :: selection
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSelectedCusCom_D1_CK5(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSelectedCusCom_D1_CK4(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSelectedCusCom_D1_CK3(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSelectedCusCom_D1_CK2(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                , intent(out)                   :: selection
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSelectedCusCom_D1_CK1(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        complex(CKG)                , intent(out)                   :: selection
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSelectedCusCom_D1_RK5(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                   , intent(out)                   :: selection
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSelectedCusCom_D1_RK4(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                   , intent(out)                   :: selection
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSelectedCusCom_D1_RK3(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                   , intent(out)                   :: selection
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSelectedCusCom_D1_RK2(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                   , intent(out)                   :: selection
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSelectedCusCom_D1_RK1(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        real(RKG)                   , intent(out)                   :: selection
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setSelectedCusCom_D1_PSSK5(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))          , intent(out)                   :: selection
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSelectedCusCom_D1_PSSK4(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))          , intent(out)                   :: selection
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSelectedCusCom_D1_PSSK3(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))          , intent(out)                   :: selection
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSelectedCusCom_D1_PSSK2(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))          , intent(out)                   :: selection
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSelectedCusCom_D1_PSSK1(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_pdt(SKG))          , intent(out)                   :: selection
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSelectedCusCom_D1_BSSK(selection, array, rank, isSorted, lb, ub)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSelectedCusCom_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(in)                    :: rank
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                 , intent(in)    , optional      :: lb, ub
        type(css_type)              , intent(out)                   :: selection
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arraySelect ! LCOV_EXCL_LINE