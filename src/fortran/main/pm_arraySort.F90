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
!>  This module contains procedures and generic interfaces for various sorting tasks.<br>
!>
!>  \details
!>  Specifically, the generic interfaces of this module fall into three categories:<br>
!>  <ol>
!>          <li>    Testing whether an array of intrinsic type and kind is in
!>                  <ol>
!>                      <li>    [ascending](@ref pm_arraySort::isAscending),
!>                      <li>    [descending](@ref pm_arraySort::isDescending),
!>                      <li>    [sorted (i.e., either ascending or descending)](@ref pm_arraySort::isSorted),
!>                      <li>    [all ascending (i.e., strictly ascending with no duplicates)](@ref pm_arraySort::isSorted),
!>                      <li>    [all descending (i.e., strictly descending with no duplicates)](@ref pm_arraySort::isSorted),
!>                      <li>    [user-specified](@ref pm_arraySort::isSorted),
!>                  </ol>
!>                  order.<br>
!>          <li>    [Sorting the indices](@ref pm_arraySort::setSorted) of an array of intrinsic type and kind in ascending order.<br>
!>          <li>    [Sorting the elements](@ref pm_arraySort::setSorted) of an array of intrinsic type and kind in ascending order.<br>
!>          <li>    [Generating sorted indices](@ref pm_arraySort::getSorted) of an array of intrinsic type and kind in ascending order.<br>
!>          <li>    [Generating sorted elements](@ref pm_arraySort::getSorted) of an array of intrinsic type and kind in ascending order.<br>
!>  </ol>
!>
!>  There are currently twelve different sorting algorithms implemented in
!>  this module in addition to index sorting and sort checking algorithms.<br>
!>
!>  \note
!>  <b>Recommended routines for sorting arrays:</b><br>
!>  The procedures under the generic interface [setSorted](@ref pm_arraySort::setSorted) with default
!<  method are among the fastest sorting algorithms, particularly for fully random input arrays.<br>
!>
!>  \benchmarks
!>
!>  \benchmark{sorting}
!>  The following is program to test the performance of the different sorting algorithms in this module.<br>
!>  \include{lineno} benchmark/pm_arraySort/sorting/main.F90
!>  \compilefb{sorting}
!>  \postprocb{sorting}
!>  \include{lineno} benchmark/pm_arraySort/sorting/main.py
!>  \visb{sorting}
!>  \image html benchmark/pm_arraySort/sorting/benchmark.sorting.random.png width=1000
!>  \image html benchmark/pm_arraySort/sorting/benchmark.sorting.sorted.png width=1000
!>  \image html benchmark/pm_arraySort/sorting/benchmark.sorting.random.sorted.ratio.png width=1000
!>
!>  \benchmark{sorting_vs_indexing, sorting vs. indexing}
!>  The following program generates a performance comparison of the [setSorted](@ref pm_arraySort::setSorted) algorithm for various sorting methods.
!>  \include{lineno} benchmark/pm_arraySort/sorting_vs_indexing/main.F90
!>  \compilefb{sorting_vs_indexing}
!>  \postprocb{sorting_vs_indexing}
!>  \include{lineno} benchmark/pm_arraySort/sorting_vs_indexing/main.py
!>  \visb{sorting_vs_indexing}
!>  \image html benchmark/pm_arraySort/sorting_vs_indexing/benchmark.sorting_vs_indexing.random.png width=1000
!>  \image html benchmark/pm_arraySort/sorting_vs_indexing/benchmark.sorting_vs_indexing.random.ratio.png width=1000
!>  \image html benchmark/pm_arraySort/sorting_vs_indexing/benchmark.sorting_vs_indexing.sorted.png width=1000
!>  \image html benchmark/pm_arraySort/sorting_vs_indexing/benchmark.sorting_vs_indexing.sorted.ratio.png width=1000
!>
!>  \todo
!>  \pvlow
!>  An equivalent functional versions of [setSorted](@ref pm_arraySort::setSorted) and
!>  [setSorted](@ref pm_arraySort::setSorted) could be added along with the relevant benchmarks.<br>
!>
!>  \remark
!>  The sorting routines of this module are inspired by (although substantially different from),
!>  +   the works of <a href="https://www.mjr19.org.uk/IT/sorts/" target="_blank">Dr. Michael Rutter</a> and,
!>  +   the [Numerical Recipes in Fortran](https://www.amazon.com/Numerical-Recipes-Fortran-Scientific-Computing/dp/052143064X/) by Press et al., 1992.
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arraySort

    use pm_kind, only: IK, LK, RK, SK

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_arraySort"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type :: isort_type; end type
    type, abstract :: sort_type; end type
    type, extends(sort_type) :: qsort_type; end type
    type, extends(qsort_type) :: qsorti_type; end type
    type, extends(sort_type) :: heap_type; end type
    type, extends(sort_type) :: insertion_type; end type

    type, extends(qsort_type) :: qsortr_type; end type
    type, extends(qsort_type) :: qsortrdp_type; end type
    type, extends(sort_type) :: bubble_type; end type
    type, extends(heap_type):: heapi_type; end type
    type, extends(heap_type):: heapr_type; end type
    type, extends(insertion_type) :: insertionl_type; end type
    type, extends(insertion_type) :: insertionb_type; end type
    type, extends(sort_type) :: merger_type; end type
    type, extends(sort_type) :: selection_type; end type
    type, extends(sort_type) :: shell_type; end type

    type(isort_type), parameter :: isort = isort_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: isort
#endif

    type(qsorti_type), parameter :: qsorti = qsorti_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: qsorti
#endif
    type(qsortr_type), parameter :: qsortr = qsortr_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: qsortr
#endif
    type(qsortrdp_type), parameter :: qsortrdp = qsortrdp_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: qsortrdp
#endif
    type(bubble_type), parameter :: bubble = bubble_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: bubble
#endif
    type(heapi_type), parameter :: heapi = heapi_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: heapi
#endif
    type(heapr_type), parameter :: heapr = heapr_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: heapr
#endif
    type(insertionl_type), parameter :: insertionl = insertionl_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: insertionl
#endif
    type(insertionb_type), parameter :: insertionb = insertionb_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: insertionb
#endif
    type(merger_type), parameter :: merger = merger_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: merger
#endif
    type(selection_type), parameter :: selection = selection_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: selection
#endif
    type(shell_type), parameter :: shell = shell_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: shell
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input array is sorted in **strictly** all ascending order (**without** equal elements),
    !>  otherwise, generate and return `.false.`.
    !>
    !>  \param[in]  array       :   The input `contiguous` array of rank `1` of either<br>
    !>                              <ol>
    !>                                  <li>    type [css_pdt](@ref pm_container::css_pdt) or,<br>
    !>                                  <li>    type [css_type](@ref pm_container::css_type) or,<br>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter or,<br>
    !>                                  <li>    type `integer` of kind \IKALL or,<br>
    !>                                  <li>    type `logical` of kind \LKALL or,<br>
    !>                                  <li>    type `complex` of kind \CKALL or,<br>
    !>                                  <li>    type `real` of kind \RKALL or,<br>
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    a *scalar* of type `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                              </ol>
    !>                              whose elements will be checked for a **strictly** all ascending order (**without** equal elements).
    !>
    !>  \return
    !>  `ascendingAll`          :   The output scalar `logical` of default kind \LK that is `.true.` if elements of the input array are
    !>                              in **strictly** all ascending order (**without** equal elements), Otherwise, it is `.false.`.
    !>
    !>  \interface{isAscendingAll}
    !>  \code{.F90}
    !>
    !>      use pm_arraySort, only: isAscendingAll
    !>
    !>      ascendingAll = isAscendingAll(array)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The output of this procedure is `.true.` when the input array has zero length.
    !>
    !>  \warnpure
    !>  The procedures under this generic interface are always `impure` when the input argument `isSorted()` is present.
    !>
    !>  \see
    !>  [isSorted](@ref pm_arraySort::isSorted)<br>
    !>  [isDescending](@ref pm_arraySort::isDescending)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>
    !>  \example{isAscendingAll}
    !>  \include{lineno} example/pm_arraySort/isAscendingAll/main.F90
    !>  \compilef{isAscendingAll}
    !>  \output{isAscendingAll}
    !>  \include{lineno} example/pm_arraySort/isAscendingAll/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arraySort](@ref test_pm_arraySort)
    !>
    !>  \todo
    !>  This interface can be extended to scalar containers of strings.
    !>
    !>  \final{isAscendingAll}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 3:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    interface isAscendingAll

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function isAscendingAllDefCom_D0_SK5(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if SK4_ENABLED
    PURE module function isAscendingAllDefCom_D0_SK4(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if SK3_ENABLED
    PURE module function isAscendingAllDefCom_D0_SK3(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if SK2_ENABLED
    PURE module function isAscendingAllDefCom_D0_SK2(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if SK1_ENABLED
    PURE module function isAscendingAllDefCom_D0_SK1(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: ascendingAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function isAscendingAllDefCom_D1_SK5(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if SK4_ENABLED
    PURE module function isAscendingAllDefCom_D1_SK4(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if SK3_ENABLED
    PURE module function isAscendingAllDefCom_D1_SK3(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if SK2_ENABLED
    PURE module function isAscendingAllDefCom_D1_SK2(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if SK1_ENABLED
    PURE module function isAscendingAllDefCom_D1_SK1(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function isAscendingAllDefCom_D1_IK5(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if IK4_ENABLED
    PURE module function isAscendingAllDefCom_D1_IK4(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if IK3_ENABLED
    PURE module function isAscendingAllDefCom_D1_IK3(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if IK2_ENABLED
    PURE module function isAscendingAllDefCom_D1_IK2(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if IK1_ENABLED
    PURE module function isAscendingAllDefCom_D1_IK1(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function isAscendingAllDefCom_D1_LK5(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if LK4_ENABLED
    PURE module function isAscendingAllDefCom_D1_LK4(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if LK3_ENABLED
    PURE module function isAscendingAllDefCom_D1_LK3(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if LK2_ENABLED
    PURE module function isAscendingAllDefCom_D1_LK2(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if LK1_ENABLED
    PURE module function isAscendingAllDefCom_D1_LK1(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function isAscendingAllDefCom_D1_RK5(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if RK4_ENABLED
    PURE module function isAscendingAllDefCom_D1_RK4(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if RK3_ENABLED
    PURE module function isAscendingAllDefCom_D1_RK3(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if RK2_ENABLED
    PURE module function isAscendingAllDefCom_D1_RK2(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if RK1_ENABLED
    PURE module function isAscendingAllDefCom_D1_RK1(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function isAscendingAllDefCom_D1_CK5(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if CK4_ENABLED
    PURE module function isAscendingAllDefCom_D1_CK4(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if CK3_ENABLED
    PURE module function isAscendingAllDefCom_D1_CK3(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if CK2_ENABLED
    PURE module function isAscendingAllDefCom_D1_CK2(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if CK1_ENABLED
    PURE module function isAscendingAllDefCom_D1_CK1(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module function isAscendingAllDefCom_D1_PSSK5(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if SK4_ENABLED
    PURE module function isAscendingAllDefCom_D1_PSSK4(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if SK3_ENABLED
    PURE module function isAscendingAllDefCom_D1_PSSK3(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if SK2_ENABLED
    PURE module function isAscendingAllDefCom_D1_PSSK2(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#if SK1_ENABLED
    PURE module function isAscendingAllDefCom_D1_PSSK1(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function isAscendingAllDefCom_D1_BSSK(array) result(ascendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingAllDefCom_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascendingAll
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input array is sorted in **strictly** all descending order (**without** equal elements),
    !>  otherwise, generate and return `.false.`.
    !>
    !>  \param[in]  array       :   The input `contiguous` array of rank `1` of either<br>
    !>                              <ol>
    !>                                  <li>    type [css_pdt](@ref pm_container::css_pdt) or,<br>
    !>                                  <li>    type [css_type](@ref pm_container::css_type) or,<br>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter or,<br>
    !>                                  <li>    type `integer` of kind \IKALL or,<br>
    !>                                  <li>    type `logical` of kind \LKALL or,<br>
    !>                                  <li>    type `complex` of kind \CKALL or,<br>
    !>                                  <li>    type `real` of kind \RKALL or,<br>
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    a *scalar* of type `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                              </ol>
    !>                              whose elements will be checked for a **strictly** all descending order (**without** equal elements).
    !>
    !>  \return
    !>  `descendingAll`          :   The output scalar `logical` of default kind \LK that is `.true.` if elements of the input array are
    !>                              in **strictly** all descending order (**without** equal elements), Otherwise, it is `.false.`.
    !>
    !>  \interface{isDescendingAll}
    !>  \code{.F90}
    !>
    !>      use pm_arraySort, only: isDescendingAll
    !>
    !>      descendingAll = isDescendingAll(array)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The output of this procedure is `.true.` when the input array has zero length.
    !>
    !>  \warnpure
    !>  The procedures under this generic interface are always `impure` when the input argument `isSorted()` is present.
    !>
    !>  \see
    !>  [isSorted](@ref pm_arraySort::isSorted)<br>
    !>  [isDescending](@ref pm_arraySort::isDescending)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>
    !>  \example{isDescendingAll}
    !>  \include{lineno} example/pm_arraySort/isDescendingAll/main.F90
    !>  \compilef{isDescendingAll}
    !>  \output{isDescendingAll}
    !>  \include{lineno} example/pm_arraySort/isDescendingAll/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arraySort](@ref test_pm_arraySort)
    !>
    !>  \todo
    !>  This interface can be extended to scalar containers of strings.
    !>
    !>  \final{isDescendingAll}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 3:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    interface isDescendingAll

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function isDescendingAllDefCom_D0_SK5(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: descendingAll
    end function
#endif

#if SK4_ENABLED
    PURE module function isDescendingAllDefCom_D0_SK4(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: descendingAll
    end function
#endif

#if SK3_ENABLED
    PURE module function isDescendingAllDefCom_D0_SK3(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: descendingAll
    end function
#endif

#if SK2_ENABLED
    PURE module function isDescendingAllDefCom_D0_SK2(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: descendingAll
    end function
#endif

#if SK1_ENABLED
    PURE module function isDescendingAllDefCom_D0_SK1(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: descendingAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function isDescendingAllDefCom_D1_SK5(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if SK4_ENABLED
    PURE module function isDescendingAllDefCom_D1_SK4(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if SK3_ENABLED
    PURE module function isDescendingAllDefCom_D1_SK3(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if SK2_ENABLED
    PURE module function isDescendingAllDefCom_D1_SK2(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if SK1_ENABLED
    PURE module function isDescendingAllDefCom_D1_SK1(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function isDescendingAllDefCom_D1_IK5(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if IK4_ENABLED
    PURE module function isDescendingAllDefCom_D1_IK4(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if IK3_ENABLED
    PURE module function isDescendingAllDefCom_D1_IK3(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if IK2_ENABLED
    PURE module function isDescendingAllDefCom_D1_IK2(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if IK1_ENABLED
    PURE module function isDescendingAllDefCom_D1_IK1(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function isDescendingAllDefCom_D1_LK5(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if LK4_ENABLED
    PURE module function isDescendingAllDefCom_D1_LK4(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if LK3_ENABLED
    PURE module function isDescendingAllDefCom_D1_LK3(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if LK2_ENABLED
    PURE module function isDescendingAllDefCom_D1_LK2(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if LK1_ENABLED
    PURE module function isDescendingAllDefCom_D1_LK1(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function isDescendingAllDefCom_D1_RK5(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if RK4_ENABLED
    PURE module function isDescendingAllDefCom_D1_RK4(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if RK3_ENABLED
    PURE module function isDescendingAllDefCom_D1_RK3(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if RK2_ENABLED
    PURE module function isDescendingAllDefCom_D1_RK2(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if RK1_ENABLED
    PURE module function isDescendingAllDefCom_D1_RK1(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function isDescendingAllDefCom_D1_CK5(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if CK4_ENABLED
    PURE module function isDescendingAllDefCom_D1_CK4(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if CK3_ENABLED
    PURE module function isDescendingAllDefCom_D1_CK3(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if CK2_ENABLED
    PURE module function isDescendingAllDefCom_D1_CK2(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if CK1_ENABLED
    PURE module function isDescendingAllDefCom_D1_CK1(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module function isDescendingAllDefCom_D1_PSSK5(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if SK4_ENABLED
    PURE module function isDescendingAllDefCom_D1_PSSK4(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if SK3_ENABLED
    PURE module function isDescendingAllDefCom_D1_PSSK3(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if SK2_ENABLED
    PURE module function isDescendingAllDefCom_D1_PSSK2(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#if SK1_ENABLED
    PURE module function isDescendingAllDefCom_D1_PSSK1(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function isDescendingAllDefCom_D1_BSSK(array) result(descendingAll)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingAllDefCom_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descendingAll
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input array is sorted in ascending order (with the possibility of elements being equal),
    !>  otherwise, generate and return `.false.`.
    !>
    !>  \param[in]  array       :   The input `contiguous` array of rank `1` of either<br>
    !>                              <ol>
    !>                                  <li>    type [css_pdt](@ref pm_container::css_pdt) or,<br>
    !>                                  <li>    type [css_type](@ref pm_container::css_type) or,<br>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter or,<br>
    !>                                  <li>    type `integer` of kind \IKALL or,<br>
    !>                                  <li>    type `logical` of kind \LKALL or,<br>
    !>                                  <li>    type `complex` of kind \CKALL or,<br>
    !>                                  <li>    type `real` of kind \RKALL or,<br>
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    a *scalar* of type `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                              </ol>
    !>                              whose elements will be checked for an all descending order (with the possibility of elements being equal).
    !>
    !>  \return
    !>  `ascending`             :   An output `logical` of default kind \LK that is `.true.` if all elements of the input array are
    !>                              in descending order (with the possibility of elements being equal) or all equal. Otherwise, it is `.false.`.
    !>
    !>  \interface{isAscending}
    !>  \code{.F90}
    !>
    !>      use pm_arraySort, only: isAscending
    !>
    !>      ascending = isAscending(array)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The output of this procedure is `.true.` when the input array has zero length.
    !>
    !>  \warnpure
    !>  The procedures under this generic interface are always `impure` when the input argument `isSorted()` is present.
    !>
    !>  \see
    !>  [isSorted](@ref pm_arraySort::isSorted)<br>
    !>  [isDescending](@ref pm_arraySort::isDescending)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>
    !>  \example{isAscending}
    !>  \include{lineno} example/pm_arraySort/isAscending/main.F90
    !>  \compilef{isAscending}
    !>  \output{isAscending}
    !>  \include{lineno} example/pm_arraySort/isAscending/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arraySort](@ref test_pm_arraySort)
    !>
    !>  \todo
    !>  This interface can be extended to scalar containers of strings.
    !>
    !>  \final{isAscending}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 3:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    interface isAscending

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function isAscendingDefCom_D0_SK5(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: ascending
    end function
#endif

#if SK4_ENABLED
    PURE module function isAscendingDefCom_D0_SK4(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: ascending
    end function
#endif

#if SK3_ENABLED
    PURE module function isAscendingDefCom_D0_SK3(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: ascending
    end function
#endif

#if SK2_ENABLED
    PURE module function isAscendingDefCom_D0_SK2(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: ascending
    end function
#endif

#if SK1_ENABLED
    PURE module function isAscendingDefCom_D0_SK1(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: ascending
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function isAscendingDefCom_D1_SK5(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if SK4_ENABLED
    PURE module function isAscendingDefCom_D1_SK4(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if SK3_ENABLED
    PURE module function isAscendingDefCom_D1_SK3(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if SK2_ENABLED
    PURE module function isAscendingDefCom_D1_SK2(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if SK1_ENABLED
    PURE module function isAscendingDefCom_D1_SK1(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function isAscendingDefCom_D1_IK5(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if IK4_ENABLED
    PURE module function isAscendingDefCom_D1_IK4(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if IK3_ENABLED
    PURE module function isAscendingDefCom_D1_IK3(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if IK2_ENABLED
    PURE module function isAscendingDefCom_D1_IK2(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if IK1_ENABLED
    PURE module function isAscendingDefCom_D1_IK1(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function isAscendingDefCom_D1_LK5(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if LK4_ENABLED
    PURE module function isAscendingDefCom_D1_LK4(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if LK3_ENABLED
    PURE module function isAscendingDefCom_D1_LK3(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if LK2_ENABLED
    PURE module function isAscendingDefCom_D1_LK2(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if LK1_ENABLED
    PURE module function isAscendingDefCom_D1_LK1(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function isAscendingDefCom_D1_RK5(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if RK4_ENABLED
    PURE module function isAscendingDefCom_D1_RK4(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if RK3_ENABLED
    PURE module function isAscendingDefCom_D1_RK3(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if RK2_ENABLED
    PURE module function isAscendingDefCom_D1_RK2(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if RK1_ENABLED
    PURE module function isAscendingDefCom_D1_RK1(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function isAscendingDefCom_D1_CK5(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if CK4_ENABLED
    PURE module function isAscendingDefCom_D1_CK4(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if CK3_ENABLED
    PURE module function isAscendingDefCom_D1_CK3(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if CK2_ENABLED
    PURE module function isAscendingDefCom_D1_CK2(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if CK1_ENABLED
    PURE module function isAscendingDefCom_D1_CK1(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module function isAscendingDefCom_D1_PSSK5(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if SK4_ENABLED
    PURE module function isAscendingDefCom_D1_PSSK4(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if SK3_ENABLED
    PURE module function isAscendingDefCom_D1_PSSK3(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if SK2_ENABLED
    PURE module function isAscendingDefCom_D1_PSSK2(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#if SK1_ENABLED
    PURE module function isAscendingDefCom_D1_PSSK1(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function isAscendingDefCom_D1_BSSK(array) result(ascending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAscendingDefCom_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: ascending
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input array is sorted in descending order (with the possibility of elements being equal),
    !>  otherwise, generate and return `.false.`.
    !>
    !>  \param[in]  array       :   The input `contiguous` array of rank `1` of either<br>
    !>                              <ol>
    !>                                  <li>    type [css_pdt](@ref pm_container::css_pdt) or,<br>
    !>                                  <li>    type [css_type](@ref pm_container::css_type) or,<br>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter or,<br>
    !>                                  <li>    type `integer` of kind \IKALL or,<br>
    !>                                  <li>    type `logical` of kind \LKALL or,<br>
    !>                                  <li>    type `complex` of kind \CKALL or,<br>
    !>                                  <li>    type `real` of kind \RKALL or,<br>
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    a *scalar* of type `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                              </ol>
    !>                              whose elements will be checked for an all descending order<br>
    !>                              (with the possibility of elements being equal).
    !>
    !>  \return
    !>  `sorted`                :   An output `logical` of default kind \LK that is `.true.` if all elements of the input array are
    !>                              in descending order (with the possibility of elements being equal) or all equal. Otherwise, it is `.false.`.
    !>
    !>  \interface{isDescending}
    !>  \code{.F90}
    !>
    !>      use pm_arraySort, only: isDescending
    !>
    !>      descending = isDescending(array)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The output of this procedure is `.true.` when the input array has zero length.
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [isSorted](@ref isSorted)<br>
    !>  [isAscending](@ref pm_arraySort::isAscending)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>
    !>  \example{isDescending}
    !>  \include{lineno} example/pm_arraySort/isDescending/main.F90
    !>  \compilef{isDescending}
    !>  \output{isDescending}
    !>  \include{lineno} example/pm_arraySort/isDescending/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arraySort](@ref test_pm_arraySort)
    !>
    !>  \final{isDescending}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 3:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isDescending

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function isDescendingDefCom_D0_SK5(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: descending
    end function
#endif

#if SK4_ENABLED
    PURE module function isDescendingDefCom_D0_SK4(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: descending
    end function
#endif

#if SK3_ENABLED
    PURE module function isDescendingDefCom_D0_SK3(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: descending
    end function
#endif

#if SK2_ENABLED
    PURE module function isDescendingDefCom_D0_SK2(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: descending
    end function
#endif

#if SK1_ENABLED
    PURE module function isDescendingDefCom_D0_SK1(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: descending
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function isDescendingDefCom_D1_SK5(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if SK4_ENABLED
    PURE module function isDescendingDefCom_D1_SK4(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if SK3_ENABLED
    PURE module function isDescendingDefCom_D1_SK3(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if SK2_ENABLED
    PURE module function isDescendingDefCom_D1_SK2(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if SK1_ENABLED
    PURE module function isDescendingDefCom_D1_SK1(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function isDescendingDefCom_D1_IK5(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if IK4_ENABLED
    PURE module function isDescendingDefCom_D1_IK4(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if IK3_ENABLED
    PURE module function isDescendingDefCom_D1_IK3(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if IK2_ENABLED
    PURE module function isDescendingDefCom_D1_IK2(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if IK1_ENABLED
    PURE module function isDescendingDefCom_D1_IK1(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function isDescendingDefCom_D1_LK5(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if LK4_ENABLED
    PURE module function isDescendingDefCom_D1_LK4(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if LK3_ENABLED
    PURE module function isDescendingDefCom_D1_LK3(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if LK2_ENABLED
    PURE module function isDescendingDefCom_D1_LK2(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if LK1_ENABLED
    PURE module function isDescendingDefCom_D1_LK1(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function isDescendingDefCom_D1_RK5(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if RK4_ENABLED
    PURE module function isDescendingDefCom_D1_RK4(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if RK3_ENABLED
    PURE module function isDescendingDefCom_D1_RK3(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if RK2_ENABLED
    PURE module function isDescendingDefCom_D1_RK2(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if RK1_ENABLED
    PURE module function isDescendingDefCom_D1_RK1(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function isDescendingDefCom_D1_CK5(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if CK4_ENABLED
    PURE module function isDescendingDefCom_D1_CK4(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if CK3_ENABLED
    PURE module function isDescendingDefCom_D1_CK3(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if CK2_ENABLED
    PURE module function isDescendingDefCom_D1_CK2(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if CK1_ENABLED
    PURE module function isDescendingDefCom_D1_CK1(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module function isDescendingDefCom_D1_PSSK5(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if SK4_ENABLED
    PURE module function isDescendingDefCom_D1_PSSK4(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if SK3_ENABLED
    PURE module function isDescendingDefCom_D1_PSSK3(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if SK2_ENABLED
    PURE module function isDescendingDefCom_D1_PSSK2(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#if SK1_ENABLED
    PURE module function isDescendingDefCom_D1_PSSK1(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function isDescendingDefCom_D1_BSSK(array) result(descending)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isDescendingDefCom_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: descending
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return `.true.` if the input array is sorted, either ascending or descending, or all equal.
    !>
    !>  \param[in]  array       :   The input `contiguous` array of rank `1` of either<br>
    !>                              <ol>
    !>                                  <li>    type [css_pdt](@ref pm_container::css_pdt) or,<br>
    !>                                  <li>    type [css_type](@ref pm_container::css_type) or,<br>
    !>                                  <li>    type `character` of kind \SKALL of arbitrary length type parameter or,<br>
    !>                                  <li>    type `integer` of kind \IKALL or,<br>
    !>                                  <li>    type `logical` of kind \LKALL or,<br>
    !>                                  <li>    type `complex` of kind \CKALL or,<br>
    !>                                  <li>    type `real` of kind \RKALL or,<br>
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    a *scalar* of type `character` of kind \SKALL of arbitrary length type parameter,<br>
    !>                              </ol>
    !>                              whose elements will be checked to have an all ascending, or all descending order with the possibility of elements being equal.
    !>  \param      isSorted    :   The `external` user-specified function that takes two input **scalar** arguments of the same type and kind as the input `array`.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if the first
    !>                              input scalar argument is sorted with respect to the second input argument according to the user-defined condition
    !>                              within `isSorted`, otherwise, it is `.false.`.<br>
    !>                              If `array` is a scalar string (i.e., an assumed-length scalar `character`),
    !>                              then both input arguments to `isSorted()` are scalar characters of length `1` of kind \SKALL.<br>
    !>                              The following illustrates the generic interface of `isSorted()`,
    !>                              \code{.F90}
    !>                                  function isSorted(lhs, rhs) result(sorted)
    !>                                      use pm_kind, only: LK
    !>                                      TYPE(KIND)  , intent(in)    :: lhs, rhs
    !>                                      logical(LK)                 :: sorted
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` is the same as the type and kind of the input argument `array`, which can be one of the following.
    !>                              \code{.F90}
    !>                                  use pm_container, only: css_type, css_pdt
    !>                                  character(*, SK), intent(in) :: lhs, rhs
    !>                                  character(1, SK), intent(in) :: lhs, rhs
    !>                                  type(css_type)  , intent(in) :: lhs, rhs
    !>                                  type(css_pdt)   , intent(in) :: lhs, rhs
    !>                                  integer(IK)     , intent(in) :: lhs, rhs
    !>                                  logical(LK)     , intent(in) :: lhs, rhs
    !>                                  complex(CK)     , intent(in) :: lhs, rhs
    !>                                  real(RK)        , intent(in) :: lhs, rhs
    !>                              \endcode
    !>                              where the specified kind type parameters (`SK`, `IK`, `LK`, `CK`, `RK`) can refer to any of the supported kinds by the processor.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined sorting criterion other than simple ascending order
    !>                              is needed, for example, when the case-sensitivity of an input string or array of strings is irrelevant or when sorting of
    !>                              the absolute values matters excluding the signs of the numbers, or when descending order is desired.<br>
    !>                              In such cases, user can define a custom sorting condition within the user-defined external function `isSorted` to achieve the goal.<br>
    !>                              (**optional**, the default sorting condition is ascending order, that is `a < b`.)
    !>
    !>  \return
    !>  `sorted`                :   An output `logical` of default kind \LK that is `.true.` if all elements of the input array are
    !>                              in ascending order, in descending order, or all equal. Otherwise, it is `.false.`.
    !>
    !>  \interface{isSorted}
    !>  \code{.F90}
    !>
    !>      use pm_arraySort, only: isSorted
    !>
    !>      sorted = isSorted(array)
    !>      sorted = isSorted(array, isSorted)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The output of this procedure is `.true.` when the input array has zero length.
    !>
    !>  \warnpure
    !>  The procedures under this generic interface are always `impure` when the input argument `isSorted()` is present.
    !>
    !>  \see
    !>  [isAscending](@ref pm_arraySort::isAscending)<br>
    !>  [isDescending](@ref pm_arraySort::isDescending)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>
    !>  \example{isSorted}
    !>  \include{lineno} example/pm_arraySort/isSorted/main.F90
    !>  \compilef{isSorted}
    !>  \output{isSorted}
    !>  \include{lineno} example/pm_arraySort/isSorted/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arraySort](@ref test_pm_arraySort)
    !>
    !>  \final{isSorted}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 3:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface isSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function isSortedDefCom_D0_SK5(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: sorted
    end function
#endif

#if SK4_ENABLED
    PURE module function isSortedDefCom_D0_SK4(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: sorted
    end function
#endif

#if SK3_ENABLED
    PURE module function isSortedDefCom_D0_SK3(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: sorted
    end function
#endif

#if SK2_ENABLED
    PURE module function isSortedDefCom_D0_SK2(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: sorted
    end function
#endif

#if SK1_ENABLED
    PURE module function isSortedDefCom_D0_SK1(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        logical(LK)                                             :: sorted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function isSortedDefCom_D1_SK5(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if SK4_ENABLED
    PURE module function isSortedDefCom_D1_SK4(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if SK3_ENABLED
    PURE module function isSortedDefCom_D1_SK3(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if SK2_ENABLED
    PURE module function isSortedDefCom_D1_SK2(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if SK1_ENABLED
    PURE module function isSortedDefCom_D1_SK1(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function isSortedDefCom_D1_IK5(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if IK4_ENABLED
    PURE module function isSortedDefCom_D1_IK4(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if IK3_ENABLED
    PURE module function isSortedDefCom_D1_IK3(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if IK2_ENABLED
    PURE module function isSortedDefCom_D1_IK2(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if IK1_ENABLED
    PURE module function isSortedDefCom_D1_IK1(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function isSortedDefCom_D1_LK5(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if LK4_ENABLED
    PURE module function isSortedDefCom_D1_LK4(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if LK3_ENABLED
    PURE module function isSortedDefCom_D1_LK3(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if LK2_ENABLED
    PURE module function isSortedDefCom_D1_LK2(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if LK1_ENABLED
    PURE module function isSortedDefCom_D1_LK1(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function isSortedDefCom_D1_RK5(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if RK4_ENABLED
    PURE module function isSortedDefCom_D1_RK4(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if RK3_ENABLED
    PURE module function isSortedDefCom_D1_RK3(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if RK2_ENABLED
    PURE module function isSortedDefCom_D1_RK2(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if RK1_ENABLED
    PURE module function isSortedDefCom_D1_RK1(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function isSortedDefCom_D1_CK5(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if CK4_ENABLED
    PURE module function isSortedDefCom_D1_CK4(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if CK3_ENABLED
    PURE module function isSortedDefCom_D1_CK3(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if CK2_ENABLED
    PURE module function isSortedDefCom_D1_CK2(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if CK1_ENABLED
    PURE module function isSortedDefCom_D1_CK1(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module function isSortedDefCom_D1_PSSK5(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if SK4_ENABLED
    PURE module function isSortedDefCom_D1_PSSK4(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if SK3_ENABLED
    PURE module function isSortedDefCom_D1_PSSK3(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if SK2_ENABLED
    PURE module function isSortedDefCom_D1_PSSK2(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#if SK1_ENABLED
    PURE module function isSortedDefCom_D1_PSSK1(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function isSortedDefCom_D1_BSSK(array) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedDefCom_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(in)    , contiguous    :: array(:)
        logical(LK)                                             :: sorted
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function isSortedCusCom_D0_SK5(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if SK4_ENABLED
    module function isSortedCusCom_D0_SK4(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if SK3_ENABLED
    module function isSortedCusCom_D0_SK3(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if SK2_ENABLED
    module function isSortedCusCom_D0_SK2(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if SK1_ENABLED
    module function isSortedCusCom_D0_SK1(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function isSortedCusCom_D1_SK5(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if SK4_ENABLED
    module function isSortedCusCom_D1_SK4(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if SK3_ENABLED
    module function isSortedCusCom_D1_SK3(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if SK2_ENABLED
    module function isSortedCusCom_D1_SK2(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if SK1_ENABLED
    module function isSortedCusCom_D1_SK1(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function isSortedCusCom_D1_IK5(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if IK4_ENABLED
    module function isSortedCusCom_D1_IK4(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if IK3_ENABLED
    module function isSortedCusCom_D1_IK3(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if IK2_ENABLED
    module function isSortedCusCom_D1_IK2(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if IK1_ENABLED
    module function isSortedCusCom_D1_IK1(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function isSortedCusCom_D1_LK5(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if LK4_ENABLED
    module function isSortedCusCom_D1_LK4(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if LK3_ENABLED
    module function isSortedCusCom_D1_LK3(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if LK2_ENABLED
    module function isSortedCusCom_D1_LK2(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if LK1_ENABLED
    module function isSortedCusCom_D1_LK1(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function isSortedCusCom_D1_RK5(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if RK4_ENABLED
    module function isSortedCusCom_D1_RK4(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if RK3_ENABLED
    module function isSortedCusCom_D1_RK3(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if RK2_ENABLED
    module function isSortedCusCom_D1_RK2(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if RK1_ENABLED
    module function isSortedCusCom_D1_RK1(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function isSortedCusCom_D1_CK5(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if CK4_ENABLED
    module function isSortedCusCom_D1_CK4(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if CK3_ENABLED
    module function isSortedCusCom_D1_CK3(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if CK2_ENABLED
    module function isSortedCusCom_D1_CK2(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if CK1_ENABLED
    module function isSortedCusCom_D1_CK1(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function isSortedCusCom_D1_PSSK5(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if SK4_ENABLED
    module function isSortedCusCom_D1_PSSK4(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if SK3_ENABLED
    module function isSortedCusCom_D1_PSSK3(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if SK2_ENABLED
    module function isSortedCusCom_D1_PSSK2(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#if SK1_ENABLED
    module function isSortedCusCom_D1_PSSK1(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function isSortedCusCom_D1_BSSK(array, isSorted) result(sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isSortedCusCom_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        logical(LK)                                             :: sorted
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the sorted elements of the input scalar string or `contiguous` vector in ascending order, or
    !>  the **sorted indices** of the input scalar string or `contiguous` `array` of rank `1` in **ascending order**
    !>  or in the user-specified order.
    !>
    !>  \details
    !>  This generic interface either sorts the contents of an array or its indices in ascending order or via the
    !>  user-specified custom input procedure `isSorted()`.<br>
    !>
    !>  The resulting output array or its `sorting` will be in ascending order (or in the requested order as specified by `isSorted()`.<br>
    !>
    !>  \note
    !>  This generic interface is merely a convenient functional wrapper for the
    !>  lower-level subroutine interface [setSorted](@ref pm_arraySort::setSorted).<br>
    !>
    !>  \param[in]      array       :   The `contiguous` array of rank `1` of either<br>
    !>                                  <ol>
    !>                                      <li>    type [css_pdt](@ref pm_container::css_pdt) (string container) or,<br>
    !>                                      <li>    type [css_type](@ref pm_container::css_type) (string container of default kind) or,<br>
    !>                                      <li>    type `character` of kind \SKALL of arbitrary length type parameter or,
    !>                                      <li>    type `integer` of kind \IKALL or,<br>
    !>                                      <li>    type `logical` of kind \LKALL or,<br>
    !>                                      <li>    type `complex` of kind \CKALL or,<br>
    !>                                      <li>    type `real` of kind \RKALL,<br>
    !>                                  </ol>
    !>                                  or,
    !>                                  <ol>
    !>                                      <li>    a **scalar** of type `character` of kind \SKALL of arbitrary length type parameter.<br>
    !>                                  </ol>
    !>  \param          isSorted    :   The `external` user-specified function that takes two input **scalar** arguments of the same type and kind as the input `array`.<br>
    !>                                  It returns a scalar `logical` of default kind \LK that is `.true.` if the first
    !>                                  input scalar argument is sorted with respect to the second input argument according to the user-defined condition
    !>                                  within `isSorted`, otherwise, it is `.false.`.<br>
    !>                                  If `array` is a scalar string (i.e., an assumed-length scalar `character`),
    !>                                  then both input arguments to `isSorted()` are scalar characters of length `1` of kind \SKALL.<br>
    !>                                  The following illustrates the generic interface of `isSorted()`,
    !>                                  \code{.F90}
    !>                                      function isSorted(lhs, rhs) result(sorted)
    !>                                          use pm_kind, only: LK
    !>                                          TYPE(KIND)  , intent(in)    :: lhs, rhs
    !>                                          logical(LK)                 :: sorted
    !>                                      end function
    !>                                  \endcode
    !>                                  where `TYPE(KIND)` is the same as the type and kind of the input argument `array`, which can be one of the following.
    !>                                  \code{.F90}
    !>                                      use pm_container, only: css_type, css_pdt
    !>                                      character(*, SK), intent(in) :: lhs, rhs
    !>                                      character(1, SK), intent(in) :: lhs, rhs
    !>                                      type(css_type)  , intent(in) :: lhs, rhs
    !>                                      type(css_pdt)   , intent(in) :: lhs, rhs
    !>                                      integer(IK)     , intent(in) :: lhs, rhs
    !>                                      logical(LK)     , intent(in) :: lhs, rhs
    !>                                      complex(CK)     , intent(in) :: lhs, rhs
    !>                                      real(RK)        , intent(in) :: lhs, rhs
    !>                                  \endcode
    !>                                  where the specified kind type parameters (`SK`, `IK`, `LK`, `CK`, `RK`) can refer to any of the supported kinds by the processor.<br>
    !>                                  This user-defined equivalence check is extremely useful where a user-defined sorting criterion other than simple ascending order
    !>                                  is needed, for example, when the case-sensitivity of an input string or array of strings is irrelevant or when sorting of
    !>                                  the absolute values matters excluding the signs of the numbers, or when descending order is desired.<br>
    !>                                  In such cases, user can define a custom sorting condition within the user-defined external function `isSorted` to achieve the goal.<br>
    !>                                  (**optional**, the default sorting condition is ascending order, that is `a < b`.)
    !>  \param          method      :   The input scalar constant that can be any of the following:<br>
    !>                                  <ol>
    !>                                      <li>    The constant [isort](@ref pm_arraySort::isort) or equivalently, an object of type [isort_type](@ref pm_arraySort::isort_type),
    !>                                              implying that the **sorted indices** of the input array must be returned (instead of the sorted elements of array).<br>
    !>                                      <li>    The constant [qsorti](@ref pm_arraySort::qsorti) or equivalently, an object of type [qsorti_type](@ref pm_arraySort::qsorti_type),
    !>                                              implying that the **iterative** version of the **QuickSort** sorting algorithm should be used.<br>
    !>                                              This method sorts the input array by a mixture of Quicksort and Selection sorting methods.<br>
    !>                                              When the size the array to be sorted reaches `30` or less, the algorithm switches from Quicksort to Selection sorting.<br>
    !>                                              This algorithm is typically of order \f$ N \log_2 N\f$, and the worst-case order of \f$N^2\f$.<br>
    !>                                              The worst case performance occurs for completely sorted input arrays.<br>
    !>                                      <li>    The constant [qsortr](@ref pm_arraySort::qsortr) or equivalently, an object of type [qsortr_type](@ref pm_arraySort::qsortr_type),
    !>                                              implying that the **recursive** version of the **QuickSort** sorting algorithm should be used.<br>
    !>                                      <li>    The constant [qsortrdp](@ref pm_arraySort::qsortrdp) or equivalently, an object of type [qsortrdp_type](@ref pm_arraySort::qsortrdp_type),
    !>                                              implying that the **recursive** version of the **Dual-Pivot QuickSort** sorting algorithm should be used.<br>
    !>                                              The Dual-Pivot Quicksort algorithm can be slightly faster than the default Quicksort algorithms above.<br>
    !>                                              However, [performance benchmarks](@ref pm_arraySort) indicate that the efficiency gain is marginal and insignificant.<br>
    !>                                              This algorithm is typically of order \f$ N \log_2 N\f$, and the worst-case order of \f$N^2\f$.
    !>                                              The worst case performance occurs for completely sorted input arrays.<br>
    !>                                      <li>    The constant [bubble](@ref pm_arraySort::bubble) or equivalently, an object of type [bubble_type](@ref pm_arraySort::bubble_type),
    !>                                              implying that the **Bubble** sorting algorithm should be used.<br>
    !>                                              This algorithm is of order \f$N^2\f$.<br>
    !>                                      <li>    The constant [heapi](@ref pm_arraySort::heapi) or equivalently, an object of type [heapi_type](@ref pm_arraySort::heapi_type),
    !>                                              implying that the **iterative** version of the **Heap** sorting algorithm should be used.<br>
    !>                                              This algorithm is typically of order \f$ N \log_2 N\f$.<br>
    !>                                      <li>    The constant [heapr](@ref pm_arraySort::heapr) or equivalently, an object of type [heapr_type](@ref pm_arraySort::heapr_type),
    !>                                              implying that the **recursive** version of the **Heap** sorting algorithm should be used.<br>
    !>                                              This algorithm is typically of order \f$ N \log_2 N\f$.<br>
    !>                                      <li>    The constant [insertionl](@ref pm_arraySort::insertionl) or equivalently, an object of type [insertionl_type](@ref pm_arraySort::insertionl_type),
    !>                                              implying that the **linear-search** version of the **insertion** sorting algorithm should be used.<br>
    !>                                              The complexity of this algorithm is typically of order \f$N^2\f$.<br>
    !>                                      <li>    The constant [insertionb](@ref pm_arraySort::insertionb) or equivalently, an object of type [insertionb_type](@ref pm_arraySort::insertionb_type),
    !>                                              implying that the **binary-search** version of the **insertion** sorting algorithm should be used.<br>
    !>                                              The complexity of this algorithm is typically of order \f$N^2\f$.<br>
    !>                                      <li>    The constant [merger](@ref pm_arraySort::merger) or equivalently, an object of type [merger_type](@ref pm_arraySort::merger_type),
    !>                                              implying that the **recursive** version of the **Merge** sorting algorithm should be used.<br>
    !>                                              This algorithm is typically of order \f$ N \log_2 N\f$.<br>
    !>                                      <li>    The constant [selection](@ref pm_arraySort::selection) or equivalently, an object of type [selection_type](@ref pm_arraySort::selection_type),
    !>                                              implying that the **Selection** sorting algorithm should be used.<br>
    !>                                              This algorithm is of order \f$N^2\f$.<br>
    !>                                      <li>    The constant [shell](@ref pm_arraySort::shell) or equivalently, an object of type [shell_type](@ref pm_arraySort::shell_type),
    !>                                              implying that the **Shell** sorting algorithm should be used.<br>
    !>                                              This algorithm is of order \f$N\log(N)\f$.<br>
    !>                                  </ol>
    !>                                  The presence of this argument is merely for compile-time resolution of the procedures of this generic interface.<br>
    !>                                  (**optional**. default = [qsorti](@ref pm_arraySort::qsorti). It can be present only if the input argument `sorting` is missing.)
    !>
    !>  \return
    !>  `sorting`                   :   The output object whose type and kind depends on the specified input `method`:<br>
    !>                                  <ol>
    !>                                      <li>    If `method` is set to [isort](@ref pm_arraySort::isort) or
    !>                                              an object of type set to [isort_type](@ref pm_arraySort::isort_type),
    !>                                              then the output `sorting` is of type `integer` of default kind \IK of rank `1`
    !>                                              of the same size as the number of elements in the input `array` containing the
    !>                                              **sorted indices** of the elements of the input array`.<br>
    !>                                              **Read `sorting(i)` as the index of the element of `array` that contains the `i`th smallest (or ranked) value in the `array`.**<br>
    !>                                              This kind of sorting of indices is also widely known as **ordinal ranking**.<br>
    !>                                              In ordinal ranking, all items receive distinct ordinal numbers, including items that compare equal.<br>
    !>                                              The assignment of distinct ordinal numbers to items that compare equal can be done at random, or arbitrarily.<br>
    !>                                              But it is generally preferable to use a system that is arbitrary but consistent,
    !>                                              as this gives stable results if the ranking is done multiple times.<br>
    !>                                              In computer data processing, ordinal ranking is also referred to as **row numbering**.<br>
    !>                                      <li>    If `method` is missing or set to any other value than the above,
    !>                                              then the output `sorting` is of the same type, kind, rank, and shape as
    !>                                              the input `array`, containing the **sorted elements** of the input array.<br>
    !>                                  </ol>
    !>
    !>  \interface{getSorted}
    !>  \code{.F90}
    !>
    !>      use pm_arraySort, only: getSorted
    !>
    !>      ! Sorting the indices of an array.
    !>
    !>      sorting(1:len(array)) = getSorted(array, isort, isSorted = isSorted) ! scalar characer `array`.
    !>      sorting(1:size(array)) = getSorted(array(:), isort, isSorted = isSorted) ! all other intrinsic data types.
    !>
    !>      ! Sorting the contents of an array.
    !>
    !>      sorting(1:len(array)) = getSorted(array, method = method, isSorted = isSorted) ! scalar characer `array`.
    !>      sorting(1:size(array)) = getSorted(array(:), method = method, isSorted = isSorted) ! all other intrinsic data types.
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>  The procedures under this generic interface are always `impure` when the input argument `isSorted()` is present.<br>
    !>
    !>  \devnote
    !>  The use of the Fortran intrinsic `minloc` in the sorting routines appears to lead a slightly better performance than the manual search.<br>
    !>
    !>  \note
    !>  To rank the input array in descending order using the output `sorting` indices,
    !>  simply call the output `sorting` indices from the end to the beginning or,
    !>  rewrite the array with `array = array(sorting(ubound(array, dim = 1) : lbound(array, dim = 1) : -1) )`.<br>
    !>
    !>  \note
    !>  To sort the input array in descending order, simply call the output array elements from the end to the beginning or,
    !>  rewrite the array with `array = array(ubound(array, dim = 1) : lbound(array, dim = 1) : -1)`.<br>
    !>  Alternatively, supply an external comparison function `isSorted()` with the appropriate comparison.<br>
    !>
    !>  \note
    !>  A rooter array can be sorted along with a leader array with the help of [getSorted](@ref pm_arraySort::getSorted).<br>
    !>
    !>  \see
    !>  [getLoc](@ref pm_arrayFind::getLoc)<br>
    !>  [setLoc](@ref pm_arrayFind::setLoc)<br>
    !>  [getBin](@ref pm_arraySearch::getBin)<br>
    !>  [getSorted](@ref pm_arraySort::getSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [getRankDense](@ref pm_arrayRank::getRankDense)<br>
    !>  [setRankDense](@ref pm_arrayRank::setRankDense)<br>
    !>  [getRankOrdinal](@ref pm_arrayRank::getRankOrdinal)<br>
    !>  [setRankOrdinal](@ref pm_arrayRank::setRankOrdinal)<br>
    !>  [getRankModified](@ref pm_arrayRank::getRankModified)<br>
    !>  [setRankModified](@ref pm_arrayRank::setRankModified)<br>
    !>  [getRankStandard](@ref pm_arrayRank::getRankStandard)<br>
    !>  [setRankStandard](@ref pm_arrayRank::setRankStandard)<br>
    !>  [getRankFractional](@ref pm_arrayRank::getRankFractional)<br>
    !>  [setRankFractional](@ref pm_arrayRank::setRankFractional)<br>
    !>  [getSelected](@ref pm_arraySelect::getSelected)<br>
    !>  [setSelected](@ref pm_arraySelect::setSelected)<br>
    !>
    !>  \example{getSorted}
    !>  \include{lineno} example/pm_arraySort/getSorted/main.F90
    !>  \compilef{getSorted}
    !>  \output{getSorted}
    !>  \include{lineno} example/pm_arraySort/getSorted/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arraySort](@ref test_pm_arraySort)
    !>
    !>  \final{getSorted}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! getSortedIndCusComDef

    interface getSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getSortedIndCusComDef_D0_SK5(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(len(array, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK4_ENABLED
    module function getSortedIndCusComDef_D0_SK4(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(len(array, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK3_ENABLED
    module function getSortedIndCusComDef_D0_SK3(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(len(array, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK2_ENABLED
    module function getSortedIndCusComDef_D0_SK2(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(len(array, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK1_ENABLED
    module function getSortedIndCusComDef_D0_SK1(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(len(array, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getSortedIndCusComDef_D1_SK5(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK4_ENABLED
    module function getSortedIndCusComDef_D1_SK4(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK3_ENABLED
    module function getSortedIndCusComDef_D1_SK3(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK2_ENABLED
    module function getSortedIndCusComDef_D1_SK2(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK1_ENABLED
    module function getSortedIndCusComDef_D1_SK1(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getSortedIndCusComDef_D1_IK5(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if IK4_ENABLED
    module function getSortedIndCusComDef_D1_IK4(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if IK3_ENABLED
    module function getSortedIndCusComDef_D1_IK3(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if IK2_ENABLED
    module function getSortedIndCusComDef_D1_IK2(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if IK1_ENABLED
    module function getSortedIndCusComDef_D1_IK1(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getSortedIndCusComDef_D1_LK5(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if LK4_ENABLED
    module function getSortedIndCusComDef_D1_LK4(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if LK3_ENABLED
    module function getSortedIndCusComDef_D1_LK3(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if LK2_ENABLED
    module function getSortedIndCusComDef_D1_LK2(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if LK1_ENABLED
    module function getSortedIndCusComDef_D1_LK1(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getSortedIndCusComDef_D1_CK5(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if CK4_ENABLED
    module function getSortedIndCusComDef_D1_CK4(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if CK3_ENABLED
    module function getSortedIndCusComDef_D1_CK3(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if CK2_ENABLED
    module function getSortedIndCusComDef_D1_CK2(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if CK1_ENABLED
    module function getSortedIndCusComDef_D1_CK1(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getSortedIndCusComDef_D1_RK5(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    module function getSortedIndCusComDef_D1_RK4(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    module function getSortedIndCusComDef_D1_RK3(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    module function getSortedIndCusComDef_D1_RK2(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    module function getSortedIndCusComDef_D1_RK1(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getSortedIndCusComDef_D1_PSSK5(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK4_ENABLED
    module function getSortedIndCusComDef_D1_PSSK4(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK3_ENABLED
    module function getSortedIndCusComDef_D1_PSSK3(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK2_ENABLED
    module function getSortedIndCusComDef_D1_PSSK2(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK1_ENABLED
    module function getSortedIndCusComDef_D1_PSSK1(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getSortedIndCusComDef_D1_BSSK(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndCusComDef_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        procedure(logical(LK))                                      :: isSorted
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! getSortedArrCusComDef

    interface getSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getSortedArrCusComDef_D0_SK5(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        character(len(array, IK),SKG)                           :: sorting
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK4_ENABLED
    module function getSortedArrCusComDef_D0_SK4(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        character(len(array, IK),SKG)                           :: sorting
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK3_ENABLED
    module function getSortedArrCusComDef_D0_SK3(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        character(len(array, IK),SKG)                           :: sorting
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK2_ENABLED
    module function getSortedArrCusComDef_D0_SK2(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        character(len(array, IK),SKG)                           :: sorting
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK1_ENABLED
    module function getSortedArrCusComDef_D0_SK1(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        character(len(array, IK),SKG)                           :: sorting
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getSortedArrCusComDef_D1_SK5(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(len(array, IK),SKG)                           :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK4_ENABLED
    module function getSortedArrCusComDef_D1_SK4(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(len(array, IK),SKG)                           :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK3_ENABLED
    module function getSortedArrCusComDef_D1_SK3(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(len(array, IK),SKG)                           :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK2_ENABLED
    module function getSortedArrCusComDef_D1_SK2(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(len(array, IK),SKG)                           :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK1_ENABLED
    module function getSortedArrCusComDef_D1_SK1(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(len(array, IK),SKG)                           :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getSortedArrCusComDef_D1_IK5(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                            :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if IK4_ENABLED
    module function getSortedArrCusComDef_D1_IK4(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                            :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if IK3_ENABLED
    module function getSortedArrCusComDef_D1_IK3(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                            :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if IK2_ENABLED
    module function getSortedArrCusComDef_D1_IK2(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                            :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if IK1_ENABLED
    module function getSortedArrCusComDef_D1_IK1(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                            :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getSortedArrCusComDef_D1_LK5(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                            :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if LK4_ENABLED
    module function getSortedArrCusComDef_D1_LK4(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                            :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if LK3_ENABLED
    module function getSortedArrCusComDef_D1_LK3(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                            :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if LK2_ENABLED
    module function getSortedArrCusComDef_D1_LK2(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                            :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if LK1_ENABLED
    module function getSortedArrCusComDef_D1_LK1(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                            :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getSortedArrCusComDef_D1_CK5(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                            :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if CK4_ENABLED
    module function getSortedArrCusComDef_D1_CK4(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                            :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if CK3_ENABLED
    module function getSortedArrCusComDef_D1_CK3(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                            :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if CK2_ENABLED
    module function getSortedArrCusComDef_D1_CK2(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                            :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if CK1_ENABLED
    module function getSortedArrCusComDef_D1_CK1(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                            :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getSortedArrCusComDef_D1_RK5(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)                                               :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if RK4_ENABLED
    module function getSortedArrCusComDef_D1_RK4(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)                                               :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if RK3_ENABLED
    module function getSortedArrCusComDef_D1_RK3(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)                                               :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if RK2_ENABLED
    module function getSortedArrCusComDef_D1_RK2(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)                                               :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if RK1_ENABLED
    module function getSortedArrCusComDef_D1_RK1(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)                                               :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getSortedArrCusComDef_D1_PSSK5(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                      :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK4_ENABLED
    module function getSortedArrCusComDef_D1_PSSK4(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                      :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK3_ENABLED
    module function getSortedArrCusComDef_D1_PSSK3(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                      :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK2_ENABLED
    module function getSortedArrCusComDef_D1_PSSK2(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                      :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK1_ENABLED
    module function getSortedArrCusComDef_D1_PSSK1(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                      :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getSortedArrCusComDef_D1_BSSK(array, isSorted, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrCusComDef_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(in)    , contiguous    :: array(:)
        type(css_type)                                          :: sorting(size(array, 1, IK))
        procedure(logical(LK))                                  :: isSorted
        class(sort_type)        , intent(in)    , optional      :: method
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! getSortedIndDefComDef

    interface getSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getSortedIndDefComDef_D0_SK5(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                                                 :: sorting(len(array, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK4_ENABLED
    module function getSortedIndDefComDef_D0_SK4(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                                                 :: sorting(len(array, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK3_ENABLED
    module function getSortedIndDefComDef_D0_SK3(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                                                 :: sorting(len(array, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK2_ENABLED
    module function getSortedIndDefComDef_D0_SK2(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                                                 :: sorting(len(array, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK1_ENABLED
    module function getSortedIndDefComDef_D0_SK1(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                                                 :: sorting(len(array, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getSortedIndDefComDef_D1_SK5(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK4_ENABLED
    module function getSortedIndDefComDef_D1_SK4(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK3_ENABLED
    module function getSortedIndDefComDef_D1_SK3(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK2_ENABLED
    module function getSortedIndDefComDef_D1_SK2(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK1_ENABLED
    module function getSortedIndDefComDef_D1_SK1(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getSortedIndDefComDef_D1_IK5(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if IK4_ENABLED
    module function getSortedIndDefComDef_D1_IK4(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if IK3_ENABLED
    module function getSortedIndDefComDef_D1_IK3(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if IK2_ENABLED
    module function getSortedIndDefComDef_D1_IK2(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if IK1_ENABLED
    module function getSortedIndDefComDef_D1_IK1(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getSortedIndDefComDef_D1_LK5(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if LK4_ENABLED
    module function getSortedIndDefComDef_D1_LK4(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if LK3_ENABLED
    module function getSortedIndDefComDef_D1_LK3(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if LK2_ENABLED
    module function getSortedIndDefComDef_D1_LK2(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if LK1_ENABLED
    module function getSortedIndDefComDef_D1_LK1(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getSortedIndDefComDef_D1_CK5(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if CK4_ENABLED
    module function getSortedIndDefComDef_D1_CK4(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if CK3_ENABLED
    module function getSortedIndDefComDef_D1_CK3(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if CK2_ENABLED
    module function getSortedIndDefComDef_D1_CK2(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if CK1_ENABLED
    module function getSortedIndDefComDef_D1_CK1(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getSortedIndDefComDef_D1_RK5(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if RK4_ENABLED
    module function getSortedIndDefComDef_D1_RK4(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if RK3_ENABLED
    module function getSortedIndDefComDef_D1_RK3(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if RK2_ENABLED
    module function getSortedIndDefComDef_D1_RK2(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if RK1_ENABLED
    module function getSortedIndDefComDef_D1_RK1(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getSortedIndDefComDef_D1_PSSK5(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK4_ENABLED
    module function getSortedIndDefComDef_D1_PSSK4(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK3_ENABLED
    module function getSortedIndDefComDef_D1_PSSK3(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK2_ENABLED
    module function getSortedIndDefComDef_D1_PSSK2(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#if SK1_ENABLED
    module function getSortedIndDefComDef_D1_PSSK1(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getSortedIndDefComDef_D1_BSSK(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedIndDefComDef_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        integer(IK)                                                 :: sorting(size(array, 1, IK))
        type(isort_type)            , intent(in)                    :: method
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! getSortedArrDefComDef

    interface getSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getSortedArrDefComDef_D0_SK5(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: array
        character(len(array, IK),SKG)                           :: sorting
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK4_ENABLED
    module function getSortedArrDefComDef_D0_SK4(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: array
        character(len(array, IK),SKG)                           :: sorting
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK3_ENABLED
    module function getSortedArrDefComDef_D0_SK3(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: array
        character(len(array, IK),SKG)                           :: sorting
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK2_ENABLED
    module function getSortedArrDefComDef_D0_SK2(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: array
        character(len(array, IK),SKG)                           :: sorting
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK1_ENABLED
    module function getSortedArrDefComDef_D0_SK1(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: array
        character(len(array, IK),SKG)                           :: sorting
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getSortedArrDefComDef_D1_SK5(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(len(array, IK),SKG)                           :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK4_ENABLED
    module function getSortedArrDefComDef_D1_SK4(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(len(array, IK),SKG)                           :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK3_ENABLED
    module function getSortedArrDefComDef_D1_SK3(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(len(array, IK),SKG)                           :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK2_ENABLED
    module function getSortedArrDefComDef_D1_SK2(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(len(array, IK),SKG)                           :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK1_ENABLED
    module function getSortedArrDefComDef_D1_SK1(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: array(:)
        character(len(array, IK),SKG)                           :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getSortedArrDefComDef_D1_IK5(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                            :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if IK4_ENABLED
    module function getSortedArrDefComDef_D1_IK4(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                            :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if IK3_ENABLED
    module function getSortedArrDefComDef_D1_IK3(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                            :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if IK2_ENABLED
    module function getSortedArrDefComDef_D1_IK2(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                            :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if IK1_ENABLED
    module function getSortedArrDefComDef_D1_IK1(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: array(:)
        integer(IKG)                                            :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getSortedArrDefComDef_D1_LK5(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                            :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if LK4_ENABLED
    module function getSortedArrDefComDef_D1_LK4(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                            :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if LK3_ENABLED
    module function getSortedArrDefComDef_D1_LK3(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                            :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if LK2_ENABLED
    module function getSortedArrDefComDef_D1_LK2(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                            :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if LK1_ENABLED
    module function getSortedArrDefComDef_D1_LK1(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: array(:)
        logical(LKG)                                            :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getSortedArrDefComDef_D1_CK5(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                            :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if CK4_ENABLED
    module function getSortedArrDefComDef_D1_CK4(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                            :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if CK3_ENABLED
    module function getSortedArrDefComDef_D1_CK3(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                            :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if CK2_ENABLED
    module function getSortedArrDefComDef_D1_CK2(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                            :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if CK1_ENABLED
    module function getSortedArrDefComDef_D1_CK1(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: array(:)
        complex(CKG)                                            :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getSortedArrDefComDef_D1_RK5(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)                                               :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if RK4_ENABLED
    module function getSortedArrDefComDef_D1_RK4(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)                                               :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if RK3_ENABLED
    module function getSortedArrDefComDef_D1_RK3(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)                                               :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if RK2_ENABLED
    module function getSortedArrDefComDef_D1_RK2(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)                                               :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if RK1_ENABLED
    module function getSortedArrDefComDef_D1_RK1(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: array(:)
        real(RKG)                                               :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getSortedArrDefComDef_D1_PSSK5(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                      :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK4_ENABLED
    module function getSortedArrDefComDef_D1_PSSK4(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                      :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK3_ENABLED
    module function getSortedArrDefComDef_D1_PSSK3(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                      :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK2_ENABLED
    module function getSortedArrDefComDef_D1_PSSK2(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                      :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#if SK1_ENABLED
    module function getSortedArrDefComDef_D1_PSSK1(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: array(:)
        type(css_pdt(SKG))                                      :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getSortedArrDefComDef_D1_BSSK(array, method) result(sorting)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSortedArrDefComDef_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(in)    , contiguous    :: array(:)
        type(css_type)                                          :: sorting(size(array, 1, IK))
        class(sort_type)        , intent(in)    , optional      :: method
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Sort the input scalar string or `contiguous` vector in ascending order, or
    !>  return the **sorted indices** of the input scalar string or `contiguous` `array` of rank `1` in **ascending order**
    !>  or in the user-specified order.
    !>
    !>  \details
    !>  This generic interface either sorts the contents of an array or its indices in ascending order or via the
    !>  user-specified custom input procedure `isSorted()`.<br>
    !>
    !>  The resulting output array or its `index` will be in ascending order (or in the requested order as specified by `isSorted()`.<br>
    !>  There are currently twelve sorting algorithms implemented in this generic interface.<br>
    !>  Only the default sorting method can be used for sorting the indices of an array.<br>
    !>  However, all sorting methods can be used to sort the elements of the array.<br>
    !>  The default algorithm is generally the fastest sorting method.<br>
    !>
    !>  \param[inout]   array       :   The `contiguous` array of rank `1` of either<br>
    !>                                  <ol>
    !>                                      <li>    type [css_pdt](@ref pm_container::css_pdt) (string container) or,<br>
    !>                                      <li>    type [css_type](@ref pm_container::css_type) (string container of default kind) or,<br>
    !>                                      <li>    type `character` of kind \SKALL of arbitrary length type parameter or,
    !>                                      <li>    type `integer` of kind \IKALL or,<br>
    !>                                      <li>    type `logical` of kind \LKALL or,<br>
    !>                                      <li>    type `complex` of kind \CKALL or,<br>
    !>                                      <li>    type `real` of kind \RKALL,<br>
    !>                                  </ol>
    !>                                  or,
    !>                                  <ol>
    !>                                      <li>    a **scalar** of type `character` of kind \SKALL of arbitrary length type parameter.<br>
    !>                                  </ol>
    !>                                  On output,
    !>                                  <ol>
    !>                                      <li>    If the input argument `index` is present, then `array` has `intent(in)` and its contents remain intact upon return.<br>
    !>                                      <li>    If the input argument `index` is missing, then `array` has `intent(inout)` and its contents will be in ascending order upon return.<br>
    !>                                  </ol>
    !>  \param[out]     index       :   The output `contiguous` array of rank `1` of type `integer` of kind \IKALL of the same length as the input argument `array`
    !>                                  containing the rank (i.e., sorted indices) of `array` such that `array(index(:))` is sorted in ascending order on return.<br>
    !>                                  The size of `index` must match that of `array` (or its length type parameter if `array` is a scalar string).<br>
    !>                                  **Read `index(i)` as the index of the element of `array` that contains the `i`th smallest (or ranked) value in the `array`.**<br>
    !>                                  This kind of ranking of values is widely known as **ordinal ranking**.<br>
    !>                                  In ordinal ranking, all items receive distinct ordinal numbers, including items that compare equal.<br>
    !>                                  The assignment of distinct ordinal numbers to items that compare equal can be done at random, or arbitrarily.<br>
    !>                                  But it is generally preferable to use a system that is arbitrary but consistent,
    !>                                  as this gives stable results if the ranking is done multiple times.<br>
    !>                                  In computer data processing, ordinal ranking is also referred to as **row numbering**.<br>
    !>                                  (**optional**. If missing, the contents of the input `array` will be sorted in ascending order and returned instead.
    !>                                  It can be present only if the input argument `method` is missing.)
    !>  \param          isSorted    :   The `external` user-specified function that takes two input **scalar** arguments of the same type and kind as the input `array`.<br>
    !>                                  It returns a scalar `logical` of default kind \LK that is `.true.` if the first
    !>                                  input scalar argument is sorted with respect to the second input argument according to the user-defined condition
    !>                                  within `isSorted`, otherwise, it is `.false.`.<br>
    !>                                  If `array` is a scalar string (i.e., an assumed-length scalar `character`),
    !>                                  then both input arguments to `isSorted()` are scalar characters of length `1` of kind \SKALL.<br>
    !>                                  The following illustrates the generic interface of `isSorted()`,
    !>                                  \code{.F90}
    !>                                      function isSorted(lhs, rhs) result(sorted)
    !>                                          use pm_kind, only: LK
    !>                                          TYPE(KIND)  , intent(in)    :: lhs, rhs
    !>                                          logical(LK)                 :: sorted
    !>                                      end function
    !>                                  \endcode
    !>                                  where `TYPE(KIND)` is the same as the type and kind of the input argument `array`, which can be one of the following.
    !>                                  \code{.F90}
    !>                                      use pm_container, only: css_type, css_pdt
    !>                                      character(*, SK), intent(in) :: lhs, rhs
    !>                                      character(1, SK), intent(in) :: lhs, rhs
    !>                                      type(css_type)  , intent(in) :: lhs, rhs
    !>                                      type(css_pdt)   , intent(in) :: lhs, rhs
    !>                                      integer(IK)     , intent(in) :: lhs, rhs
    !>                                      logical(LK)     , intent(in) :: lhs, rhs
    !>                                      complex(CK)     , intent(in) :: lhs, rhs
    !>                                      real(RK)        , intent(in) :: lhs, rhs
    !>                                  \endcode
    !>                                  where the specified kind type parameters (`SK`, `IK`, `LK`, `CK`, `RK`) can refer to any of the supported kinds by the processor.<br>
    !>                                  This user-defined equivalence check is extremely useful where a user-defined sorting criterion other than simple ascending order
    !>                                  is needed, for example, when the case-sensitivity of an input string or array of strings is irrelevant or when sorting of
    !>                                  the absolute values matters excluding the signs of the numbers, or when descending order is desired.<br>
    !>                                  In such cases, user can define a custom sorting condition within the user-defined external function `isSorted` to achieve the goal.<br>
    !>                                  (**optional**, the default sorting condition is ascending order, that is `a < b`.)
    !>  \param          method      :   The input scalar constant that can be any of the following:<br>
    !>                                  <ol>
    !>                                      <li>    The constant [qsorti](@ref pm_arraySort::qsorti) or equivalently, an object of type [qsorti_type](@ref pm_arraySort::qsorti_type),
    !>                                              implying that the **iterative** version of the **QuickSort** sorting algorithm should be used.<br>
    !>                                              This method sorts the input array by a mixture of Quicksort and Selection sorting methods.<br>
    !>                                              When the size the array to be sorted reaches `30` or less, the algorithm switches from Quicksort to Selection sorting.<br>
    !>                                              This algorithm is typically of order \f$ N \log_2 N\f$, and the worst-case order of \f$N^2\f$.<br>
    !>                                              The worst case performance occurs for completely sorted input arrays.<br>
    !>                                      <li>    The constant [qsortr](@ref pm_arraySort::qsortr) or equivalently, an object of type [qsortr_type](@ref pm_arraySort::qsortr_type),
    !>                                              implying that the **recursive** version of the **QuickSort** sorting algorithm should be used.<br>
    !>                                      <li>    The constant [qsortrdp](@ref pm_arraySort::qsortrdp) or equivalently, an object of type [qsortrdp_type](@ref pm_arraySort::qsortrdp_type),
    !>                                              implying that the **recursive** version of the **Dual-Pivot QuickSort** sorting algorithm should be used.<br>
    !>                                              The Dual-Pivot Quicksort algorithm can be slightly faster than the default Quicksort algorithms above.<br>
    !>                                              However, [performance benchmarks](@ref pm_arraySort) indicate that the efficiency gain is marginal and insignificant.<br>
    !>                                              This algorithm is typically of order \f$ N \log_2 N\f$, and the worst-case order of \f$N^2\f$.
    !>                                              The worst case performance occurs for completely sorted input arrays.<br>
    !>                                      <li>    The constant [bubble](@ref pm_arraySort::bubble) or equivalently, an object of type [bubble_type](@ref pm_arraySort::bubble_type),
    !>                                              implying that the **Bubble** sorting algorithm should be used.<br>
    !>                                              This algorithm is of order \f$N^2\f$.<br>
    !>                                      <li>    The constant [heapi](@ref pm_arraySort::heapi) or equivalently, an object of type [heapi_type](@ref pm_arraySort::heapi_type),
    !>                                              implying that the **iterative** version of the **Heap** sorting algorithm should be used.<br>
    !>                                              This algorithm is typically of order \f$ N \log_2 N\f$.<br>
    !>                                      <li>    The constant [heapr](@ref pm_arraySort::heapr) or equivalently, an object of type [heapr_type](@ref pm_arraySort::heapr_type),
    !>                                              implying that the **recursive** version of the **Heap** sorting algorithm should be used.<br>
    !>                                              This algorithm is typically of order \f$ N \log_2 N\f$.<br>
    !>                                      <li>    The constant [insertionl](@ref pm_arraySort::insertionl) or equivalently, an object of type [insertionl_type](@ref pm_arraySort::insertionl_type),
    !>                                              implying that the **linear-search** version of the **insertion** sorting algorithm should be used.<br>
    !>                                              The complexity of this algorithm is typically of order \f$N^2\f$.<br>
    !>                                      <li>    The constant [insertionb](@ref pm_arraySort::insertionb) or equivalently, an object of type [insertionb_type](@ref pm_arraySort::insertionb_type),
    !>                                              implying that the **binary-search** version of the **insertion** sorting algorithm should be used.<br>
    !>                                              The complexity of this algorithm is typically of order \f$N^2\f$.<br>
    !>                                      <li>    The constant [merger](@ref pm_arraySort::merger) or equivalently, an object of type [merger_type](@ref pm_arraySort::merger_type),
    !>                                              implying that the **recursive** version of the **Merge** sorting algorithm should be used.<br>
    !>                                              This algorithm is typically of order \f$ N \log_2 N\f$.<br>
    !>                                      <li>    The constant [selection](@ref pm_arraySort::selection) or equivalently, an object of type [selection_type](@ref pm_arraySort::selection_type),
    !>                                              implying that the **Selection** sorting algorithm should be used.<br>
    !>                                              This algorithm is of order \f$N^2\f$.<br>
    !>                                      <li>    The constant [shell](@ref pm_arraySort::shell) or equivalently, an object of type [shell_type](@ref pm_arraySort::shell_type),
    !>                                              implying that the **Shell** sorting algorithm should be used.<br>
    !>                                              This algorithm is of order \f$N\log(N)\f$.<br>
    !>                                  </ol>
    !>                                  The presence of this argument is merely for compile-time resolution of the procedures of this generic interface.<br>
    !>                                  (**optional**. default = [qsorti](@ref pm_arraySort::qsorti). It can be present only if the input argument `index` is missing.)
    !>
    !>  \interface{setSorted}
    !>  \code{.F90}
    !>
    !>      use pm_arraySort, only: setSorted
    !>
    !>      ! Sorting the indices of an array.
    !>
    !>      call setSorted(array, index(:))
    !>      call setSorted(array, index(:), isSorted)
    !>
    !>      ! Sorting the contents of an array.
    !>
    !>      call setSorted(array)
    !>      call setSorted(array, method)
    !>      call setSorted(array, isSorted)
    !>      call setSorted(array, isSorted, method)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The sizes of the input `array` and `index` must be equal.<br>
    !>  \vericon
    !>
    !>  \pure
    !>  The procedures under this generic interface are always `impure` when the input argument `isSorted()` is present.<br>
    !>
    !>  \recursive
    !>  The procedures are explicitly `recursive` <b>only if</b> the input argument `method`is set to any of the `recursive` sorting algorithms.<br>
    !>
    !>  \devnote
    !>  The use of the Fortran intrinsic `minloc` in the sorting routines appears to lead a slightly better performance than the manual search.<br>
    !>
    !>  \note
    !>  To rank the input array in descending order, simply call the output `index` elements from the end to the beginning or,
    !>  rewrite the array with `array = array( index(ubound(array, dim = 1) : lbound(array, dim = 1) : -1) )`.<br>
    !>
    !>  \note
    !>  To sort the input array in descending order, simply call the output array elements from the end to the beginning or,
    !>  rewrite the array with `array = array(ubound(array, dim = 1) : lbound(array, dim = 1) : -1)`.<br>
    !>  Alternatively, supply an external comparison function `isSorted()` with the appropriate comparison.<br>
    !>
    !>  \note
    !>  A rooter array can be sorted along with a leader array with the help of [setSorted](@ref pm_arraySort::setSorted).<br>
    !>
    !>  \see
    !>  [getLoc](@ref pm_arrayFind::getLoc)<br>
    !>  [setLoc](@ref pm_arrayFind::setLoc)<br>
    !>  [getBin](@ref pm_arraySearch::getBin)<br>
    !>  [getSorted](@ref pm_arraySort::getSorted)<br>
    !>  [setSorted](@ref pm_arraySort::setSorted)<br>
    !>  [isAscending](@ref pm_arraySort::isAscending)<br>
    !>  [isDescending](@ref pm_arraySort::isDescending)<br>
    !>  [isAscendingAll](@ref pm_arraySort::isAscendingAll)<br>
    !>  [isDescendingAll](@ref pm_arraySort::isDescendingAll)<br>
    !>  [isSorted](@ref pm_arraySort::isSorted)<br>
    !>  [getRankDense](@ref pm_arrayRank::getRankDense)<br>
    !>  [setRankDense](@ref pm_arrayRank::setRankDense)<br>
    !>  [getRankOrdinal](@ref pm_arrayRank::getRankOrdinal)<br>
    !>  [setRankOrdinal](@ref pm_arrayRank::setRankOrdinal)<br>
    !>  [getRankModified](@ref pm_arrayRank::getRankModified)<br>
    !>  [setRankModified](@ref pm_arrayRank::setRankModified)<br>
    !>  [getRankStandard](@ref pm_arrayRank::getRankStandard)<br>
    !>  [setRankStandard](@ref pm_arrayRank::setRankStandard)<br>
    !>  [getRankFractional](@ref pm_arrayRank::getRankFractional)<br>
    !>  [setRankFractional](@ref pm_arrayRank::setRankFractional)<br>
    !>  [getSelected](@ref pm_arraySelect::getSelected)<br>
    !>  [setSelected](@ref pm_arraySelect::setSelected)<br>
    !>
    !>  \example{setSorted}
    !>  \include{lineno} example/pm_arraySort/setSorted/main.F90
    !>  \compilef{setSorted}
    !>  \output{setSorted}
    !>  \include{lineno} example/pm_arraySort/setSorted/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arraySort](@ref test_pm_arraySort)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \ifort{2021.5}
    !>  \desc
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the description of a relevant bug in PDT
    !>  name aliasing when compiled with \ifort{2021.5} that also applies to this module.
    !>  \remedy
    !>  See [pm_arraySplit](@ref pm_arraySplit) for the remedy.<br>
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{10.3}
    !>  \desc
    !>  The \gfortran{10.3} cannot not compile the allocation of PDT string container object `temp` within the implementation of the corresponding procedures.<br>
    !>  \code{.sh}
    !>
    !>       99 |         deallocate(Temp)
    !>          |                        ^
    !>      internal compiler error: in gimplify_var_or_parm_decl, at gimplify.c:2834
    !>      Please submit a full bug report...
    !>
    !>  \endcode
    !>  \remedy
    !>  For now, all allocatable `Temp` objects are converted to automatic arrays.
    !>
    !>  \todo
    !>  \pmed
    !>  \plow The current bypass for the PDT name aliasing bug should be reverted back to PDT name aliasing once the ifort bug is resolved.
    !>
    !>  \final{setSorted}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! setSortedIndDefComDef

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSortedIndDefComDef_D0_SK5(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSortedIndDefComDef_D0_SK4(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSortedIndDefComDef_D0_SK3(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSortedIndDefComDef_D0_SK2(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSortedIndDefComDef_D0_SK1(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_SK5(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_SK4(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_SK3(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_SK2(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_SK1(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_IK5(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_IK4(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_IK3(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_IK2(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_IK1(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_LK5(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_LK4(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_LK3(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_LK2(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_LK1(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_CK5(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_CK4(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_CK3(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_CK2(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_CK1(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_RK5(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_RK4(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_RK3(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_RK2(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_RK1(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_PSSK5(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_PSSK4(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_PSSK3(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_PSSK2(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setSortedIndDefComDef_D1_PSSK1(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setSortedIndDefComDef_D1_BSSK(array, index)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndDefComDef_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedIndCusComDef

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedIndCusComDef_D0_SK5(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedIndCusComDef_D0_SK4(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedIndCusComDef_D0_SK3(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedIndCusComDef_D0_SK2(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedIndCusComDef_D0_SK1(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)                    :: array
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedIndCusComDef_D1_SK5(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedIndCusComDef_D1_SK4(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedIndCusComDef_D1_SK3(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedIndCusComDef_D1_SK2(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedIndCusComDef_D1_SK1(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSortedIndCusComDef_D1_IK5(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSortedIndCusComDef_D1_IK4(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSortedIndCusComDef_D1_IK3(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSortedIndCusComDef_D1_IK2(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSortedIndCusComDef_D1_IK1(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSortedIndCusComDef_D1_LK5(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSortedIndCusComDef_D1_LK4(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSortedIndCusComDef_D1_LK3(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSortedIndCusComDef_D1_LK2(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSortedIndCusComDef_D1_LK1(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSortedIndCusComDef_D1_CK5(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSortedIndCusComDef_D1_CK4(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSortedIndCusComDef_D1_CK3(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSortedIndCusComDef_D1_CK2(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSortedIndCusComDef_D1_CK1(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSortedIndCusComDef_D1_RK5(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSortedIndCusComDef_D1_RK4(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSortedIndCusComDef_D1_RK3(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSortedIndCusComDef_D1_RK2(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSortedIndCusComDef_D1_RK1(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setSortedIndCusComDef_D1_PSSK5(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedIndCusComDef_D1_PSSK4(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedIndCusComDef_D1_PSSK3(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedIndCusComDef_D1_PSSK2(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedIndCusComDef_D1_PSSK1(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))          , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSortedIndCusComDef_D1_BSSK(array, index, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedIndCusComDef_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)              , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(out)   , contiguous    :: index(:)
        procedure(logical(LK))                                      :: isSorted
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! setSortedArrCusComDef

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedArrCusComDef_D0_SK5(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComDef_D0_SK4(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComDef_D0_SK3(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComDef_D0_SK2(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComDef_D0_SK1(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedArrCusComDef_D1_SK5(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComDef_D1_SK4(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComDef_D1_SK3(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComDef_D1_SK2(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComDef_D1_SK1(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSortedArrCusComDef_D1_IK5(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSortedArrCusComDef_D1_IK4(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSortedArrCusComDef_D1_IK3(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSortedArrCusComDef_D1_IK2(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSortedArrCusComDef_D1_IK1(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSortedArrCusComDef_D1_LK5(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSortedArrCusComDef_D1_LK4(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSortedArrCusComDef_D1_LK3(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSortedArrCusComDef_D1_LK2(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSortedArrCusComDef_D1_LK1(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSortedArrCusComDef_D1_CK5(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSortedArrCusComDef_D1_CK4(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSortedArrCusComDef_D1_CK3(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSortedArrCusComDef_D1_CK2(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSortedArrCusComDef_D1_CK1(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSortedArrCusComDef_D1_RK5(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSortedArrCusComDef_D1_RK4(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSortedArrCusComDef_D1_RK3(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSortedArrCusComDef_D1_RK2(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSortedArrCusComDef_D1_RK1(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setSortedArrCusComDef_D1_PSSK5(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComDef_D1_PSSK4(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComDef_D1_PSSK3(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComDef_D1_PSSK2(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComDef_D1_PSSK1(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSortedArrCusComDef_D1_BSSK(array, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComDef_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrCusComQsorti

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedArrCusComQsorti_D0_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComQsorti_D0_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComQsorti_D0_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComQsorti_D0_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComQsorti_D0_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_IK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_IK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_IK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_IK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_IK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_LK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_LK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_LK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_LK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_LK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_CK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_CK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_CK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_CK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_CK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_RK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_RK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_RK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_RK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_RK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_PSSK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_PSSK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_PSSK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_PSSK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComQsorti_D1_PSSK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSortedArrCusComQsorti_D1_BSSK(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsorti_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrCusComQsortr

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D0_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D0_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D0_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D0_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D0_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_IK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_IK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_IK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_IK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_IK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_LK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_LK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_LK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_LK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_LK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_CK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_CK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_CK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_CK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_CK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_RK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_RK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_RK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_RK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_RK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_PSSK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_PSSK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_PSSK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_PSSK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    recursive module subroutine setSortedArrCusComQsortr_D1_PSSK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    recursive module subroutine setSortedArrCusComQsortr_D1_BSSK(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortr_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrCusComQsortrdp

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D0_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D0_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D0_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D0_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D0_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_IK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_IK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_IK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_IK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_IK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_LK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_LK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_LK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_LK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_LK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_CK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_CK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_CK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_CK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_CK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_RK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_RK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_RK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_RK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_RK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_PSSK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_PSSK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_PSSK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_PSSK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    recursive module subroutine setSortedArrCusComQsortrdp_D1_PSSK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    recursive module subroutine setSortedArrCusComQsortrdp_D1_BSSK(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComQsortrdp_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrCusComBubble

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedArrCusComBubble_D0_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComBubble_D0_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComBubble_D0_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComBubble_D0_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComBubble_D0_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedArrCusComBubble_D1_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComBubble_D1_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComBubble_D1_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComBubble_D1_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComBubble_D1_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSortedArrCusComBubble_D1_IK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSortedArrCusComBubble_D1_IK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSortedArrCusComBubble_D1_IK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSortedArrCusComBubble_D1_IK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSortedArrCusComBubble_D1_IK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSortedArrCusComBubble_D1_LK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSortedArrCusComBubble_D1_LK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSortedArrCusComBubble_D1_LK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSortedArrCusComBubble_D1_LK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSortedArrCusComBubble_D1_LK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSortedArrCusComBubble_D1_CK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSortedArrCusComBubble_D1_CK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSortedArrCusComBubble_D1_CK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSortedArrCusComBubble_D1_CK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSortedArrCusComBubble_D1_CK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSortedArrCusComBubble_D1_RK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSortedArrCusComBubble_D1_RK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSortedArrCusComBubble_D1_RK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSortedArrCusComBubble_D1_RK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSortedArrCusComBubble_D1_RK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setSortedArrCusComBubble_D1_PSSK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComBubble_D1_PSSK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComBubble_D1_PSSK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComBubble_D1_PSSK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComBubble_D1_PSSK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSortedArrCusComBubble_D1_BSSK(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComBubble_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(bubble_type)       , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrCusComHeapi

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedArrCusComHeapi_D0_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComHeapi_D0_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComHeapi_D0_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComHeapi_D0_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComHeapi_D0_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_IK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_IK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_IK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_IK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_IK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_LK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_LK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_LK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_LK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_LK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_CK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_CK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_CK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_CK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_CK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_RK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_RK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_RK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_RK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_RK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_PSSK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_PSSK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_PSSK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_PSSK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComHeapi_D1_PSSK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSortedArrCusComHeapi_D1_BSSK(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapi_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapi_type)        , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrCusComHeapr

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedArrCusComHeapr_D0_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComHeapr_D0_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComHeapr_D0_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComHeapr_D0_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComHeapr_D0_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_IK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_IK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_IK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_IK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_IK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_LK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_LK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_LK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_LK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_LK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_CK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_CK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_CK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_CK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_CK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_RK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_RK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_RK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_RK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_RK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_PSSK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_PSSK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_PSSK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_PSSK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComHeapr_D1_PSSK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSortedArrCusComHeapr_D1_BSSK(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComHeapr_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(heapr_type)        , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrCusComInsertionl

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedArrCusComInsertionl_D0_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComInsertionl_D0_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComInsertionl_D0_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComInsertionl_D0_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComInsertionl_D0_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_IK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_IK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_IK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_IK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_IK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_LK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_LK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_LK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_LK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_LK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_CK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_CK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_CK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_CK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_CK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_RK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_RK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_RK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_RK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_RK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_PSSK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_PSSK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_PSSK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_PSSK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComInsertionl_D1_PSSK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSortedArrCusComInsertionl_D1_BSSK(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionl_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrCusComInsertionb

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedArrCusComInsertionb_D0_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComInsertionb_D0_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComInsertionb_D0_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComInsertionb_D0_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComInsertionb_D0_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_IK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_IK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_IK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_IK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_IK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_LK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_LK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_LK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_LK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_LK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_CK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_CK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_CK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_CK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_CK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_RK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_RK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_RK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_RK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_RK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_PSSK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_PSSK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_PSSK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_PSSK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComInsertionb_D1_PSSK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSortedArrCusComInsertionb_D1_BSSK(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComInsertionb_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrCusComMerger

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D0_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D0_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D0_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D0_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D0_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_IK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_IK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_IK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_IK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_IK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_LK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_LK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_LK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_LK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_LK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_CK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_CK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_CK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_CK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_CK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_RK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_RK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_RK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_RK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_RK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_PSSK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_PSSK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_PSSK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_PSSK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    recursive module subroutine setSortedArrCusComMerger_D1_PSSK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    recursive module subroutine setSortedArrCusComMerger_D1_BSSK(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComMerger_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(merger_type)       , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrCusComSelection

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedArrCusComSelection_D0_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComSelection_D0_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComSelection_D0_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComSelection_D0_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComSelection_D0_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedArrCusComSelection_D1_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComSelection_D1_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComSelection_D1_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComSelection_D1_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComSelection_D1_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSortedArrCusComSelection_D1_IK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSortedArrCusComSelection_D1_IK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSortedArrCusComSelection_D1_IK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSortedArrCusComSelection_D1_IK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSortedArrCusComSelection_D1_IK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSortedArrCusComSelection_D1_LK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSortedArrCusComSelection_D1_LK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSortedArrCusComSelection_D1_LK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSortedArrCusComSelection_D1_LK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSortedArrCusComSelection_D1_LK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSortedArrCusComSelection_D1_CK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSortedArrCusComSelection_D1_CK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSortedArrCusComSelection_D1_CK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSortedArrCusComSelection_D1_CK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSortedArrCusComSelection_D1_CK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSortedArrCusComSelection_D1_RK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSortedArrCusComSelection_D1_RK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSortedArrCusComSelection_D1_RK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSortedArrCusComSelection_D1_RK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSortedArrCusComSelection_D1_RK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setSortedArrCusComSelection_D1_PSSK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComSelection_D1_PSSK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComSelection_D1_PSSK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComSelection_D1_PSSK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComSelection_D1_PSSK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSortedArrCusComSelection_D1_BSSK(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComSelection_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(selection_type)    , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrCusComShell

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedArrCusComShell_D0_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComShell_D0_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComShell_D0_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComShell_D0_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComShell_D0_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setSortedArrCusComShell_D1_SK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComShell_D1_SK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComShell_D1_SK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComShell_D1_SK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComShell_D1_SK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setSortedArrCusComShell_D1_IK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setSortedArrCusComShell_D1_IK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setSortedArrCusComShell_D1_IK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setSortedArrCusComShell_D1_IK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setSortedArrCusComShell_D1_IK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setSortedArrCusComShell_D1_LK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setSortedArrCusComShell_D1_LK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setSortedArrCusComShell_D1_LK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setSortedArrCusComShell_D1_LK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setSortedArrCusComShell_D1_LK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setSortedArrCusComShell_D1_CK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setSortedArrCusComShell_D1_CK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setSortedArrCusComShell_D1_CK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setSortedArrCusComShell_D1_CK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setSortedArrCusComShell_D1_CK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setSortedArrCusComShell_D1_RK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setSortedArrCusComShell_D1_RK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setSortedArrCusComShell_D1_RK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setSortedArrCusComShell_D1_RK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setSortedArrCusComShell_D1_RK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setSortedArrCusComShell_D1_PSSK5(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setSortedArrCusComShell_D1_PSSK4(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setSortedArrCusComShell_D1_PSSK3(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setSortedArrCusComShell_D1_PSSK2(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setSortedArrCusComShell_D1_PSSK1(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setSortedArrCusComShell_D1_BSSK(array, isSorted, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrCusComShell_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        procedure(logical(LK))                                  :: isSorted
        type(shell_type)        , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! setSortedArrDefComDef

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComDef_D0_SK5(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComDef_D0_SK4(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComDef_D0_SK3(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComDef_D0_SK2(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComDef_D0_SK1(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_SK5(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_SK4(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_SK3(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_SK2(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_SK1(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_IK5(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if IK4_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_IK4(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if IK3_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_IK3(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if IK2_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_IK2(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if IK1_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_IK1(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_LK5(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if LK4_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_LK4(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if LK3_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_LK3(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if LK2_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_LK2(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if LK1_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_LK1(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_CK5(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if CK4_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_CK4(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if CK3_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_CK3(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if CK2_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_CK2(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if CK1_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_CK1(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_RK5(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if RK4_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_RK4(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if RK3_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_RK3(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if RK2_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_RK2(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if RK1_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_RK1(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_PSSK5(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_PSSK4(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_PSSK3(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_PSSK2(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComDef_D1_PSSK1(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure module subroutine setSortedArrDefComDef_D1_BSSK(array)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComDef_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrDefComQsorti

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D0_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D0_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D0_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D0_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D0_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_IK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_IK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_IK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_IK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_IK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_LK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_LK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_LK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_LK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_LK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_CK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_CK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_CK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_CK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_CK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_RK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_RK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_RK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_RK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_RK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_PSSK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_PSSK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_PSSK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_PSSK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComQsorti_D1_PSSK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure module subroutine setSortedArrDefComQsorti_D1_BSSK(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsorti_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        type(qsorti_type)       , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrDefComQsortr

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D0_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D0_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D0_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D0_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D0_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_IK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_IK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_IK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_IK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_IK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_LK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_LK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_LK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_LK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_LK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_CK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_CK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_CK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_CK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_CK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_RK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_RK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_RK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_RK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_RK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_PSSK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_PSSK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_PSSK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_PSSK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortr_D1_PSSK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure recursive module subroutine setSortedArrDefComQsortr_D1_BSSK(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortr_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        type(qsortr_type)       , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrDefComQsortrdp

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D0_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D0_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D0_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D0_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D0_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_IK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_IK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_IK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_IK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_IK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_LK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_LK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_LK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_LK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_LK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_CK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_CK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_CK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_CK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_CK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_RK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_RK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_RK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_RK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_RK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_PSSK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_PSSK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_PSSK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_PSSK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_PSSK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure recursive module subroutine setSortedArrDefComQsortrdp_D1_BSSK(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComQsortrdp_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        type(qsortrdp_type)     , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrDefComBubble

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComBubble_D0_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComBubble_D0_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComBubble_D0_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComBubble_D0_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComBubble_D0_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_IK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_IK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_IK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_IK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_IK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_LK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_LK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_LK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_LK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_LK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_CK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_CK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_CK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_CK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_CK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_RK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_RK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_RK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_RK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_RK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_PSSK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_PSSK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_PSSK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_PSSK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComBubble_D1_PSSK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure module subroutine setSortedArrDefComBubble_D1_BSSK(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComBubble_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        type(bubble_type)       , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrDefComHeapi

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D0_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D0_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D0_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D0_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D0_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_IK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_IK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_IK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_IK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_IK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_LK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_LK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_LK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_LK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_LK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_CK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_CK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_CK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_CK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_CK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_RK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_RK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_RK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_RK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_RK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_PSSK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_PSSK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_PSSK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_PSSK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComHeapi_D1_PSSK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure module subroutine setSortedArrDefComHeapi_D1_BSSK(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapi_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        type(heapi_type)        , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrDefComHeapr

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D0_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D0_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D0_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D0_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D0_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_IK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_IK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_IK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_IK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_IK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_LK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_LK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_LK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_LK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_LK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_CK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_CK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_CK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_CK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_CK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_RK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_RK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_RK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_RK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_RK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_PSSK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_PSSK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_PSSK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_PSSK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComHeapr_D1_PSSK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure module subroutine setSortedArrDefComHeapr_D1_BSSK(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComHeapr_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        type(heapr_type)        , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrDefComInsertionl

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D0_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D0_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D0_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D0_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D0_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_IK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_IK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_IK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_IK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_IK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_LK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_LK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_LK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_LK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_LK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_CK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_CK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_CK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_CK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_CK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_RK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_RK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_RK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_RK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_RK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_PSSK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_PSSK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_PSSK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_PSSK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComInsertionl_D1_PSSK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure module subroutine setSortedArrDefComInsertionl_D1_BSSK(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionl_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        type(insertionl_type)   , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrDefComInsertionb

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D0_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D0_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D0_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D0_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D0_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_IK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_IK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_IK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_IK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_IK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_LK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_LK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_LK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_LK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_LK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_CK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_CK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_CK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_CK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_CK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_RK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_RK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_RK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_RK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_RK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_PSSK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_PSSK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_PSSK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_PSSK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComInsertionb_D1_PSSK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure module subroutine setSortedArrDefComInsertionb_D1_BSSK(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComInsertionb_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        type(insertionb_type)   , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrDefComMerger

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D0_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D0_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D0_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D0_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D0_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_IK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_IK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_IK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_IK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_IK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_LK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_LK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_LK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_LK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_LK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_CK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_CK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_CK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_CK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_CK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_RK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_RK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_RK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_RK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_RK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_PSSK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_PSSK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_PSSK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_PSSK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    PURE recursive module subroutine setSortedArrDefComMerger_D1_PSSK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE recursive module subroutine setSortedArrDefComMerger_D1_BSSK(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComMerger_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        type(merger_type)       , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrDefComSelection

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComSelection_D0_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComSelection_D0_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComSelection_D0_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComSelection_D0_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComSelection_D0_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_IK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_IK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_IK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_IK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_IK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_LK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_LK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_LK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_LK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_LK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_CK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_CK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_CK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_CK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_CK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_RK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_RK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_RK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_RK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_RK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_PSSK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_PSSK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_PSSK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_PSSK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComSelection_D1_PSSK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure module subroutine setSortedArrDefComSelection_D1_BSSK(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComSelection_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        type(selection_type)    , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setSortedArrDefComShell

    interface setSorted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComShell_D0_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout)                 :: array
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComShell_D0_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout)                 :: array
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComShell_D0_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout)                 :: array
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComShell_D0_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout)                 :: array
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComShell_D0_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout)                 :: array
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_SK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_SK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_SK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_SK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_SK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_IK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK4_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_IK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK3_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_IK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK2_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_IK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if IK1_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_IK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_LK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK4_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_LK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK3_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_LK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK2_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_LK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if LK1_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_LK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_CK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK4_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_CK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK3_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_CK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK2_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_CK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if CK1_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_CK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_RK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK4_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_RK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK3_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_RK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK2_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_RK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if RK1_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_RK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_PSSK5(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK4_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_PSSK4(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK3_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_PSSK3(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK2_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_PSSK2(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#if SK1_ENABLED
    pure module subroutine setSortedArrDefComShell_D1_PSSK1(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure module subroutine setSortedArrDefComShell_D1_BSSK(array, method)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setSortedArrDefComShell_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(inout) , contiguous    :: array(:)
        type(shell_type)        , intent(in)                    :: method
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arraySort