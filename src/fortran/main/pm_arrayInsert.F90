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
!>  This module contains procedures and generic interfaces for inserting an insertion into the specified locations of an input arrays of various types.
!>
!>  \benchmarks
!>
!>  \benchmark{scalarInsertion_vs_vectorInsertion, The runtime performance of [setInserted](@ref pm_arrayInsert::setInserted) for scalar vs. vector input `insertion` argument.}
!>  \include{lineno} benchmark/pm_arrayInsert/scalarInsertion_vs_vectorInsertion/main.F90
!>  \compilefb{scalarInsertion_vs_vectorInsertion}
!>  \postprocb{scalarInsertion_vs_vectorInsertion}
!>  \include{lineno} benchmark/pm_arrayInsert/scalarInsertion_vs_vectorInsertion/main.py
!>  \visb{scalarInsertion_vs_vectorInsertion}
!>  \image html benchmark/pm_arrayInsert/scalarInsertion_vs_vectorInsertion/benchmark.scalarInsertion_vs_vectorInsertion.runtime.png width=1000
!>  \image html benchmark/pm_arrayInsert/scalarInsertion_vs_vectorInsertion/benchmark.scalarInsertion_vs_vectorInsertion.runtime.ratio.png width=1000
!>  \moralb{scalarInsertion_vs_vectorInsertion}
!>      -#  The procedures under the generic interface [setInserted](@ref pm_arrayInsert::setInserted) take both scalar and vector `insertion` arguments.<br>
!>          As evidenced by the above benchmark, when the input `insertion` is vector of length `1`, it is much faster, by **4X** or more,
!>          to pass `insertion` as a scalar instead of a whole array of length `1`.<br>
!>          This benchmark represents the worst-case scenario.<br>
!>          Note that this benchmark is likely irrelevant to inserting substrings to Fortran strings.<br>
!>
!>  \benchmark{getInserted_vs_setInserted, The runtime performance of [getInserted](@ref pm_arrayInsert::getInserted) vs. [setInserted](@ref pm_arrayInsert::setInserted)}
!>  \include{lineno} benchmark/pm_arrayInsert/getInserted_vs_setInserted/main.F90
!>  \compilefb{getInserted_vs_setInserted}
!>  \postprocb{getInserted_vs_setInserted}
!>  \include{lineno} benchmark/pm_arrayInsert/getInserted_vs_setInserted/main.py
!>  \visb{getInserted_vs_setInserted}
!>  \image html benchmark/pm_arrayInsert/getInserted_vs_setInserted/benchmark.getInserted_vs_setInserted.runtime.png width=1000
!>  \image html benchmark/pm_arrayInsert/getInserted_vs_setInserted/benchmark.getInserted_vs_setInserted.runtime.ratio.png width=1000
!>  \moralb{getInserted_vs_setInserted}
!>      -#  The procedures under the generic interface [getInserted](@ref pm_arrayInsert::getInserted) are functions while
!>          the procedures under the generic interface [setInserted](@ref pm_arrayInsert::setInserted) are subroutines.<br>
!>          From the benchmark results, it appears that the functional interface performs slightly less efficiently than the subroutine interface.<br>
!>          The sole apparent reason for this performance loss seems to be the extra copy of the result to the allocatable `arrayNew` on return from the function.<br>
!>          Note that this benchmark does not even include the cost of repeated reallcations, that is, the allocation of `arrayNew` happen only once in all tests.<br>
!>      -#  Note that this benchmark considers the worst-case scenario where `insertion` must be inserted at all positions of the input `array`.
!>
!>  \test
!>  [test_pm_arrayInsert](@ref test_pm_arrayInsert)
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:20 AM, October 13, 2021, Dallas, TX

module pm_arrayInsert

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_arrayInsert"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a new array containing the original array within which the input `insertion` has been inserted at
    !>  the specified indices `index` of the original array.
    !>
    !>  \param[in]  array       :   The output `contiguous` array of shape `(:)` of either <br>
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL, or <br>
    !>                                  <li>    type `logical` of kind \LKALL, or <br>
    !>                                  <li>    type `integer` of kind \IKALL, or <br>
    !>                                  <li>    type `complex` of kind \CKALL, or <br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              or,
    !>                              <ul>
    !>                                  <li>    a scalar assumed-length `character` of kind \SKALL, <br>
    !>                              </ul>
    !>                              containing the sequence of values within which the input `insertion` must be inserted at the specified indices.
    !>  \param[in]  insertion   :   The input scalar or `contiguous` array of shape `(:)` of the same type and kind as the input `array`
    !>                              containing the insertion that must be inserted at the specified indices of the input `array`.
    !>  \param[in]  index       :   The input `contiguous` array of shape `(:)` of type `integer` of default kind \IK,
    !>                              containing the indices of the input `array` where `insertion` must be inserted to construct `arrayNew`.<br>
    !>                              <ul>
    !>                                  <li>    All elements of `index` must have values between `1` and the length of `array` plus `1`.<br>
    !>                                  <li>    Any insertion of `index` that is negatively valued will be counted from end of the input `array`.<br>
    !>                                  <li>    For example, `index = [2,-1]` requests inserting `insertion` at `array(2)` and `array(lenArray)`.<br>
    !>                                  <li>    To append `insertion` to `array`, specify an index value of `lenArray + 1`, for example, `index = [lenArray + 1]`.<br>
    !>                                  <li>    If any value appears repeatedly and sequentially within `index`, then the corresponding number of instances of `insertion`
    !>                                          will be inserted sequentially at the specified position within `array`.
    !>                              </ul>
    !>  \param[in]  positive    :   The input `logical` of default kind \LK indicating whether the elements of `index` are all positive (indicating counts from the beginning of `array`).<br>
    !>                              Setting `positive = .true.` will lead to a slightly better runtime performance of the algorithm since a conversion of potentially-negative index
    !>                              values to the corresponding positive values (counting from the beginning of `array` will be avoided.<br>
    !>                              (**optional**, default = `.false.`)
    !>  \param[in]  sorted      :   The input `logical` of default kind \LK indicating whether the insertions of the specified input `index` are all in ascending-order.<br>
    !>                              This includes the negative values of `index` **after** they are converted to the corresponding **positive** indices from the beginning of the input `array`.<br>
    !>                              Setting `sorted = .true.` will lead to a better runtime performance of the algorithm since a call to the
    !>                              [setSorted](@ref pm_arraySort::setSorted) to sort the index values in ascending order will be avoided.<br>
    !>                              However, the onus will be on the user to guarantee the ascending order of the elements of the input index.<br>
    !>                              (**optional**, default = `.false.`)
    !>
    !>  \return
    !>  `arrayNew`              :   The output array of the same type, kind, and rank as the input `array`, and containing the original array
    !>                              within which the input `insertion` has been inserted at the requested indices of the `array`.<br>
    !>                              The size of `arrayNew` will be,
    !>                              <ul>
    !>                                  <li>    `size(array,kind=IK) + size(index,kind=IK) * size(insertion)` if `array` is a non-scalar-character and `insertion` is a vector,<br>
    !>                                  <li>    `size(array,kind=IK) + size(index,kind=IK) * 1` if `array` is a non-scalar-character and `insertion` is a scalar,<br>
    !>                                  <li>    `len(array) + size(index,kind=IK) * len(insertion)` if both `array` and `insertion` are scalar characters.<br>
    !>                              </ul>
    !>
    !>  \interface{getInserted}
    !>  \code{.F90}
    !>
    !>      use pm_arrayInsert, only: LK, getInserted
    !>
    !>      arrayNew = getInserted(array, insertion, index)
    !>      arrayNew = getInserted(array, insertion, index, positive)
    !>      arrayNew = getInserted(array, insertion, index, sorted = sorted)
    !>      arrayNew = getInserted(array, insertion, index, positive, sorted)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  All elements of `index` must have values between `1` and the length of `array` plus `1`.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [setInserted](@ref pm_arrayInsert::setInserted)<br>
    !>  [getRemoved](@ref pm_arrayRemove::getRemoved)<br>
    !>  [getReplaced](@ref pm_arrayReplace::getReplaced)<br>
    !>  [setReplaced](@ref pm_arrayReplace::setReplaced)<br>
    !>  [setSplit](@ref pm_arraySplit::setSplit)<br>
    !>
    !>  \example{getInserted}
    !>  \include{lineno} example/pm_arrayInsert/getInserted/main.F90
    !>  \compilef{getInserted}
    !>  \output{getInserted}
    !>  \include{lineno} example/pm_arrayInsert/getInserted/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayInsert](@ref test_pm_arrayInsert)
    !>
    !>  \todo
    !>  \plow This generic interface can be extended to 2D input objects.
    !>
    !>  \todo
    !>  \pmed A benchmark comparing the performance of [getInserted](@ref pm_arrayInsert::getInserted) with and without `positive, sorted` should be added.
    !>
    !>  \final{getInserted}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getInserted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getInserted_D0_D0_SK5(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                                :: array
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(len(array,IK)+size(index,kind=IK)*len(insertion,IK),SKG)  :: arrayNew
    end function
#endif

#if SK4_ENABLED
    PURE module function getInserted_D0_D0_SK4(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                                :: array
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(len(array,IK)+size(index,kind=IK)*len(insertion,IK),SKG)  :: arrayNew
    end function
#endif

#if SK3_ENABLED
    PURE module function getInserted_D0_D0_SK3(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                                :: array
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(len(array,IK)+size(index,kind=IK)*len(insertion,IK),SKG)  :: arrayNew
    end function
#endif

#if SK2_ENABLED
    PURE module function getInserted_D0_D0_SK2(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                                :: array
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(len(array,IK)+size(index,kind=IK)*len(insertion,IK),SKG)  :: arrayNew
    end function
#endif

#if SK1_ENABLED
    PURE module function getInserted_D0_D0_SK1(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                                :: array
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(len(array,IK)+size(index,kind=IK)*len(insertion,IK),SKG)  :: arrayNew
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getInserted_D1_D0_SK5(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(len(array,IK),SKG)                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if SK4_ENABLED
    PURE module function getInserted_D1_D0_SK4(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(len(array,IK),SKG)                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if SK3_ENABLED
    PURE module function getInserted_D1_D0_SK3(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(len(array,IK),SKG)                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if SK2_ENABLED
    PURE module function getInserted_D1_D0_SK2(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(len(array,IK),SKG)                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if SK1_ENABLED
    PURE module function getInserted_D1_D0_SK1(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(len(array,IK),SKG)                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getInserted_D1_D0_IK5(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if IK4_ENABLED
    PURE module function getInserted_D1_D0_IK4(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if IK3_ENABLED
    PURE module function getInserted_D1_D0_IK3(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if IK2_ENABLED
    PURE module function getInserted_D1_D0_IK2(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if IK1_ENABLED
    PURE module function getInserted_D1_D0_IK1(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getInserted_D1_D0_LK5(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if LK4_ENABLED
    PURE module function getInserted_D1_D0_LK4(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if LK3_ENABLED
    PURE module function getInserted_D1_D0_LK3(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if LK2_ENABLED
    PURE module function getInserted_D1_D0_LK2(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if LK1_ENABLED
    PURE module function getInserted_D1_D0_LK1(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getInserted_D1_D0_CK5(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getInserted_D1_D0_CK4(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getInserted_D1_D0_CK3(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getInserted_D1_D0_CK2(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getInserted_D1_D0_CK1(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getInserted_D1_D0_RK5(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)                                                           :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getInserted_D1_D0_RK4(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)                                                           :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getInserted_D1_D0_RK3(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)                                                           :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getInserted_D1_D0_RK2(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)                                                           :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getInserted_D1_D0_RK1(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)                                                           :: arrayNew(size(array,kind=IK)+size(index,kind=IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getInserted_D1_D1_SK5(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(len(array,IK),SKG)                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion,kind=IK))
    end function
#endif

#if SK4_ENABLED
    PURE module function getInserted_D1_D1_SK4(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(len(array,IK),SKG)                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion,kind=IK))
    end function
#endif

#if SK3_ENABLED
    PURE module function getInserted_D1_D1_SK3(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(len(array,IK),SKG)                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion,kind=IK))
    end function
#endif

#if SK2_ENABLED
    PURE module function getInserted_D1_D1_SK2(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(len(array,IK),SKG)                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion,kind=IK))
    end function
#endif

#if SK1_ENABLED
    PURE module function getInserted_D1_D1_SK1(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(len(array,IK),SKG)                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion,kind=IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getInserted_D1_D1_IK5(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

#if IK4_ENABLED
    PURE module function getInserted_D1_D1_IK4(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

#if IK3_ENABLED
    PURE module function getInserted_D1_D1_IK3(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

#if IK2_ENABLED
    PURE module function getInserted_D1_D1_IK2(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

#if IK1_ENABLED
    PURE module function getInserted_D1_D1_IK1(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getInserted_D1_D1_LK5(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

#if LK4_ENABLED
    PURE module function getInserted_D1_D1_LK4(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

#if LK3_ENABLED
    PURE module function getInserted_D1_D1_LK3(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

#if LK2_ENABLED
    PURE module function getInserted_D1_D1_LK2(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

#if LK1_ENABLED
    PURE module function getInserted_D1_D1_LK1(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getInserted_D1_D1_CK5(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

#if CK4_ENABLED
    PURE module function getInserted_D1_D1_CK4(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

#if CK3_ENABLED
    PURE module function getInserted_D1_D1_CK3(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

#if CK2_ENABLED
    PURE module function getInserted_D1_D1_CK2(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

#if CK1_ENABLED
    PURE module function getInserted_D1_D1_CK1(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)                                                        :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getInserted_D1_D1_RK5(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)                                                           :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

#if RK4_ENABLED
    PURE module function getInserted_D1_D1_RK4(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)                                                           :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

#if RK3_ENABLED
    PURE module function getInserted_D1_D1_RK3(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)                                                           :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

#if RK2_ENABLED
    PURE module function getInserted_D1_D1_RK2(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)                                                           :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

#if RK1_ENABLED
    PURE module function getInserted_D1_D1_RK1(array, insertion, index, positive, sorted) result(arrayNew)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInserted_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)                                                           :: arrayNew(size(array,kind=IK)+size(index,kind=IK)*size(insertion))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a new array `arrayNew` containing the original `array` within which the input `insertion` has been inserted at
    !>  the specified indices `index` of the original array.
    !>
    !>  \param[out] arrayNew    :   The output `contiguous` array of shape `(:)` of either <br>
    !>                              <ul>
    !>                                  <li>    type `character` of kind \SKALL, or <br>
    !>                                  <li>    type `logical` of kind \LKALL, or <br>
    !>                                  <li>    type `integer` of kind \IKALL, or <br>
    !>                                  <li>    type `complex` of kind \CKALL, or <br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              or,
    !>                              <ul>
    !>                                  <li>    a scalar assumed-length `character` of kind \SKALL, <br>
    !>                              </ul>
    !>                              containing the original array within which the input `insertion` has been inserted at the requested indices of the `array`.<br>
    !>                              The size of `arrayNew` must be,<br>
    !>                              <ul>
    !>                                  <li>    `size(array,kind=IK) + size(index,kind=IK) * size(insertion)` if `array` is a non-scalar-character and `insertion` is a vector,<br>
    !>                                  <li>    `size(array,kind=IK) + size(index,kind=IK) * 1` if `array` is a non-scalar-character and `insertion` is a scalar,<br>
    !>                                  <li>    `len(array) + size(index,kind=IK) * len(insertion)` if both `array` and `insertion` are scalar characters.<br>
    !>                              </ul>
    !>  \param[in]  array       :   The input `contiguous` array of shape `(:)` of the same type and kind as `arrayNew` containing the sequence of values
    !>                              within which the input `insertion` must be inserted at the specified indices.
    !>  \param[in]  insertion   :   The input scalar or `contiguous` array of shape `(:)` of the same type and kind as the input `array`
    !>                              containing the insertion that must be inserted at the specified indices of the input `array`.
    !>  \param[in]  index       :   The input `contiguous` array of shape `(:)` of type `integer` of default kind \IK,
    !>                              containing the indices of the input `array` where `insertion` must be inserted to construct `arrayNew`.<br>
    !>                              <ul>
    !>                                  <li>    All elements of `index` must have values between `1` and the length of `array` plus `1`.<br>
    !>                                  <li>    Any insertion of `index` that is negatively valued will be counted from end of the input `array`.<br>
    !>                                  <li>    For example, `index = [2,-1]` requests inserting `insertion` at `array(2)` and `array(lenArray)`.<br>
    !>                                  <li>    To append `insertion` to `array`, specify an index value of `lenArray + 1`, for example, `index = [lenArray + 1]`.<br>
    !>                                  <li>    If any value appears repeatedly and sequentially within `index`, then the corresponding number of instances of `insertion` will
    !>                              </ul>
    !>                              be inserted sequentially at the specified position within `array`.
    !>  \param[in]  positive    :   The input `logical` of default kind \LK indicating whether the elements of `index` are all positive (indicating counts from the beginning of `array`).<br>
    !>                              Setting `positive = .true.` will lead to a slightly better runtime performance of the algorithm since a conversion of potentially-negative index
    !>                              values to the corresponding positive values (counting from the beginning of `array` will be avoided.<br>
    !>                              (**optional**, default = `.false.`)
    !>  \param[in]  sorted      :   The input `logical` of default kind \LK indicating whether the insertions of the specified input `index` are all in ascending-order.<br>
    !>                              This includes the negative values of `index` **after** they are converted to the corresponding **positive** indices from the beginning of the input `array`.<br>
    !>                              Setting `sorted = .true.` will lead to a better runtime performance of the algorithm since a call to the
    !>                              [setSorted](@ref pm_arraySort::setSorted) to sort the index values in ascending order will be avoided.<br>
    !>                              However, the onus will be on the user to guarantee the ascending order of the elements of the input index.<br>
    !>                              (**optional**, default = `.false.`)
    !>
    !>  \interface{setInserted}
    !>  \code{.F90}
    !>
    !>      use pm_arrayInsert, only: setInserted
    !>
    !>      call setInserted(arrayNew, array, insertion, index)
    !>      call setInserted(arrayNew, array, insertion, index, positive)
    !>      call setInserted(arrayNew, array, insertion, index, sorted = sorted)
    !>      call setInserted(arrayNew, array, insertion, index, positive, sorted)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  All elements of `index` must have values between `1` and the length of `array` plus `1`.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getRemoved](@ref pm_arrayRemove::getRemoved)<br>
    !>  [getInserted](@ref pm_arrayInsert::getInserted)<br>
    !>  [getReplaced](@ref pm_arrayReplace::getReplaced)<br>
    !>  [setReplaced](@ref pm_arrayReplace::setReplaced)<br>
    !>  [setSplit](@ref pm_arraySplit::setSplit)<br>
    !>
    !>  \example{setInserted}
    !>  \include{lineno} example/pm_arrayInsert/setInserted/main.F90
    !>  \compilef{setInserted}
    !>  \output{setInserted}
    !>  \include{lineno} example/pm_arrayInsert/setInserted/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayInsert](@ref test_pm_arrayInsert)
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to 2D input objects.<br>
    !>
    !>  \todo
    !>  \pvlow This generic interface can be extended to accept an `intent(inout), allocatable :: array` instead of `ArraNew`
    !>  to simplify in-place insertion. However, the potential resulting code bloat outweigh the slightly improved calling syntax benefits.<br>
    !>
    !>  \todo
    !>  \pmed A benchmark comparing the performance of [setInserted](@ref pm_arrayInsert::setInserted) with and without `positive, sorted` should be added.<br>
    !>
    !>  \final{setInserted}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setInserted

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setInserted_D0_D0_SK5(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D0_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                                :: array
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(*,SKG)        , intent(out)                               :: arrayNew
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setInserted_D0_D0_SK4(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D0_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                                :: array
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(*,SKG)        , intent(out)                               :: arrayNew
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setInserted_D0_D0_SK3(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D0_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                                :: array
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(*,SKG)        , intent(out)                               :: arrayNew
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setInserted_D0_D0_SK2(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D0_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                                :: array
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(*,SKG)        , intent(out)                               :: arrayNew
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setInserted_D0_D0_SK1(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D0_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                                :: array
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(*,SKG)        , intent(out)                               :: arrayNew
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setInserted_D1_D0_SK5(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(*,SKG)        , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setInserted_D1_D0_SK4(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(*,SKG)        , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setInserted_D1_D0_SK3(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(*,SKG)        , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setInserted_D1_D0_SK2(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(*,SKG)        , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setInserted_D1_D0_SK1(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(*,SKG)        , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setInserted_D1_D0_IK5(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setInserted_D1_D0_IK4(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setInserted_D1_D0_IK3(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setInserted_D1_D0_IK2(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setInserted_D1_D0_IK1(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setInserted_D1_D0_LK5(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setInserted_D1_D0_LK4(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setInserted_D1_D0_LK3(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setInserted_D1_D0_LK2(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setInserted_D1_D0_LK1(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setInserted_D1_D0_CK5(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setInserted_D1_D0_CK4(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setInserted_D1_D0_CK3(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setInserted_D1_D0_CK2(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setInserted_D1_D0_CK1(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setInserted_D1_D0_RK5(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)               , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setInserted_D1_D0_RK4(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)               , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setInserted_D1_D0_RK3(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)               , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setInserted_D1_D0_RK2(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)               , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setInserted_D1_D0_RK1(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D0_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)                                :: insertion
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)               , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setInserted_D1_D1_SK5(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(*,SKG)        , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setInserted_D1_D1_SK4(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(*,SKG)        , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setInserted_D1_D1_SK3(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(*,SKG)        , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setInserted_D1_D1_SK2(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(*,SKG)        , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setInserted_D1_D1_SK1(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous                :: array(:)
        character(*,SKG)        , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        character(*,SKG)        , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setInserted_D1_D1_IK5(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setInserted_D1_D1_IK4(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setInserted_D1_D1_IK3(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setInserted_D1_D1_IK2(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setInserted_D1_D1_IK1(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous                :: array(:)
        integer(IKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        integer(IKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setInserted_D1_D1_LK5(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setInserted_D1_D1_LK4(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setInserted_D1_D1_LK3(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setInserted_D1_D1_LK2(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setInserted_D1_D1_LK1(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous                :: array(:)
        logical(LKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        logical(LKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setInserted_D1_D1_CK5(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setInserted_D1_D1_CK4(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setInserted_D1_D1_CK3(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setInserted_D1_D1_CK2(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setInserted_D1_D1_CK1(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous                :: array(:)
        complex(CKG)            , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        complex(CKG)            , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setInserted_D1_D1_RK5(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)               , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setInserted_D1_D1_RK4(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)               , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setInserted_D1_D1_RK3(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)               , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setInserted_D1_D1_RK2(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)               , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setInserted_D1_D1_RK1(arrayNew, array, insertion, index, positive, sorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInserted_D1_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous                :: array(:)
        real(RKG)               , intent(in)    , contiguous                :: insertion(:)
        integer(IK)             , intent(in)    , contiguous                :: index(:)
        logical(LK)             , intent(in)    , optional                  :: positive, sorted
        real(RKG)               , intent(out)   , contiguous                :: arrayNew(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayInsert ! LCOV_EXCL_LINE