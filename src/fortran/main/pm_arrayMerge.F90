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
!>  This module contains procedures and generic interfaces for sorting and merging two **previously-sorted** arrays.
!>
!>  \test
!>  [test_pm_arrayMerge](@ref test_pm_arrayMerge)
!>
!>  \todo
!>  \plow
!>  The examples of this module require an update.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayMerge

    use pm_kind, only: IK, LK, RK, SK

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_arrayMerge"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate an ascending-sorted merger of the two ascending-sorted input arrays.
    !>
    !>  \param[in]  sortedArray1    :   The input scalar of either<br>
    !>                                  <ul>
    !>                                      <li>    type `character` of kind \SKALL,<br>
    !>                                  </ul>
    !>                                  or the input `contiguous` array of rank `1` of either<br>
    !>                                  <ul>
    !>                                      <li>    type `character` of kind \SKALL or<br>
    !>                                      <li>    type `logical` of kind \LKALL or<br>
    !>                                      <li>    type `integer` of kind \IKALL or<br>
    !>                                      <li>    type `complex` of kind \CKALL or<br>
    !>                                      <li>    type `real` of kind \RKALL or<br>
    !>                                      <li>    type [css_pdt](@ref pm_container::css_pdt) of kind \SKALL or<br>
    !>                                      <li>    type [css_type](@ref pm_container::css_type) of default kind \SK,<br>
    !>                                  </ul>
    !>                                  whose elements **are already sorted** in ascending order.<br>
    !>                                  When the type of the input argument is `complex`, the array must be sorted only based on the real component of the values.<br>
    !>  \param[in]  sortedArray2    :   The input `contiguous` array of the same type, kind, and rank as the input `sortedArray1` whose elements **are already sorted** in ascending order.<br>
    !>                                  When the type of the input argument is `complex`, the array must be sorted only based on the real component of the values.<br>
    !>  \param          isSorted    :   The `external` user-specified function that takes two input **scalar** arguments of the same type and kind as the input `Array`.<br>
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
    !>                                  where `TYPE(KIND)` is the same as the type and kind of the input argument `Array`, which can be one of the following.
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
    !>
    !>  \return
    !>  `mergedSortedArray`         :   The output object of the same type, kind, and rank as the input `sortedArray1`
    !>                                  whose length is the sum of the lengths of the input `sortedArray1` and `sortedArray2` and
    !>                                  whose elements are the combined elements of `sortedArray1` and `sortedArray2` in ascending order.
    !>
    !>  \interface{getMerged}
    !>  \code{.F90}
    !>
    !>      use pm_arrayMerge, only: getMerged
    !>
    !>      mergedSortedArray(1 : size(sortedArray1) + size(sortedArray2)) = getMerged(sortedArray1, sortedArray2)
    !>      mergedSortedArray(1 : size(sortedArray1) + size(sortedArray2)) = getMerged(sortedArray1, sortedArray2, isSorted)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The input arguments `sortedArray1` and `sortedArray2` must be sorted in ascending order when `isSorted()` is missing, or properly sorted when `isSorted()` is present.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  This generic interface is a functional version of the subroutine implementation [pm_arrayMerge::setMerged](@ref pm_arrayMerge::setMerged).
    !>
    !>  \see
    !>  [pm_arrayMerge::setMerged](@ref pm_arrayMerge::setMerged)<br>
    !>
    !>  \example{getMerged}
    !>  \include{lineno} example/pm_arrayMerge/getMerged/main.F90
    !>  \compilef{getMerged}
    !>  \output{getMerged}
    !>  \include{lineno} example/pm_arrayMerge/getMerged/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayMerge](@ref test_pm_arrayMerge)
    !>
    !>  \final{getMerged}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! DefCom

    interface getMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getMergedDefCom_D0_SK5(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(len(sortedArray1,IK) + len(sortedArray2,IK))  :: mergedSortedArray
    end function
#endif

#if SK4_ENABLED
    PURE module function getMergedDefCom_D0_SK4(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(len(sortedArray1,IK) + len(sortedArray2,IK))  :: mergedSortedArray
    end function
#endif

#if SK3_ENABLED
    PURE module function getMergedDefCom_D0_SK3(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(len(sortedArray1,IK) + len(sortedArray2,IK))  :: mergedSortedArray
    end function
#endif

#if SK2_ENABLED
    PURE module function getMergedDefCom_D0_SK2(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(len(sortedArray1,IK) + len(sortedArray2,IK))  :: mergedSortedArray
    end function
#endif

#if SK1_ENABLED
    PURE module function getMergedDefCom_D0_SK1(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(len(sortedArray1,IK) + len(sortedArray2,IK))  :: mergedSortedArray
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getMergedDefCom_D1_SK5(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(len(sortedArray1,IK))                         :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if SK4_ENABLED
    PURE module function getMergedDefCom_D1_SK4(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(len(sortedArray1,IK))                         :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if SK3_ENABLED
    PURE module function getMergedDefCom_D1_SK3(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(len(sortedArray1,IK))                         :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if SK2_ENABLED
    PURE module function getMergedDefCom_D1_SK2(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(len(sortedArray1,IK))                         :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if SK1_ENABLED
    PURE module function getMergedDefCom_D1_SK1(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(len(sortedArray1,IK))                         :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getMergedDefCom_D1_IK5(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if IK4_ENABLED
    PURE module function getMergedDefCom_D1_IK4(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if IK3_ENABLED
    PURE module function getMergedDefCom_D1_IK3(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if IK2_ENABLED
    PURE module function getMergedDefCom_D1_IK2(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if IK1_ENABLED
    PURE module function getMergedDefCom_D1_IK1(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getMergedDefCom_D1_LK5(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if LK4_ENABLED
    PURE module function getMergedDefCom_D1_LK4(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if LK3_ENABLED
    PURE module function getMergedDefCom_D1_LK3(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if LK2_ENABLED
    PURE module function getMergedDefCom_D1_LK2(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if LK1_ENABLED
    PURE module function getMergedDefCom_D1_LK1(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMergedDefCom_D1_CK5(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if CK4_ENABLED
    PURE module function getMergedDefCom_D1_CK4(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if CK3_ENABLED
    PURE module function getMergedDefCom_D1_CK3(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if CK2_ENABLED
    PURE module function getMergedDefCom_D1_CK2(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if CK1_ENABLED
    PURE module function getMergedDefCom_D1_CK1(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMergedDefCom_D1_RK5(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)                                               :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if RK4_ENABLED
    PURE module function getMergedDefCom_D1_RK4(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)                                               :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if RK3_ENABLED
    PURE module function getMergedDefCom_D1_RK3(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)                                               :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if RK2_ENABLED
    PURE module function getMergedDefCom_D1_RK2(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)                                               :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if RK1_ENABLED
    PURE module function getMergedDefCom_D1_RK1(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)                                               :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module function getMergedDefCom_D1_PSSK5(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))                                      :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if SK4_ENABLED
    PURE module function getMergedDefCom_D1_PSSK4(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))                                      :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if SK3_ENABLED
    PURE module function getMergedDefCom_D1_PSSK3(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))                                      :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if SK2_ENABLED
    PURE module function getMergedDefCom_D1_PSSK2(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))                                      :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#if SK1_ENABLED
    PURE module function getMergedDefCom_D1_PSSK1(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))                                      :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module function getMergedDefCom_D1_BSSK(sortedArray1, sortedArray2) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedDefCom_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_type)                                          :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! CusCom

    interface getMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getMergedCusCom_D0_SK5(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(len(sortedArray1,IK) + len(sortedArray2,IK))  :: mergedSortedArray
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if SK4_ENABLED
    module function getMergedCusCom_D0_SK4(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(len(sortedArray1,IK) + len(sortedArray2,IK))  :: mergedSortedArray
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if SK3_ENABLED
    module function getMergedCusCom_D0_SK3(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(len(sortedArray1,IK) + len(sortedArray2,IK))  :: mergedSortedArray
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if SK2_ENABLED
    module function getMergedCusCom_D0_SK2(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(len(sortedArray1,IK) + len(sortedArray2,IK))  :: mergedSortedArray
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if SK1_ENABLED
    module function getMergedCusCom_D0_SK1(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(len(sortedArray1,IK) + len(sortedArray2,IK))  :: mergedSortedArray
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getMergedCusCom_D1_SK5(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(len(sortedArray1,IK))                         :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if SK4_ENABLED
    module function getMergedCusCom_D1_SK4(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(len(sortedArray1,IK))                         :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if SK3_ENABLED
    module function getMergedCusCom_D1_SK3(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(len(sortedArray1,IK))                         :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if SK2_ENABLED
    module function getMergedCusCom_D1_SK2(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(len(sortedArray1,IK))                         :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if SK1_ENABLED
    module function getMergedCusCom_D1_SK1(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(len(sortedArray1,IK))                         :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getMergedCusCom_D1_IK5(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if IK4_ENABLED
    module function getMergedCusCom_D1_IK4(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if IK3_ENABLED
    module function getMergedCusCom_D1_IK3(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if IK2_ENABLED
    module function getMergedCusCom_D1_IK2(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if IK1_ENABLED
    module function getMergedCusCom_D1_IK1(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getMergedCusCom_D1_LK5(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if LK4_ENABLED
    module function getMergedCusCom_D1_LK4(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if LK3_ENABLED
    module function getMergedCusCom_D1_LK3(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if LK2_ENABLED
    module function getMergedCusCom_D1_LK2(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if LK1_ENABLED
    module function getMergedCusCom_D1_LK1(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getMergedCusCom_D1_CK5(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if CK4_ENABLED
    module function getMergedCusCom_D1_CK4(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if CK3_ENABLED
    module function getMergedCusCom_D1_CK3(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if CK2_ENABLED
    module function getMergedCusCom_D1_CK2(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if CK1_ENABLED
    module function getMergedCusCom_D1_CK1(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)                                            :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getMergedCusCom_D1_RK5(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)                                               :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if RK4_ENABLED
    module function getMergedCusCom_D1_RK4(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)                                               :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if RK3_ENABLED
    module function getMergedCusCom_D1_RK3(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)                                               :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if RK2_ENABLED
    module function getMergedCusCom_D1_RK2(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)                                               :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if RK1_ENABLED
    module function getMergedCusCom_D1_RK1(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)                                               :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module function getMergedCusCom_D1_PSSK5(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))                                      :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if SK4_ENABLED
    module function getMergedCusCom_D1_PSSK4(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))                                      :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if SK3_ENABLED
    module function getMergedCusCom_D1_PSSK3(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))                                      :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if SK2_ENABLED
    module function getMergedCusCom_D1_PSSK2(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))                                      :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#if SK1_ENABLED
    module function getMergedCusCom_D1_PSSK1(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))                                      :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module function getMergedCusCom_D1_BSSK(sortedArray1, sortedArray2, isSorted) result(mergedSortedArray)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMergedCusCom_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_type)                                          :: mergedSortedArray(size(sortedArray1) + size(sortedArray2))
        procedure(logical(LK))                                  :: isSorted
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Merge two ascending-sorted arrays such that the resulting merged array contains all elements of the two arrays in ascending order.
    !>
    !>  \param[in]  mergedSortedArray   :   The output scalar of either<br>
    !>                                      <ul>
    !>                                          <li>    type `character` of kind \SKALL or<br>
    !>                                      </ul>
    !>                                      the output `contiguous` array of rank `1` of either<br>
    !>                                      <ul>
    !>                                          <li>    type `character` of kind \SKALL or<br>
    !>                                          <li>    type `logical` of kind \LKALL or<br>
    !>                                          <li>    type `integer` of kind \IKALL or<br>
    !>                                          <li>    type `complex` of kind \CKALL or<br>
    !>                                          <li>    type `real` of kind \RKALL or<br>
    !>                                          <li>    type [css_pdt](@ref pm_container::css_pdt) of kind \SKALL or<br>
    !>                                          <li>    type [css_type](@ref pm_container::css_type) of default kind \SK,<br>
    !>                                      </ul>
    !>                                      whose size/length is the sum of the sizes/lengths of `sortedArray1` and `sortedArray2` and
    !>                                      whose elements are the combined elements of the input `sortedArray1` and `sortedArray2` in ascending order.
    !>  \param[in]  sortedArray1        :   The input `contiguous` array of the same type, kind, and rank as the
    !>                                      output `mergedSortedArray` argument whose elements **are already sorted** in ascending order.<br>
    !>                                      When the type of the input argument is `complex`, the array must be sorted only based on the real component of the values.
    !>  \param[in]  sortedArray2        :   The input `contiguous` array of the same type, kind, and rank as the
    !>                                      output `mergedSortedArray` argument whose elements **are already sorted** in ascending order.<br>
    !>                                      When the type of the input argument is `complex`, the array must be sorted only based on the real component of the values.
    !>  \param          isSorted        :   The `external` user-specified function that takes two input **scalar** arguments of the same type and kind as the input `Array`.<br>
    !>                                      It returns a scalar `logical` of default kind \LK that is `.true.` if the first
    !>                                      input scalar argument is sorted with respect to the second input argument according to the user-defined condition
    !>                                      within `isSorted`, otherwise, it is `.false.`.<br>
    !>                                      If `array` is a scalar string (i.e., an assumed-length scalar `character`),
    !>                                      then both input arguments to `isSorted()` are scalar characters of length `1` of kind \SKALL.<br>
    !>                                      The following illustrates the generic interface of `isSorted()`,
    !>                                      \code{.F90}
    !>                                          function isSorted(lhs, rhs) result(sorted)
    !>                                              use pm_kind, only: LK
    !>                                              TYPE(KIND)  , intent(in)    :: lhs, rhs
    !>                                              logical(LK)                 :: sorted
    !>                                          end function
    !>                                      \endcode
    !>                                      where `TYPE(KIND)` is the same as the type and kind of the input argument `Array`, which can be one of the following.
    !>                                      \code{.F90}
    !>                                          use pm_container, only: css_type, css_pdt
    !>                                          character(*, SK), intent(in) :: lhs, rhs
    !>                                          character(1, SK), intent(in) :: lhs, rhs
    !>                                          type(css_type)  , intent(in) :: lhs, rhs
    !>                                          type(css_pdt)   , intent(in) :: lhs, rhs
    !>                                          integer(IK)     , intent(in) :: lhs, rhs
    !>                                          logical(LK)     , intent(in) :: lhs, rhs
    !>                                          complex(CK)     , intent(in) :: lhs, rhs
    !>                                          real(RK)        , intent(in) :: lhs, rhs
    !>                                      \endcode
    !>                                      where the specified kind type parameters (`SK`, `IK`, `LK`, `CK`, `RK`) can refer to any of the supported kinds by the processor.<br>
    !>                                      This user-defined equivalence check is extremely useful where a user-defined sorting criterion other than simple ascending order
    !>                                      is needed, for example, when the case-sensitivity of an input string or array of strings is irrelevant or when sorting of
    !>                                      the absolute values matters excluding the signs of the numbers, or when descending order is desired.<br>
    !>                                      In such cases, user can define a custom sorting condition within the user-defined external function `isSorted` to achieve the goal.<br>
    !>                                      (**optional**, the default sorting condition is ascending order, that is `a < b`.)
    !>
    !>  \interface{setMerged}
    !>  \code{.F90}
    !>
    !>      use pm_arrayMerge, only: setMerged
    !>
    !>      call setMerged(mergedSortedArray, sortedArray1, sortedArray2)
    !>      call setMerged(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \warning
    !>  The condition `size(mergedSortedArray) == size(sortedArray1) + size(sortedArray2)` must hold for all interfaces with array-like arguments.<br>
    !>  The condition `len(mergedSortedArray) == len(sortedArray1) + len(sortedArray2)` must hold for all interfaces with scalar arguments of type `character`.<br>
    !>  The input arguments `sortedArray1` and `sortedArray2` must be sorted in ascending order when `isSorted()` is missing, or properly sorted when `isSorted()` is present.<br>
    !>  \vericons
    !>
    !>  \see
    !>  [getMerged](@ref pm_arrayMerge::getMerged)<br>
    !>
    !>  \example{setMerged}
    !>  \include{lineno} example/pm_arrayMerge/setMerged/main.F90
    !>  \compilef{setMerged}
    !>  \output{setMerged}
    !>  \include{lineno} example/pm_arrayMerge/setMerged/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayMerge](@ref test_pm_arrayMerge)
    !>
    !>  \final{setMerged}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! DefCom

    interface setMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMergedDefCom_D0_SK5(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(out)                   :: mergedSortedArray
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMergedDefCom_D0_SK4(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(*,SKG)        , intent(out)                   :: mergedSortedArray
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMergedDefCom_D0_SK3(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(*,SKG)        , intent(out)                   :: mergedSortedArray
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMergedDefCom_D0_SK2(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(*,SKG)        , intent(out)                   :: mergedSortedArray
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMergedDefCom_D0_SK1(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(*,SKG)        , intent(out)                   :: mergedSortedArray
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setMergedDefCom_D1_SK5(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(*,SKG)        , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMergedDefCom_D1_SK4(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(*,SKG)        , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMergedDefCom_D1_SK3(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(*,SKG)        , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMergedDefCom_D1_SK2(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(*,SKG)        , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMergedDefCom_D1_SK1(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(*,SKG)        , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMergedDefCom_D1_IK5(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMergedDefCom_D1_IK4(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMergedDefCom_D1_IK3(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMergedDefCom_D1_IK2(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMergedDefCom_D1_IK1(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setMergedDefCom_D1_LK5(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setMergedDefCom_D1_LK4(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setMergedDefCom_D1_LK3(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setMergedDefCom_D1_LK2(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setMergedDefCom_D1_LK1(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMergedDefCom_D1_CK5(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMergedDefCom_D1_CK4(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMergedDefCom_D1_CK3(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMergedDefCom_D1_CK2(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMergedDefCom_D1_CK1(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMergedDefCom_D1_RK5(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)               , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMergedDefCom_D1_RK4(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)               , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMergedDefCom_D1_RK3(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)               , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMergedDefCom_D1_RK2(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)               , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMergedDefCom_D1_RK1(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)               , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    PURE module subroutine setMergedDefCom_D1_PSSK5(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))      , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setMergedDefCom_D1_PSSK4(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))      , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setMergedDefCom_D1_PSSK3(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))      , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setMergedDefCom_D1_PSSK2(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))      , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setMergedDefCom_D1_PSSK1(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))      , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    PURE module subroutine setMergedDefCom_D1_BSSK(mergedSortedArray, sortedArray1, sortedArray2)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedDefCom_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_type)          , intent(out)   , contiguous    :: mergedSortedArray(:)
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! CusCom

    interface setMerged

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setMergedCusCom_D0_SK5(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D0_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(*,SKG)        , intent(out)                   :: mergedSortedArray
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setMergedCusCom_D0_SK4(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D0_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(*,SKG)        , intent(out)                   :: mergedSortedArray
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setMergedCusCom_D0_SK3(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D0_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(*,SKG)        , intent(out)                   :: mergedSortedArray
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setMergedCusCom_D0_SK2(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D0_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(*,SKG)        , intent(out)                   :: mergedSortedArray
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setMergedCusCom_D0_SK1(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D0_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)                    :: sortedArray1, sortedArray2
        character(*,SKG)        , intent(out)                   :: mergedSortedArray
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module subroutine setMergedCusCom_D1_SK5(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(*,SKG)        , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setMergedCusCom_D1_SK4(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(*,SKG)        , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setMergedCusCom_D1_SK3(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(*,SKG)        , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setMergedCusCom_D1_SK2(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(*,SKG)        , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setMergedCusCom_D1_SK1(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        character(*,SKG)        , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module subroutine setMergedCusCom_D1_IK5(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if IK4_ENABLED
    module subroutine setMergedCusCom_D1_IK4(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if IK3_ENABLED
    module subroutine setMergedCusCom_D1_IK3(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if IK2_ENABLED
    module subroutine setMergedCusCom_D1_IK2(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if IK1_ENABLED
    module subroutine setMergedCusCom_D1_IK1(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        integer(IKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module subroutine setMergedCusCom_D1_LK5(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if LK4_ENABLED
    module subroutine setMergedCusCom_D1_LK4(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if LK3_ENABLED
    module subroutine setMergedCusCom_D1_LK3(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if LK2_ENABLED
    module subroutine setMergedCusCom_D1_LK2(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if LK1_ENABLED
    module subroutine setMergedCusCom_D1_LK1(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        logical(LKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setMergedCusCom_D1_CK5(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setMergedCusCom_D1_CK4(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setMergedCusCom_D1_CK3(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setMergedCusCom_D1_CK2(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setMergedCusCom_D1_CK1(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        complex(CKG)            , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setMergedCusCom_D1_RK5(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)               , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setMergedCusCom_D1_RK4(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)               , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setMergedCusCom_D1_RK3(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)               , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setMergedCusCom_D1_RK2(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)               , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setMergedCusCom_D1_RK1(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        real(RKG)               , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if PDT_ENABLED

#if SK5_ENABLED
    module subroutine setMergedCusCom_D1_PSSK5(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_PSSK5
#endif
        use pm_kind, only: SKG => SK5
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))      , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK4_ENABLED
    module subroutine setMergedCusCom_D1_PSSK4(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_PSSK4
#endif
        use pm_kind, only: SKG => SK4
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))      , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK3_ENABLED
    module subroutine setMergedCusCom_D1_PSSK3(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_PSSK3
#endif
        use pm_kind, only: SKG => SK3
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))      , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK2_ENABLED
    module subroutine setMergedCusCom_D1_PSSK2(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_PSSK2
#endif
        use pm_kind, only: SKG => SK2
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))      , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#if SK1_ENABLED
    module subroutine setMergedCusCom_D1_PSSK1(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_PSSK1
#endif
        use pm_kind, only: SKG => SK1
        use pm_container, only: css_pdt
        type(css_pdt(SKG))      , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_pdt(SKG))      , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine
#endif

#endif
!PDT_ENABLED

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    module subroutine setMergedCusCom_D1_BSSK(mergedSortedArray, sortedArray1, sortedArray2, isSorted)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMergedCusCom_D1_BSSK
#endif
        use pm_kind, only: SKG => SK
        use pm_container, only: css_type
        type(css_type)          , intent(in)    , contiguous    :: sortedArray1(:), sortedArray2(:)
        type(css_type)          , intent(out)   , contiguous    :: mergedSortedArray(:)
        procedure(logical(LK))                                  :: isSorted
    end subroutine

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayMerge