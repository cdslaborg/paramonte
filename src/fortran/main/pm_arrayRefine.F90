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

!>  \cond excluded
!   \bug
!   The following bypasses the bug reported below that creates a conflict between Intel and gfortran.
#if     __INTEL_COMPILER && 0
#define LEN_STR :
#else
#define LEN_STR len(array,IK)
#endif
!>  \endcond excluded

!>  \brief
!>  This module contains procedures and generic interfaces for refining (thinning) (weighted) arrays of arbitrary intrinsic types.<br>
!>
!>  \details
!>  Refinement in the context of this module means skipping through (weighted) array elements by a certain `skip` amount.<br>
!>  Refining unweighted arrays is straightforward in Fortran as there is an intrinsic slicing syntax for it.<br>
!>  However, the task can become cumbersome for weighted arrays.<br>
!>  This module aims to facilitate refinement of weighted arrays.<br>
!>
!>  \note
!>  See [pm_arrayCopy](@ref pm_arrayCopy) for refining unweighted strings and arrays.<br>
!>
!>  \see
!>  [pm_sampleWeight](@ref pm_sampleWeight)<br>
!>  [pm_arrayCompact](@ref pm_arrayCompact)<br>
!>  [pm_arrayVerbose](@ref pm_arrayVerbose)<br>
!>  [pm_arrayRefine](@ref pm_arrayRefine)<br>
!>  [pm_arrayCopy](@ref pm_arrayCopy)<br>
!>
!>  \test
!>  [test_pm_arrayRefine](@ref test_pm_arrayRefine)
!>
!>  \todo
!>  \pmed
!>  Interfaces for `real` weights and without weights should be added in future.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Saturday 1:48 AM, August 20, 2016, Institute for Computational Engineering and Sciences, UT Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayRefine

    use pm_kind, only: SK, IK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_arrayRefine"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate a refined version of the input `array` by the specified `weight` and `skip`.<br>
    !>
    !>  \param[in]  array       :   The input `contiguous` array of shape `(:)` or `(:,:)` of either
    !>                              <ol>
    !>                                  <li>    type `character` of kind \SKALL or
    !>                                  <li>    type `integer` of kind \IKALL or
    !>                                  <li>    type `logical` of kind \LKALL or
    !>                                  <li>    type `complex` of kind \CKALL or
    !>                                  <li>    type `real` of kind \RKALL
    !>                              </ol>
    !>                              or,
    !>                              <ol>
    !>                                  <li>    scalar `character` of kind \SKALL,
    !>                              </ol>
    !>                              containing the array that has to be refined.
    !>  \param[in]  dim         :   The input scalar of type `integer` of default kind \IK representing the axis of `array(:,:)` along which `array` must be refined.<br>
    !>                              (**optional**, it must be present <b>if and only if</b> `array` is of shape `(:,:)`.)
    !>  \param[in]  weight      :   The input vector of<br>
    !>                              <ol>
    !>                                  <li>    type `integer` of default kind \IK,
    !>                              </ol>
    !>                              of size `size(array, dim)` containing the weights of individual elements of the input array.<br>
    !>  \param[in]  skip        :   The input scalar of type `integer` of default kind \IK representing the number of elements to skip in the input sequence.<br>
    !>
    !>  \return
    !>  `arrayRefined`          :   The output `allocatable` array of the same type, kind, and shape as the input `array`,
    !>                              containing the refined **unweighted** version of the input `array`.<br>
    !>                              The returned array is unweighted to preserve the purity of the procedure.<br>
    !>                              See [setRefined](@ref pm_arrayRefine::setRefined) for an alternative interface.<br>
    !>
    !>  \interface{getRefined}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRefine, only: getRefined
    !>
    !>      arrayRefined = getRefined(array, weight, skip) ! scalar character objects.
    !>      arrayRefined(:) = getRefined(array(:), weight, skip) ! all intrinsic array objects.
    !>      arrayRefined(:,:) = getRefined(array(:,:), dim, weight, skip) ! all intrinsic array objects.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < skip` must hold.<br>
    !>  The condition `dim == 1 .or. dim == 2` must hold.<br>
    !>  The condition `size(weight) == size(array)` must hold when `array` is rank `1`.<br>
    !>  The condition `size(weight) == size(array, dim)` must hold when `array` is rank `2`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  See [pm_arrayCopy](@ref pm_arrayCopy) for refining unweighted strings and arrays.<br>
    !>
    !>  \see
    !>  [setRefined](@ref pm_arrayRefine::setRefined)<br>
    !>  [getCompact](@ref pm_arrayCompact::getCompact)<br>
    !>  [getVerbose](@ref pm_arrayVerbose::getVerbose)<br>
    !>  [pm_arrayCopy](@ref pm_arrayCopy)<br>
    !>
    !>  \example{getRefined}
    !>  \include{lineno} example/pm_arrayRefine/getRefined/main.F90
    !>  \compilef{getRefined}
    !>  \output{getRefined}
    !>  \include{lineno} example/pm_arrayRefine/getRefined/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRefine](@ref test_pm_arrayRefine)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{10.3-12}, \ifort{2021-2022}
    !>  \desc
    !>  \ifort and \gfortran share a common bug with opposing behavior.<br>
    !>  \ifort cannot handle assumed-length `allocatable` output arguments of type `character`.<br>
    !>  \gfortran cannot handle deferred-length `allocatable` output arguments of type `character`.<br>
    !>  \remedy
    !>  For now, a preprocessor macro defines two separate interfaces for the two compilers so that both compilers can compile this file.<br>
    !>  This minor interface difference should not impact the usage of this module with different compilers.<br>
    !>
    !>  \final{getRefined}
    !>
    !>  \author
    !>  \AmirShahmoradi, Saturday 1:48 AM, August 20, 2016, Institute for Computational Engineering and Sciences, UT Austin, TX
    interface getRefined

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRefined_WTI_D0_SK5(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)                    :: array
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        character(:,SKC)                            , allocatable   :: arrayRefined
    end function
#endif

#if SK4_ENABLED
    PURE module function getRefined_WTI_D0_SK4(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)                    :: array
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        character(:,SKC)                            , allocatable   :: arrayRefined
    end function
#endif

#if SK3_ENABLED
    PURE module function getRefined_WTI_D0_SK3(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)                    :: array
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        character(:,SKC)                            , allocatable   :: arrayRefined
    end function
#endif

#if SK2_ENABLED
    PURE module function getRefined_WTI_D0_SK2(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)                    :: array
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        character(:,SKC)                            , allocatable   :: arrayRefined
    end function
#endif

#if SK1_ENABLED
    PURE module function getRefined_WTI_D0_SK1(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)                    :: array
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        character(:,SKC)                            , allocatable   :: arrayRefined
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRefined_WTI_D1_SK5(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        character(LEN_STR,SKC)                      , allocatable   :: arrayRefined(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getRefined_WTI_D1_SK4(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        character(LEN_STR,SKC)                      , allocatable   :: arrayRefined(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getRefined_WTI_D1_SK3(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        character(LEN_STR,SKC)                      , allocatable   :: arrayRefined(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getRefined_WTI_D1_SK2(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        character(LEN_STR,SKC)                      , allocatable   :: arrayRefined(:)
    end function
#endif

#if SK1_ENABLED
    PURE module function getRefined_WTI_D1_SK1(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        character(LEN_STR,SKC)                      , allocatable   :: arrayRefined(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getRefined_WTI_D1_IK5(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IKC)                                , allocatable   :: arrayRefined(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getRefined_WTI_D1_IK4(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IKC)                                , allocatable   :: arrayRefined(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getRefined_WTI_D1_IK3(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IKC)                                , allocatable   :: arrayRefined(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getRefined_WTI_D1_IK2(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IKC)                                , allocatable   :: arrayRefined(:)
    end function
#endif

#if IK1_ENABLED
    PURE module function getRefined_WTI_D1_IK1(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IKC)                                , allocatable   :: arrayRefined(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getRefined_WTI_D1_LK5(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        logical(LKC)                                , allocatable   :: arrayRefined(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getRefined_WTI_D1_LK4(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        logical(LKC)                                , allocatable   :: arrayRefined(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getRefined_WTI_D1_LK3(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        logical(LKC)                                , allocatable   :: arrayRefined(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getRefined_WTI_D1_LK2(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        logical(LKC)                                , allocatable   :: arrayRefined(:)
    end function
#endif

#if LK1_ENABLED
    PURE module function getRefined_WTI_D1_LK1(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        logical(LKC)                                , allocatable   :: arrayRefined(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getRefined_WTI_D1_CK5(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        complex(CKC)                                , allocatable   :: arrayRefined(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getRefined_WTI_D1_CK4(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        complex(CKC)                                , allocatable   :: arrayRefined(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getRefined_WTI_D1_CK3(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        complex(CKC)                                , allocatable   :: arrayRefined(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getRefined_WTI_D1_CK2(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        complex(CKC)                                , allocatable   :: arrayRefined(:)
    end function
#endif

#if CK1_ENABLED
    PURE module function getRefined_WTI_D1_CK1(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        complex(CKC)                                , allocatable   :: arrayRefined(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getRefined_WTI_D1_RK5(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        real(RKC)                                   , allocatable   :: arrayRefined(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getRefined_WTI_D1_RK4(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        real(RKC)                                   , allocatable   :: arrayRefined(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getRefined_WTI_D1_RK3(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        real(RKC)                                   , allocatable   :: arrayRefined(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getRefined_WTI_D1_RK2(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        real(RKC)                                   , allocatable   :: arrayRefined(:)
    end function
#endif

#if RK1_ENABLED
    PURE module function getRefined_WTI_D1_RK1(array, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(in)    , contiguous    :: array(:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        real(RKC)                                   , allocatable   :: arrayRefined(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getRefined_WTI_D2_SK5(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        character(LEN_STR,SKC)                      , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getRefined_WTI_D2_SK4(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        character(LEN_STR,SKC)                      , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getRefined_WTI_D2_SK3(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        character(LEN_STR,SKC)                      , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getRefined_WTI_D2_SK2(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        character(LEN_STR,SKC)                      , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if SK1_ENABLED
    PURE module function getRefined_WTI_D2_SK1(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        character(LEN_STR,SKC)                      , allocatable   :: arrayRefined(:,:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getRefined_WTI_D2_IK5(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        integer(IKC)                                , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getRefined_WTI_D2_IK4(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        integer(IKC)                                , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getRefined_WTI_D2_IK3(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        integer(IKC)                                , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getRefined_WTI_D2_IK2(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        integer(IKC)                                , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if IK1_ENABLED
    PURE module function getRefined_WTI_D2_IK1(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        integer(IKC)                                , allocatable   :: arrayRefined(:,:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getRefined_WTI_D2_LK5(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        logical(LKC)                                , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getRefined_WTI_D2_LK4(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        logical(LKC)                                , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getRefined_WTI_D2_LK3(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        logical(LKC)                                , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getRefined_WTI_D2_LK2(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        logical(LKC)                                , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if LK1_ENABLED
    PURE module function getRefined_WTI_D2_LK1(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        logical(LKC)                                , allocatable   :: arrayRefined(:,:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getRefined_WTI_D2_CK5(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        complex(CKC)                                , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getRefined_WTI_D2_CK4(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        complex(CKC)                                , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getRefined_WTI_D2_CK3(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        complex(CKC)                                , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getRefined_WTI_D2_CK2(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        complex(CKC)                                , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if CK1_ENABLED
    PURE module function getRefined_WTI_D2_CK1(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        complex(CKC)                                , allocatable   :: arrayRefined(:,:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getRefined_WTI_D2_RK5(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        real(RKC)                                   , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getRefined_WTI_D2_RK4(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        real(RKC)                                   , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getRefined_WTI_D2_RK3(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        real(RKC)                                   , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getRefined_WTI_D2_RK2(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        real(RKC)                                   , allocatable   :: arrayRefined(:,:)
    end function
#endif

#if RK1_ENABLED
    PURE module function getRefined_WTI_D2_RK1(array, dim, weight, skip) result(arrayRefined)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRefined_WTI_D2_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(in)    , contiguous    :: array(:,:)
        integer(IK)                 , intent(in)    , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(in)                    :: dim
        real(RKC)                                   , allocatable   :: arrayRefined(:,:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface getRefined

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate a refined version of the input `array` where the **sequentially** unweighted entries
    !>  along the specified dimension of `array` are skipped every `skip` to create a refined weighted output array of size `rsize`.<br>
    !>
    !>  \param[inout]   array       :   The input/output `contiguous` array of shape `(:)` or `(:,:)` of either
    !>                                  <ol>
    !>                                      <li>    type `character` of kind \SKALL or
    !>                                      <li>    type `integer` of kind \IKALL or
    !>                                      <li>    type `logical` of kind \LKALL or
    !>                                      <li>    type `complex` of kind \CKALL or
    !>                                      <li>    type `real` of kind \RKALL
    !>                                  </ol>
    !>                                  or scalar `character` of kind \SKALL.<br>
    !>                                  On input, it contains the elements to be refined.<br>
    !>                                  On output, the first `rsize` entries of array contain the condensed elements of the original input `array`,
    !>                                  such that **no** two adjacent entries in `array` from entry `1` to entry `rsize` are duplicates.
    !>  \param[in]      dim         :   The input scalar of type `integer` of default kind \IK representing the axis of `array(:,:)` along which `array` must be refined.<br>
    !>                                  (**optional**, it must be present <b>if and only if</b> `array` is of shape `(:,:)`.)
    !>  \param[inout]   weight      :   The input/output vector of<br>
    !>                                  <ol>
    !>                                      <li>    type `integer` of default kind \IK,
    !>                                  </ol>
    !>                                  of size `size(array, dim)` containing the weights of individual elements of the input array.<br>
    !>                                  On output, the first `rsize` elements of `weight` are written with the new refined weights of the corresponding elements in the output `array`.<br>
    !>  \param[in]      skip        :   The input scalar of type `integer` of default kind \IK representing the number of elements to skip in the input sequence.<br>
    !>  \param[out]     rsize       :   The output scalar of type `integer` of default kind \IK such that,<br>
    !>                                  <ol>
    !>                                      <li>    `weight(1:rsize)` contains the weights of the refined elements of the input `array`, and<br>
    !>                                      <li>    `array(1:rsize)` or `array(:,1:rsize)` or `array(1:rsize,:)` contains the refined elements of the input `array`.<br>
    !>                                  </ol>
    !>
    !>  \interface{setRefined}
    !>  \code{.F90}
    !>
    !>      use pm_arrayRefine, only: setRefined
    !>
    !>      call setRefined(array, weight(:), rsize) ! scalar character objects.
    !>      call setRefined(array(:), weight(:), rsize) ! all intrinsic array objects.
    !>      call setRefined(array(:,:), dim, weight(:), rsize) ! all intrinsic array objects.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < skip` must hold.<br>
    !>  The condition `dim == 1 .or. dim == 2` must hold.<br>
    !>  The condition `size(weight) == size(array)` must hold when `array` is rank `1`.<br>
    !>  The condition `size(weight) == size(array, dim)` must hold when `array` is rank `2`.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  See [pm_arrayCopy](@ref pm_arrayCopy) for refining unweighted strings and arrays.<br>
    !>
    !>  \see
    !>  [setRefined](@ref pm_arrayRefine::setRefined)<br>
    !>  [getCompact](@ref pm_arrayCompact::getCompact)<br>
    !>  [getVerbose](@ref pm_arrayVerbose::getVerbose)<br>
    !>  [pm_arrayCopy](@ref pm_arrayCopy)<br>
    !>
    !>  \example{setRefined}
    !>  \include{lineno} example/pm_arrayRefine/setRefined/main.F90
    !>  \compilef{setRefined}
    !>  \output{setRefined}
    !>  \include{lineno} example/pm_arrayRefine/setRefined/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayRefine](@ref test_pm_arrayRefine)
    !>
    !>  \final{setRefined}
    !>
    !>  \author
    !>  \AmirShahmoradi, Saturday 1:48 AM, August 20, 2016, Institute for Computational Engineering and Sciences, UT Austin, TX

    interface setRefined

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRefined_WTI_D0_SK5(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(inout)                 :: array
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefined_WTI_D0_SK4(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(inout)                 :: array
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefined_WTI_D0_SK3(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(inout)                 :: array
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefined_WTI_D0_SK2(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(inout)                 :: array
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefined_WTI_D0_SK1(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(inout)                 :: array
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRefined_WTI_D1_SK5(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefined_WTI_D1_SK4(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefined_WTI_D1_SK3(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefined_WTI_D1_SK2(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefined_WTI_D1_SK1(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRefined_WTI_D1_IK5(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRefined_WTI_D1_IK4(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRefined_WTI_D1_IK3(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRefined_WTI_D1_IK2(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRefined_WTI_D1_IK1(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRefined_WTI_D1_LK5(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRefined_WTI_D1_LK4(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRefined_WTI_D1_LK3(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRefined_WTI_D1_LK2(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRefined_WTI_D1_LK1(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRefined_WTI_D1_CK5(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRefined_WTI_D1_CK4(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRefined_WTI_D1_CK3(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRefined_WTI_D1_CK2(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRefined_WTI_D1_CK1(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRefined_WTI_D1_RK5(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRefined_WTI_D1_RK4(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRefined_WTI_D1_RK3(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRefined_WTI_D1_RK2(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRefined_WTI_D1_RK1(array, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(inout) , contiguous    :: array(:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module subroutine setRefined_WTI_D2_SK5(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)            , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if SK4_ENABLED
    PURE module subroutine setRefined_WTI_D2_SK4(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)            , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if SK3_ENABLED
    PURE module subroutine setRefined_WTI_D2_SK3(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)            , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if SK2_ENABLED
    PURE module subroutine setRefined_WTI_D2_SK2(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)            , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if SK1_ENABLED
    PURE module subroutine setRefined_WTI_D2_SK1(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)            , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setRefined_WTI_D2_IK5(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setRefined_WTI_D2_IK4(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setRefined_WTI_D2_IK3(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setRefined_WTI_D2_IK2(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setRefined_WTI_D2_IK1(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module subroutine setRefined_WTI_D2_LK5(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if LK4_ENABLED
    PURE module subroutine setRefined_WTI_D2_LK4(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if LK3_ENABLED
    PURE module subroutine setRefined_WTI_D2_LK3(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if LK2_ENABLED
    PURE module subroutine setRefined_WTI_D2_LK2(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if LK1_ENABLED
    PURE module subroutine setRefined_WTI_D2_LK1(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setRefined_WTI_D2_CK5(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setRefined_WTI_D2_CK4(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setRefined_WTI_D2_CK3(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setRefined_WTI_D2_CK2(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setRefined_WTI_D2_CK1(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)                , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setRefined_WTI_D2_RK5(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)                   , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setRefined_WTI_D2_RK4(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)                   , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setRefined_WTI_D2_RK3(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)                   , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setRefined_WTI_D2_RK2(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)                   , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setRefined_WTI_D2_RK1(array, dim, weight, skip, rsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setRefined_WTI_D2_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)                   , intent(inout) , contiguous    :: array(:,:)
        integer(IK)                 , intent(inout) , contiguous    :: weight(:)
        integer(IK)                 , intent(in)                    :: skip
        integer(IK)                 , intent(out)                   :: rsize
        integer(IK)                 , intent(in)                    :: dim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setRefined

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayRefine ! LCOV_EXCL_LINE