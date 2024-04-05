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

!       \bug
!       The following bypasses the bug reported below that creates a conflict between Intel and gfortran.
#if     __INTEL_COMPILER
#define LEN_STR :
#else
#define LEN_STR len(setA,IK)
#endif

!>  \brief
!>  This module contains procedures and generic interfaces for computing the **absolute or relative complement** of one set in another set.
!>
!>  \details
!>  In mathematical terminology, the **complement** of set \f$A\f$ in set \f$B\f$
!>  is denoted by \f$B \backslash A\f$, or \f$B \cap A^{C}\f$, or less frequently by \f$B - A\f$.<br>
!>  Read \f$B \backslash A\f$ as *the set of elements in B that do not appear in A*.<br>
!>  <ol>
!>      <li>    If the set \f$B\f$ **contains** all elements of the set \f$A\f$, then \f$B \backslash A\f$ is the **absolute complement** of \f$A\f$.
!>      <li>    If the set \f$B\f$ **does not contain** all elements of the set \f$A\f$, then \f$B \backslash A\f$ is the **relative complement** of \f$A\f$ (with respect to \f$B\f$).
!>  </ol>
!>  An example Venn diagram illustration of \f$B \backslash A\f$ is the red region in the following graph,<br>
!>
!>  \image html pm_arrayComplement@RelativeComplement.png width=500
!>
!>  \note
!>  This module contains two generic interfaces for computing the difference of two unsorted or two similarly-sorted sets.<br>
!>  If the sets are similarly sorted, use [getComplement](@ref pm_arrayComplement::getComplement) with `sorted = .true.` input argument.<br>
!>  Computing the difference of two sorted sets is significantly faster than computing the difference of two similarly-sorted sets.<br>
!>
!>  \devnote
!>  The two generic interfaces of this module [getComplement](@ref pm_arrayComplement::getComplement)
!>  and [getComplementRange](@ref pm_arrayComplement::getComplementRange) are
!>  intentionally kept separately from each other under different names.<br>
!>
!>  \test
!>  [test_pm_arrayComplement](@ref test_pm_arrayComplement)
!>
!>  \todo
!>  \phigh
!>  The two function interfaces of this module can be merged into a single generic interface.<br>
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:24 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayComplement

    use pm_kind, only: SK, IK, LK
    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_arrayComplement"

    !>  \brief
    !>  Generate and return an array of elements of `setB` that are not in `setA`.
    !>
    !>  \param[in]  setA        :   The input `contiguous` vector of either,
    !>                              <ul>
    !>                                  <li> type `character` of kind \SKALL, or
    !>                                  <li> type `integer` of kind \IKALL, or
    !>                                  <li> type `logical` of kind \LKALL, or
    !>                                  <li> type `complex` of kind \CKALL, or
    !>                                  <li> type `real` of kind \RKALL, or
    !>                              </ul>
    !>                              or,
    !>                              <ul>
    !>                                  <li>    scalar `character` of arbitrary length of kind \SKALL,<br>
    !>                              </ul>
    !>                              whose complement in `setB` will be returned as the output `complement`.
    !>  \param[in]  setB        :   The input object of the same type, kind, and rank as `setA` representing the set with respect to which the complement of `setA` is computed.
    !>  \param[in]  sorted      :   The input scalar `logical` of default kind \LK.
    !>                              <ol>
    !>                                  <li> If `.false.`, the input sets `setA` and `setB` are assumed to be dissimilarly sorted or not sorted at all.
    !>                                  <li> If `.true.`, the input sets `setA` and `setB` are assumed to be **similarly sorted** (e.g., both ascending or both descending order according to an arbitrary criterion).
    !>                              </ol>
    !>                              If the input sets are similarly sorted, then specifying `sorted = .true._LK` can lead to significantly better runtime performance.<br>
    !>                              (**optional**, default = `.false._LK`)
    !>  \param[in]  unique      :   The input scalar `logical` of default kind \LK.
    !>                              <ol>
    !>                                  <li> If `.false.`, each of the input sets `setA` and `setB` are assumed to possibly contain duplicate elements individually.
    !>                                  <li> If `.true.`, all elements within each of the input sets `setA` and `setB` are assumed to be unique.
    !>                              </ol>
    !>                              If the elements of each of the input sets are unique, then specifying `unique = .true._LK` can lead to significantly better runtime performance.<br>
    !>                              The specified value for `unique` becomes relevant only if `sorted = .true.`. Its value is ignored when `sorted = .false.`<br>
    !>                              (**optional**, default = `.false._LK`).
    !>  \param      iseq        :   The `external` user-specified function that takes two input **scalar** arguments of the same type and kind as the input `setA`.<br>
    !>                              The first input argument is an element of `setA` to be compared with the second argument that is an element from `setB`.<br>
    !>                              If `setA` and `setB` are scalars of type `character`, the `len` type-parameter of both input arguments to `iseq` is `1`.<br>
    !>                              It returns a scalar `logical` of default kind \LK that is `.true.` if the two input arguments are equivalent (e.g., equal)
    !>                              according to the user-defined criterion, otherwise, it is `.false.`.<br>
    !>                              The following illustrates the generic interface of `iseq`,
    !>                              \code{.F90}
    !>                                  function iseq(elementA, elementB) result(equivalent)
    !>                                      use pm_kind, only: LK
    !>                                      TYPE(KIND)  , intent(in)    :: elementA, elementB
    !>                                      logical(LK)                 :: equivalent
    !>                                  end function
    !>                              \endcode
    !>                              where `TYPE(KIND)` represents the type and kind of the input argument `setA`, which can be one of the following,
    !>                              \code{.F90}
    !>                                  use pm_kind, only: SK, IK, CK, RK
    !>                                  character(*, SK), intent(in)    :: elementA, elementB ! when `array` is a string vector.
    !>                                  character(1, SK), intent(in)    :: elementA, elementB ! when `array` is a string scalar.
    !>                                  integer(IK)     , intent(in)    :: elementA, elementB
    !>                                  complex(CK)     , intent(in)    :: elementA, elementB
    !>                                  real(RK)        , intent(in)    :: elementA, elementB
    !>                              \endcode
    !>                              where the kinds `SK`, `IK`, `CK`, `RK`, can refer to any kind type parameter that is supported by the processor.<br>
    !>                              This user-defined equivalence check is extremely useful where a user-defined equivalence test other than exact equality
    !>                              or identity is needed, for example, when the set elements should match only within a given threshold or,
    !>                              when the case-sensitivity in character comparisons do not matter.<br>
    !>                              In such cases, the user can define a custom equivalence criterion within the user-defined external function `iseq` to achieve the goal.<br>
    !>                              (**optional**, the default equivalence operator is `==`)
    !>
    !>  \return
    !>  `complement`            :   The output `allocatable` object of the same type, kind, and rank as `setA` containing
    !>                              the complement of `setA` in `setB` (i.e., the elements of `setB` that are not in `setA`).
    !>
    !>  \interface{getComplement}
    !>  \code{.F90}
    !>
    !>      use pm_arrayComplement, only: getComplement
    !>
    !>      complement = getComplement(setA, setB) ! scalar assumed-length string arguments and output.
    !>      complement(:) = getComplement(setA(:), setB(:))
    !>      complement(:) = getComplement(setA(:), setB(:), iseq)
    !>      complement(:) = getComplement(setA(:), setB(:), sorted, unique)
    !>      complement(:) = getComplement(setA(:), setB(:), sorted, unique, iseq)
    !>
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  The performance of the procedures under this generic interface can dramatically improve if the input sets contain unique values (`unique = .true._LK`) and both sets are similarly-sorted (`sorted = .true._LK`).<br>
    !>  The unique elements of an arbitrary set can be obtained via the procedures of [pm_arrayUnique](@ref pm_arrayUnique).
    !>
    !>  \see
    !>  [getComplementRange](@ref pm_arrayComplement::getComplementRange)<br>
    !>
    !>  \example{getComplement}
    !>  \include{lineno} example/pm_arrayComplement/getComplement/main.F90
    !>  \compilef{getComplement}
    !>  \output{getComplement}
    !>  \include{lineno} example/pm_arrayComplement/getComplement/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayComplement](@ref test_pm_arrayComplement)
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
    !>  \finmain{getComplement}
    !>
    !>  \author
    !>  \FatemehBagheri, Wednesday 1:35 PM, August 11, 2021, Dallas, TX
    interface getComplement

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getComplementRandomDefCom_D0_SK5(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

#if SK4_ENABLED
    PURE module function getComplementRandomDefCom_D0_SK4(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

#if SK3_ENABLED
    PURE module function getComplementRandomDefCom_D0_SK3(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

#if SK2_ENABLED
    PURE module function getComplementRandomDefCom_D0_SK2(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

#if SK1_ENABLED
    PURE module function getComplementRandomDefCom_D0_SK1(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getComplementRandomDefCom_D1_SK5(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getComplementRandomDefCom_D1_SK4(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getComplementRandomDefCom_D1_SK3(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getComplementRandomDefCom_D1_SK2(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

#if SK1_ENABLED
    PURE module function getComplementRandomDefCom_D1_SK1(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getComplementRandomDefCom_D1_IK5(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getComplementRandomDefCom_D1_IK4(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getComplementRandomDefCom_D1_IK3(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getComplementRandomDefCom_D1_IK2(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK1_ENABLED
    PURE module function getComplementRandomDefCom_D1_IK1(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getComplementRandomDefCom_D1_LK5(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getComplementRandomDefCom_D1_LK4(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getComplementRandomDefCom_D1_LK3(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getComplementRandomDefCom_D1_LK2(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

#if LK1_ENABLED
    PURE module function getComplementRandomDefCom_D1_LK1(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getComplementRandomDefCom_D1_CK5(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getComplementRandomDefCom_D1_CK4(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getComplementRandomDefCom_D1_CK3(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getComplementRandomDefCom_D1_CK2(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

#if CK1_ENABLED
    PURE module function getComplementRandomDefCom_D1_CK1(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getComplementRandomDefCom_D1_RK5(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getComplementRandomDefCom_D1_RK4(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getComplementRandomDefCom_D1_RK3(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getComplementRandomDefCom_D1_RK2(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

#if RK1_ENABLED
    PURE module function getComplementRandomDefCom_D1_RK1(setA, setB) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomDefCom_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getComplementRandomCusCom_D0_SK5(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        procedure(logical(LK))                              :: iseq
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

#if SK4_ENABLED
    module function getComplementRandomCusCom_D0_SK4(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        procedure(logical(LK))                              :: iseq
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

#if SK3_ENABLED
    module function getComplementRandomCusCom_D0_SK3(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        procedure(logical(LK))                              :: iseq
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

#if SK2_ENABLED
    module function getComplementRandomCusCom_D0_SK2(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        procedure(logical(LK))                              :: iseq
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

#if SK1_ENABLED
    module function getComplementRandomCusCom_D0_SK1(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        procedure(logical(LK))                              :: iseq
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getComplementRandomCusCom_D1_SK5(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

#if SK4_ENABLED
    module function getComplementRandomCusCom_D1_SK4(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

#if SK3_ENABLED
    module function getComplementRandomCusCom_D1_SK3(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

#if SK2_ENABLED
    module function getComplementRandomCusCom_D1_SK2(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

#if SK1_ENABLED
    module function getComplementRandomCusCom_D1_SK1(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getComplementRandomCusCom_D1_IK5(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK4_ENABLED
    module function getComplementRandomCusCom_D1_IK4(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK3_ENABLED
    module function getComplementRandomCusCom_D1_IK3(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK2_ENABLED
    module function getComplementRandomCusCom_D1_IK2(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK1_ENABLED
    module function getComplementRandomCusCom_D1_IK1(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getComplementRandomCusCom_D1_LK5(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

#if LK4_ENABLED
    module function getComplementRandomCusCom_D1_LK4(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

#if LK3_ENABLED
    module function getComplementRandomCusCom_D1_LK3(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

#if LK2_ENABLED
    module function getComplementRandomCusCom_D1_LK2(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

#if LK1_ENABLED
    module function getComplementRandomCusCom_D1_LK1(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getComplementRandomCusCom_D1_CK5(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

#if CK4_ENABLED
    module function getComplementRandomCusCom_D1_CK4(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

#if CK3_ENABLED
    module function getComplementRandomCusCom_D1_CK3(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

#if CK2_ENABLED
    module function getComplementRandomCusCom_D1_CK2(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

#if CK1_ENABLED
    module function getComplementRandomCusCom_D1_CK1(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getComplementRandomCusCom_D1_RK5(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

#if RK4_ENABLED
    module function getComplementRandomCusCom_D1_RK4(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

#if RK3_ENABLED
    module function getComplementRandomCusCom_D1_RK3(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

#if RK2_ENABLED
    module function getComplementRandomCusCom_D1_RK2(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

#if RK1_ENABLED
    module function getComplementRandomCusCom_D1_RK1(setA, setB, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRandomCusCom_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        procedure(logical(LK))                              :: iseq
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getComplementSortedDefCom_D0_SK5(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        logical(LK)             , intent(in)                :: sorted, unique
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

#if SK4_ENABLED
    PURE module function getComplementSortedDefCom_D0_SK4(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        logical(LK)             , intent(in)                :: sorted, unique
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

#if SK3_ENABLED
    PURE module function getComplementSortedDefCom_D0_SK3(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        logical(LK)             , intent(in)                :: sorted, unique
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

#if SK2_ENABLED
    PURE module function getComplementSortedDefCom_D0_SK2(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        logical(LK)             , intent(in)                :: sorted, unique
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

#if SK1_ENABLED
    PURE module function getComplementSortedDefCom_D0_SK1(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        logical(LK)             , intent(in)                :: sorted, unique
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function getComplementSortedDefCom_D1_SK5(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted, unique
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

#if SK4_ENABLED
    PURE module function getComplementSortedDefCom_D1_SK4(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted, unique
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

#if SK3_ENABLED
    PURE module function getComplementSortedDefCom_D1_SK3(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted, unique
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

#if SK2_ENABLED
    PURE module function getComplementSortedDefCom_D1_SK2(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted, unique
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

#if SK1_ENABLED
    PURE module function getComplementSortedDefCom_D1_SK1(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted, unique
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getComplementSortedDefCom_D1_IK5(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted, unique
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getComplementSortedDefCom_D1_IK4(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted, unique
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getComplementSortedDefCom_D1_IK3(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getComplementSortedDefCom_D1_IK2(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK1_ENABLED
    PURE module function getComplementSortedDefCom_D1_IK1(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function getComplementSortedDefCom_D1_LK5(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted, unique
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

#if LK4_ENABLED
    PURE module function getComplementSortedDefCom_D1_LK4(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted, unique
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

#if LK3_ENABLED
    PURE module function getComplementSortedDefCom_D1_LK3(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

#if LK2_ENABLED
    PURE module function getComplementSortedDefCom_D1_LK2(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

#if LK1_ENABLED
    PURE module function getComplementSortedDefCom_D1_LK1(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getComplementSortedDefCom_D1_CK5(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getComplementSortedDefCom_D1_CK4(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getComplementSortedDefCom_D1_CK3(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getComplementSortedDefCom_D1_CK2(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

#if CK1_ENABLED
    PURE module function getComplementSortedDefCom_D1_CK1(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getComplementSortedDefCom_D1_RK5(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getComplementSortedDefCom_D1_RK4(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getComplementSortedDefCom_D1_RK3(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getComplementSortedDefCom_D1_RK2(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

#if RK1_ENABLED
    PURE module function getComplementSortedDefCom_D1_RK1(setA, setB, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedDefCom_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getComplementSortedCusCom_D0_SK5(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        logical(LK)             , intent(in)                :: sorted, unique
        procedure(logical(LK))                              :: iseq
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

#if SK4_ENABLED
    module function getComplementSortedCusCom_D0_SK4(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        logical(LK)             , intent(in)                :: sorted, unique
        procedure(logical(LK))                              :: iseq
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

#if SK3_ENABLED
    module function getComplementSortedCusCom_D0_SK3(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        logical(LK)             , intent(in)                :: sorted, unique
        procedure(logical(LK))                              :: iseq
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

#if SK2_ENABLED
    module function getComplementSortedCusCom_D0_SK2(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        logical(LK)             , intent(in)                :: sorted, unique
        procedure(logical(LK))                              :: iseq
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

#if SK1_ENABLED
    module function getComplementSortedCusCom_D0_SK1(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in)                :: setA
        character(*,SKC)        , intent(in)                :: setB
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        character(:,SKC)                    , allocatable   :: complement
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    module function getComplementSortedCusCom_D1_SK5(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

#if SK4_ENABLED
    module function getComplementSortedCusCom_D1_SK4(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

#if SK3_ENABLED
    module function getComplementSortedCusCom_D1_SK3(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

#if SK2_ENABLED
    module function getComplementSortedCusCom_D1_SK2(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

#if SK1_ENABLED
    module function getComplementSortedCusCom_D1_SK1(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in), contiguous    :: setA(:)
        character(*,SKC)        , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        character(LEN_STR,SKC)              , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    module function getComplementSortedCusCom_D1_IK5(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK4_ENABLED
    module function getComplementSortedCusCom_D1_IK4(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK3_ENABLED
    module function getComplementSortedCusCom_D1_IK3(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK2_ENABLED
    module function getComplementSortedCusCom_D1_IK2(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK1_ENABLED
    module function getComplementSortedCusCom_D1_IK1(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    module function getComplementSortedCusCom_D1_LK5(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

#if LK4_ENABLED
    module function getComplementSortedCusCom_D1_LK4(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

#if LK3_ENABLED
    module function getComplementSortedCusCom_D1_LK3(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

#if LK2_ENABLED
    module function getComplementSortedCusCom_D1_LK2(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

#if LK1_ENABLED
    module function getComplementSortedCusCom_D1_LK1(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)            , intent(in), contiguous    :: setA(:)
        logical(LKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        logical(LKC)                        , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module function getComplementSortedCusCom_D1_CK5(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

#if CK4_ENABLED
    module function getComplementSortedCusCom_D1_CK4(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

#if CK3_ENABLED
    module function getComplementSortedCusCom_D1_CK3(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

#if CK2_ENABLED
    module function getComplementSortedCusCom_D1_CK2(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

#if CK1_ENABLED
    module function getComplementSortedCusCom_D1_CK1(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in), contiguous    :: setA(:)
        complex(CKC)            , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        complex(CKC)                        , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module function getComplementSortedCusCom_D1_RK5(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

#if RK4_ENABLED
    module function getComplementSortedCusCom_D1_RK4(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

#if RK3_ENABLED
    module function getComplementSortedCusCom_D1_RK3(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

#if RK2_ENABLED
    module function getComplementSortedCusCom_D1_RK2(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

#if RK1_ENABLED
    module function getComplementSortedCusCom_D1_RK1(setA, setB, sorted, unique, iseq) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementSortedCusCom_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in), contiguous    :: setA(:)
        real(RKC)               , intent(in), contiguous    :: setB(:)
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        procedure(logical(LK))                              :: iseq
        real(RKC)                           , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the complement of the input set `setA` with respect to an array of elements of
    !>  the sorted set comprised of a range of integer values constructed by the triple `start, stop, step`.
    !>
    !>  \details
    !>  This procedures under this generic interface are particularly useful for a fast construction of the complement to a set of array indices.<br>
    !>  The additional performance originates from avoiding the unnecessary creation of `setB = getRange(start, stop, step)` by the user as the input
    !>  argument to [getComplement](@ref pm_arrayComplement::getComplement).
    !>
    !>  \param[in]  setA        :   The input `contiguous` vector of type `integer` of kind \IKALL whose complement in
    !>                              the range specified by the triple `(start, stop, step)` will be returned as `complement(:)`.
    !>  \param[in]  start       :   The input scalar of the same type, kind, and rank as `setA` representing the start of the range.
    !>  \param[in]  stop        :   The input scalar of the same type, kind, and rank as `setA` representing the stop of the range.
    !>  \param[in]  step        :   The input scalar of the same type, kind, and rank as `setA` representing the jumping step size in the range.
    !>  \param[in]  sorted      :   The input scalar `logical` of default kind \LK.
    !>                              <ol>
    !>                                  <li> If `.false.`, the input `setA` is assumed to be dissimilarly sorted or not sorted at all with respect to the range specified by `(start, stop, step)`.
    !>                                  <li> If `.true.`, the input `setA` is assumed to be **sorted similar to** the range specified by `(start, stop, step)` (e.g., both ascending or both descending order).
    !>                              </ol>
    !>                              If the input sets are similarly sorted, then specifying `sorted = .true._LK` can lead to significantly better runtime performance.<br>
    !>                              (**optional**, default = `.false._LK`. It must be present <b>if and only if</b> `unique` is also present.)
    !>  \param[in]  unique      :   The input scalar `logical` of default kind \LK.
    !>                              <ol>
    !>                                  <li> If `.false.`, the input `setA` is assumed to possibly contain duplicate elements.
    !>                                  <li> If `.true.`, all elements in the input `setA` are assumed to be unique.
    !>                              </ol>
    !>                              If the elements of `setA` are unique, then specifying `unique = .true._LK` can lead to significantly better runtime performance.<br>
    !>                              The specified value for `unique` becomes relevant only if `sorted = .true.`. Its value is ignored when `sorted = .false.`<br>
    !>                              (**optional**, default = `.false._LK`. It must be present <b>if and only if</b> `sorted` is also present.)
    !>
    !>  \return
    !>  `complement`            :   The output `allocatable` object of the same type, kind, and rank as `setA` containing the complement of `setA`
    !>                              in the range specified by the triple `(start, stop, step)` (i.e., the values in the range that are not in `setA`).
    !>
    !>  \interface{getComplementRange}
    !>  \code{.F90}
    !>
    !>      use pm_arrayComplement, only: getComplementRange
    !>
    !>      complement(:) = getComplementRange(setA(:), start, stop, step)
    !>      complement(:) = getComplementRange(setA(:), start, stop, step, sorted, unique)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Beware of the possibility of arithmetic overflow and underflow when the sum or the difference of `start` and `stop` falls out of the allowed representable range.<br>
    !>
    !>  \warning
    !>  The input `step` must be non-zero.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \note
    !>  The performane of the procedures under this generic interface can dramatically improve if the input sets contain unique values (`unique = .true._LK`) and both sets
    !>  are similarly-sorted (`sorted = .true._LK`). The unique elements of an arbitrary set can be obtained via the procedures of [pm_arrayUnique](@ref pm_arrayUnique).
    !>
    !>  \see
    !>  [getComplementRange](@ref pm_arrayComplement::getComplement)<br>
    !>
    !>  \example{getComplementRange}
    !>  \include{lineno} example/pm_arrayComplement/getComplementRange/main.F90
    !>  \compilef{getComplementRange}
    !>  \output{getComplementRange}
    !>  \include{lineno} example/pm_arrayComplement/getComplementRange/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayComplement](@ref test_pm_arrayComplement)
    !>
    !>  \finmain{getComplementRange}
    !>
    !>  \author
    !>  \FatemehBagheri, Wednesday 1:35 PM, August 11, 2021, Dallas, TX
    interface getComplementRange

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getCompRangeRandom_D1_IK5(setA, start, stop, step) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompRangeRandom_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in)                :: start, stop, step
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getCompRangeRandom_D1_IK4(setA, start, stop, step) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompRangeRandom_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in)                :: start, stop, step
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getCompRangeRandom_D1_IK3(setA, start, stop, step) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompRangeRandom_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in)                :: start, stop, step
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getCompRangeRandom_D1_IK2(setA, start, stop, step) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompRangeRandom_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in)                :: start, stop, step
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK1_ENABLED
    PURE module function getCompRangeRandom_D1_IK1(setA, start, stop, step) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCompRangeRandom_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in)                :: start, stop, step
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getCompRangeSorted_D1_IK5(setA, start, stop, step, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRangeSorted_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in)                :: start, stop, step
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getCompRangeSorted_D1_IK4(setA, start, stop, step, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRangeSorted_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in)                :: start, stop, step
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getCompRangeSorted_D1_IK3(setA, start, stop, step, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRangeSorted_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in)                :: start, stop, step
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getCompRangeSorted_D1_IK2(setA, start, stop, step, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRangeSorted_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in)                :: start, stop, step
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

#if IK1_ENABLED
    PURE module function getCompRangeSorted_D1_IK1(setA, start, stop, step, sorted, unique) result(complement)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getComplementRangeSorted_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in), contiguous    :: setA(:)
        integer(IKC)            , intent(in)                :: start, stop, step
        logical(LK)             , intent(in)                :: sorted
        logical(LK)             , intent(in)                :: unique
        integer(IKC)                        , allocatable   :: complement(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

end module pm_arrayComplement
