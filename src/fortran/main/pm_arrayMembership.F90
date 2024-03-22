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
!>  This module contains procedures and generic interfaces for assessing whether
!>  **particular** value(s) or **any** values or **all** values within a collection are members of
!>  another **collection** of values, or within a **range** of values that specifies a **mathematical set**.<br>
!>
!>  \details
!>  In mathematics, an **element** or **member** of a set is any one of the distinct objects that belong to that set.<br>
!>  The relation *is an element of*, also called **set membership**, is denoted by the symbol \f$\in\f$.<br>
!>  Writing \f$x\in A\f$ means that \f$x\f$ is an element of \f$A\f$.<br>
!>  Equivalent expressions are \f$x\f$ is a member of \f$A\f$, \f$x\f$ belongs to \f$A\f$, \f$x\f$ is in \f$A\f$ and \f$x\f$ lies in \f$A\f$.<br>
!>  The expressions \f$A\f$ includes \f$x\f$ and \f$A\f$ contains \f$x\f$ are also used to mean set membership,
!>  although it is more frequently used to mean instead \f$x\f$ is a subset of \f$A\f$.<br>
!>  Logician George Boolos strongly urged that *contains* be used for *membership* only, and *includes* for the subset relation only.<br>
!>  For the relation \f$\in\f$, the converse relation \f$\in^T\f$ may be written \f$A\ni x\f$ meaning \f$A\f$ contains or includes \f$x\f$.<br>
!>  The **negation** of set membership is denoted by the symbol \f$\notin\f$.<br>
!>  Writing \f$x\notin A\f$ means that \f$x\f$ is not an element of \f$A\f$.<br>
!>
!>  \note
!>  The functionalities of the generic interfaces [.allin.](@ref pm_arrayMembership_allin) and [.allinrange.](@ref pm_arrayMembership_allinrange)
!>  of this module are similar to that of the set-theoretic **subset** operation \f$\subset\f$.<br>
!>  In mathematics, set \f$A\f$ is a subset of a set \f$B\f$ if all elements of \f$A\f$ are also elements of \f$B\f$.<br>
!>  Then, \f$B\f$ is then a **superset** of \f$A\f$.<br>
!>  It is possible for \f$A\f$ and \f$B\f$ to be equal.<br>
!>  If they are unequal, then \f$A\f$ is a **proper subset** of \f$B\f$.<br>
!>  The relationship of one set being a subset of another is called **inclusion** (or sometimes **containment**).<br>
!>  The statement \f$A\f$ is a subset of \f$B\f$ may also be expressed as \f$B\f$ includes (or contains) \f$A\f$ or \f$A\f$ is included (or contained) in \f$B\f$.<br>
!>  A \f$k\f$-subset is a subset with \f$k\f$ elements.<br>
!>
!>  \note
!>  The functionalities of the generic interfaces [.anyin.](@ref pm_arrayMembership_anyin) and [.anyinrange.](@ref pm_arrayMembership_anyinrange)
!>  of this module are similar to testing whether the output of set-theoretic **intersection** operation \f$A\cap B\f$ is non-empty.<br>
!>  In set theory, the intersection of two sets \f$A\f$ and \f$B\f$, denoted by \f$A\cap B\f$ is the set containing any elements
!>  of \f$A\f$ that also belong to \f$B\f$ or equivalently, any elements of \f$B\f$ that also belong to \f$A\f$.<br>
!>
!>  \see
!>  [operator(.in.)](@ref pm_arrayMembership_in)<br>
!>  [operator(.allin.)](@ref pm_arrayMembership_allin)<br>
!>  [operator(.anyin.)](@ref pm_arrayMembership_anyin)<br>
!>  [operator(.inrange.)](@ref pm_arrayMembership_inrange)<br>
!>  [operator(.anyinrange.)](@ref pm_arrayMembership_anyinrange)<br>
!>  [operator(.allinrange.)](@ref pm_arrayMembership_allinrange)<br>
!>  [operator(.divmul.)](@ref pm_mathDivMul_divmul)<br>
!>  [operator(.subadd.)](@ref pm_mathSubAdd_subadd)<br>
!>  [getMinMax](@ref pm_mathMinMax::getMinMax)<br>
!>  [setMinMax](@ref pm_mathMinMax::setMinMax)<br>
!>  [pm_arraySearch](@ref pm_arraySearch)<br>
!>  [pm_arrayFind](@ref pm_arrayFind)<br>
!>
!>  \test
!>  [test_pm_arrayMembership](@ref test_pm_arrayMembership)
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Wednesday 12:24 PM, August 11, 2021, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_arrayMembership

    use pm_kind, only: SK, IK, LK
    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_arrayMembership"

    !>  \brief
    !>  \anchor pm_arrayMembership_in
    !>  Generate and return `.true.` if the input value `val` is a member of the input array-like object `set`, otherwise, return `.false.`.
    !>
    !>  \param[in]  val :   The input scalar or vector of,
    !>                      <ul>
    !>                          <li> type `character` of kind \SKALL or arbitrary length type parameter, or
    !>                          <li> type `integer` of kind \IKALL, or
    !>                          <li> type `logical` of kind \LKALL, or
    !>                          <li> type `complex` of kind \CKALL, or
    !>                          <li> type `real` of kind \RKALL, or
    !>                      </ul>
    !>                      whose (elemental) membership(s) in the `set` will be checked.
    !>  \param[in]  set :   The input scalar of the same `character` type and kind as `val` or vector of the same type and kind as `val`,
    !>                      representing the set with respect to which the (elemental) membership of `val` will be checked.
    !>
    !>  \return
    !>  `member`        :   The output scalar or vector of the same size as the input `val`, of type `logical` of default kind \LK.<br>
    !>                      It is `.true.` <b>if and only if</b> the corresponding `val` is a member of the input `set`.<br>
    !>
    !>  \interface{in}
    !>  \code{.F90}
    !>
    !>      use pm_arrayMembership, only: operator(.in.)
    !>
    !>      member = val .in. set ! scalar string arguments, but logical vector output of size `len(val)`.
    !>      member = val .in. set(:)
    !>      member(1:size(val)) = val(:) .in. set(:)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \see
    !>  [operator(.in.)](@ref pm_arrayMembership_in)<br>
    !>  [operator(.allin.)](@ref pm_arrayMembership_allin)<br>
    !>  [operator(.anyin.)](@ref pm_arrayMembership_anyin)<br>
    !>  [operator(.inrange.)](@ref pm_arrayMembership_inrange)<br>
    !>  [operator(.anyinrange.)](@ref pm_arrayMembership_anyinrange)<br>
    !>  [operator(.allinrange.)](@ref pm_arrayMembership_allinrange)<br>
    !>  [operator(.divmul.)](@ref pm_mathDivMul_divmul)<br>
    !>  [operator(.subadd.)](@ref pm_mathSubAdd_subadd)<br>
    !>  [getMinMax](@ref pm_mathMinMax::getMinMax)<br>
    !>  [setMinMax](@ref pm_mathMinMax::setMinMax)<br>
    !>  [pm_arraySearch](@ref pm_arraySearch)<br>
    !>  [pm_arrayFind](@ref pm_arrayFind)<br>
    !>
    !>  \example{in}
    !>  \include{lineno} example/pm_arrayMembership/in/main.F90
    !>  \compilef{in}
    !>  \output{in}
    !>  \include{lineno} example/pm_arrayMembership/in/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayMembership](@ref test_pm_arrayMembership)
    !>
    !>  \finmain{in}
    !>
    !>  \author
    !>  \FatemehBagheri, Wednesday 1:35 PM, August 11, 2021, Dallas, TX
    interface operator(.in.)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function in_D0_D0_SK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: set
        logical(LK)                                         :: member(len(val, IK))
    end function
#endif

#if SK4_ENABLED
    pure module function in_D0_D0_SK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: set
        logical(LK)                                         :: member(len(val, IK))
    end function
#endif

#if SK3_ENABLED
    pure module function in_D0_D0_SK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: set
        logical(LK)                                         :: member(len(val, IK))
    end function
#endif

#if SK2_ENABLED
    pure module function in_D0_D0_SK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: set
        logical(LK)                                         :: member(len(val, IK))
    end function
#endif

#if SK1_ENABLED
    pure module function in_D0_D0_SK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: set
        logical(LK)                                         :: member(len(val, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function in_D0_D1_SK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if SK4_ENABLED
    pure module function in_D0_D1_SK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if SK3_ENABLED
    pure module function in_D0_D1_SK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if SK2_ENABLED
    pure module function in_D0_D1_SK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if SK1_ENABLED
    pure module function in_D0_D1_SK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function in_D0_D1_IK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)                :: val
        integer(IKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if IK4_ENABLED
    pure module function in_D0_D1_IK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)                :: val
        integer(IKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if IK3_ENABLED
    pure module function in_D0_D1_IK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)                :: val
        integer(IKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if IK2_ENABLED
    pure module function in_D0_D1_IK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)                :: val
        integer(IKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if IK1_ENABLED
    pure module function in_D0_D1_IK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)                :: val
        integer(IKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module function in_D0_D1_LK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)            , intent(in)                :: val
        logical(LKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if LK4_ENABLED
    pure module function in_D0_D1_LK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)            , intent(in)                :: val
        logical(LKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if LK3_ENABLED
    pure module function in_D0_D1_LK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)            , intent(in)                :: val
        logical(LKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if LK2_ENABLED
    pure module function in_D0_D1_LK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)            , intent(in)                :: val
        logical(LKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if LK1_ENABLED
    pure module function in_D0_D1_LK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)            , intent(in)                :: val
        logical(LKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function in_D0_D1_CK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)                :: val
        complex(CKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if CK4_ENABLED
    pure module function in_D0_D1_CK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)                :: val
        complex(CKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if CK3_ENABLED
    pure module function in_D0_D1_CK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)                :: val
        complex(CKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if CK2_ENABLED
    pure module function in_D0_D1_CK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)                :: val
        complex(CKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if CK1_ENABLED
    pure module function in_D0_D1_CK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)                :: val
        complex(CKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function in_D0_D1_RK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)                :: val
        real(RKC)               , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if RK4_ENABLED
    pure module function in_D0_D1_RK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)                :: val
        real(RKC)               , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if RK3_ENABLED
    pure module function in_D0_D1_RK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)                :: val
        real(RKC)               , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if RK2_ENABLED
    pure module function in_D0_D1_RK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)                :: val
        real(RKC)               , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

#if RK1_ENABLED
    pure module function in_D0_D1_RK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D0_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)                :: val
        real(RKC)               , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function in_D1_D1_SK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if SK4_ENABLED
    pure module function in_D1_D1_SK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if SK3_ENABLED
    pure module function in_D1_D1_SK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if SK2_ENABLED
    pure module function in_D1_D1_SK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if SK1_ENABLED
    pure module function in_D1_D1_SK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function in_D1_D1_IK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if IK4_ENABLED
    pure module function in_D1_D1_IK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if IK3_ENABLED
    pure module function in_D1_D1_IK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if IK2_ENABLED
    pure module function in_D1_D1_IK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if IK1_ENABLED
    pure module function in_D1_D1_IK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module function in_D1_D1_LK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if LK4_ENABLED
    pure module function in_D1_D1_LK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if LK3_ENABLED
    pure module function in_D1_D1_LK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if LK2_ENABLED
    pure module function in_D1_D1_LK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if LK1_ENABLED
    pure module function in_D1_D1_LK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function in_D1_D1_CK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if CK4_ENABLED
    pure module function in_D1_D1_CK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if CK3_ENABLED
    pure module function in_D1_D1_CK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if CK2_ENABLED
    pure module function in_D1_D1_CK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if CK1_ENABLED
    pure module function in_D1_D1_CK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function in_D1_D1_RK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if RK4_ENABLED
    pure module function in_D1_D1_RK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if RK3_ENABLED
    pure module function in_D1_D1_RK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if RK2_ENABLED
    pure module function in_D1_D1_RK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if RK1_ENABLED
    pure module function in_D1_D1_RK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: in_D1_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in), contiguous    :: set(:)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_arrayMembership_inrange
    !>  Generate and return `.true.` if the input value `val` is within a range specified by the input array-like object `set(1:2)`, otherwise, return `.false.`.
    !>
    !>  \param[in]  val :   The input scalar or vector of,
    !>                      <ul>
    !>                          <li> type `character` of kind \SKALL or arbitrary length type parameter, or
    !>                          <li> type `integer` of kind \IKALL, or
    !>                          <li> type `logical` of kind \LKALL, or
    !>                          <li> type `complex` of kind \CKALL, or
    !>                          <li> type `real` of kind \RKALL, or
    !>                      </ul>
    !>                      whose (elemental) membership(s) in the `set` will be checked.
    !>  \param[in]  set :   The input scalar of the same `character` type and kind as `val` of length-type-parameter `2` or vector of the same type and kind as `val` of size `2`,
    !>                      representing the range `[set(1) : set(2)]` with respect to which the (elemental) membership of `val` will be checked.<br>
    !>                      By convention,
    !>                      <ol>
    !>                          <li>    If the input arguments are of type `character`, `integer`, `real`, then the values are compared as defined by the standard.
    !>                          <li>    If the input arguments are of type `complex`, then the input values are elementally compared [lexically](@ref pm_arrayCompareLex) (similar to string comparison).
    !>                          <li>    If the input arguments are of type `logical`, then `.false. < .true.` is assumed.
    !>                      </ol>
    !>
    !>  \return
    !>  `member`        :   The output scalar or vector of the same size as the input `val`, of type `logical` of default kind \LK.<br>
    !>                      It is `.true.` <b>if and only if</b> the corresponding `val` is within the input range specified by `set`.<br>
    !>
    !>  \interface{inrange}
    !>  \code{.F90}
    !>
    !>      use pm_arrayMembership, only: operator(.inrange.)
    !>
    !>      member = val .inrange. set(1:2)
    !>      member(1:size(val)) = val(:) .inrange. set(1:2)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \note
    !>  Normally, the condition `set(1) < set(2)` is expected to hold for the input `set` argument, otherwise, the specified range corresponds to an empty set.
    !>
    !>  \note
    !>  This generic interface offers a neat way of checking a value against a range returned by either
    !>  [getMinMax()](@ref pm_mathMinMax::getMinMax) or [operator(.subadd.)](@ref pm_mathSubAdd_subadd) or [operator(.divmul.)](@ref pm_mathDivMul_divmul).
    !>
    !>  \see
    !>  [operator(.in.)](@ref pm_arrayMembership_in)<br>
    !>  [operator(.allin.)](@ref pm_arrayMembership_allin)<br>
    !>  [operator(.anyin.)](@ref pm_arrayMembership_anyin)<br>
    !>  [operator(.inrange.)](@ref pm_arrayMembership_inrange)<br>
    !>  [operator(.anyinrange.)](@ref pm_arrayMembership_anyinrange)<br>
    !>  [operator(.allinrange.)](@ref pm_arrayMembership_allinrange)<br>
    !>  [operator(.divmul.)](@ref pm_mathDivMul_divmul)<br>
    !>  [operator(.subadd.)](@ref pm_mathSubAdd_subadd)<br>
    !>  [getMinMax](@ref pm_mathMinMax::getMinMax)<br>
    !>  [setMinMax](@ref pm_mathMinMax::setMinMax)<br>
    !>  [pm_arraySearch](@ref pm_arraySearch)<br>
    !>  [pm_arrayFind](@ref pm_arrayFind)<br>
    !>
    !>  \example{inrange}
    !>  \include{lineno} example/pm_arrayMembership/inrange/main.F90
    !>  \compilef{inrange}
    !>  \output{inrange}
    !>  \include{lineno} example/pm_arrayMembership/inrange/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayMembership](@ref test_pm_arrayMembership)
    !>
    !>  \finmain{inrange}
    !>
    !>  \author
    !>  \FatemehBagheri, Wednesday 1:35 PM, August 11, 2021, Dallas, TX
    interface operator(.inrange.)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function inrange_D0_D0_SK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in)                :: val
        character(2,SKC)        , intent(in)                :: set
        logical(LK)                                         :: member(len(val, IK))
    end function
#endif

#if SK4_ENABLED
    pure module function inrange_D0_D0_SK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in)                :: val
        character(2,SKC)        , intent(in)                :: set
        logical(LK)                                         :: member(len(val, IK))
    end function
#endif

#if SK3_ENABLED
    pure module function inrange_D0_D0_SK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in)                :: val
        character(2,SKC)        , intent(in)                :: set
        logical(LK)                                         :: member(len(val, IK))
    end function
#endif

#if SK2_ENABLED
    pure module function inrange_D0_D0_SK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in)                :: val
        character(2,SKC)        , intent(in)                :: set
        logical(LK)                                         :: member(len(val, IK))
    end function
#endif

#if SK1_ENABLED
    pure module function inrange_D0_D0_SK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in)                :: val
        character(2,SKC)        , intent(in)                :: set
        logical(LK)                                         :: member(len(val, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function inrange_D0_D1_SK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if SK4_ENABLED
    pure module function inrange_D0_D1_SK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if SK3_ENABLED
    pure module function inrange_D0_D1_SK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if SK2_ENABLED
    pure module function inrange_D0_D1_SK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if SK1_ENABLED
    pure module function inrange_D0_D1_SK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function inrange_D0_D1_IK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)                :: val
        integer(IKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if IK4_ENABLED
    pure module function inrange_D0_D1_IK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)                :: val
        integer(IKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if IK3_ENABLED
    pure module function inrange_D0_D1_IK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)                :: val
        integer(IKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if IK2_ENABLED
    pure module function inrange_D0_D1_IK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)                :: val
        integer(IKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if IK1_ENABLED
    pure module function inrange_D0_D1_IK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)                :: val
        integer(IKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module function inrange_D0_D1_LK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)            , intent(in)                :: val
        logical(LKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if LK4_ENABLED
    pure module function inrange_D0_D1_LK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)            , intent(in)                :: val
        logical(LKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if LK3_ENABLED
    pure module function inrange_D0_D1_LK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)            , intent(in)                :: val
        logical(LKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if LK2_ENABLED
    pure module function inrange_D0_D1_LK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)            , intent(in)                :: val
        logical(LKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if LK1_ENABLED
    pure module function inrange_D0_D1_LK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)            , intent(in)                :: val
        logical(LKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function inrange_D0_D1_CK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)                :: val
        complex(CKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if CK4_ENABLED
    pure module function inrange_D0_D1_CK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)                :: val
        complex(CKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if CK3_ENABLED
    pure module function inrange_D0_D1_CK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)                :: val
        complex(CKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if CK2_ENABLED
    pure module function inrange_D0_D1_CK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)                :: val
        complex(CKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if CK1_ENABLED
    pure module function inrange_D0_D1_CK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)                :: val
        complex(CKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function inrange_D0_D1_RK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)                :: val
        real(RKC)               , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if RK4_ENABLED
    pure module function inrange_D0_D1_RK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)                :: val
        real(RKC)               , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if RK3_ENABLED
    pure module function inrange_D0_D1_RK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)                :: val
        real(RKC)               , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if RK2_ENABLED
    pure module function inrange_D0_D1_RK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)                :: val
        real(RKC)               , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

#if RK1_ENABLED
    pure module function inrange_D0_D1_RK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D0_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)                :: val
        real(RKC)               , intent(in)                :: set(2)
        logical(LK)                                         :: member
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function inrange_D1_D1_SK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if SK4_ENABLED
    pure module function inrange_D1_D1_SK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if SK3_ENABLED
    pure module function inrange_D1_D1_SK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if SK2_ENABLED
    pure module function inrange_D1_D1_SK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if SK1_ENABLED
    pure module function inrange_D1_D1_SK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function inrange_D1_D1_IK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if IK4_ENABLED
    pure module function inrange_D1_D1_IK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if IK3_ENABLED
    pure module function inrange_D1_D1_IK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if IK2_ENABLED
    pure module function inrange_D1_D1_IK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if IK1_ENABLED
    pure module function inrange_D1_D1_IK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module function inrange_D1_D1_LK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if LK4_ENABLED
    pure module function inrange_D1_D1_LK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if LK3_ENABLED
    pure module function inrange_D1_D1_LK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if LK2_ENABLED
    pure module function inrange_D1_D1_LK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if LK1_ENABLED
    pure module function inrange_D1_D1_LK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function inrange_D1_D1_CK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if CK4_ENABLED
    pure module function inrange_D1_D1_CK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if CK3_ENABLED
    pure module function inrange_D1_D1_CK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if CK2_ENABLED
    pure module function inrange_D1_D1_CK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if CK1_ENABLED
    pure module function inrange_D1_D1_CK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function inrange_D1_D1_RK5(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if RK4_ENABLED
    pure module function inrange_D1_D1_RK4(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if RK3_ENABLED
    pure module function inrange_D1_D1_RK3(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if RK2_ENABLED
    pure module function inrange_D1_D1_RK2(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

#if RK1_ENABLED
    pure module function inrange_D1_D1_RK1(val, set) result(member)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: inrange_D1_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in)                :: set(2)
        logical(LK)                                         :: member(size(val, 1, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_arrayMembership_allin
    !>  Generate and return `.true.` if **all elements** of the input array-like `val` are members of the input array-like object `Set`, otherwise, return `.false.`.
    !>
    !>  \details
    !>  The functionality of this interface is equivalent to the result of the set-theoretic operation \f$\ms{val}\subset\ms{Set}\f$.
    !>
    !>  \param[in]  val :   The input scalar of
    !>                      <ul>
    !>                          <li> type `character` of kind \SKALL or arbitrary length type parameter,
    !>                      </ul>
    !>                      or vector of
    !>                      <ul>
    !>                          <li> type `character` of kind \SKALL or arbitrary length type parameter, or
    !>                          <li> type `integer` of kind \IKALL, or
    !>                          <li> type `logical` of kind \LKALL, or
    !>                          <li> type `complex` of kind \CKALL, or
    !>                          <li> type `real` of kind \RKALL, or
    !>                      </ul>
    !>                      whose elements memberships in the `Set` will be checked.
    !>  \param[in]  Set :   The input object of the same type, kind, and rank as the input `val`,
    !>                      representing the set with respect to which **all elements** of `val` will be checked.
    !>
    !>  \return
    !>  `allMember`     :   The output scalar `logical` of default kind \LK.<br>
    !>                      It is `.true.` <b>if and only if all elements</b> of the input `val` are members of the input `Set`.<br>
    !>
    !>  \interface{allin}
    !>  \code{.F90}
    !>
    !>      use pm_arrayMembership, only: operator(.allin.)
    !>
    !>      allMember = val .allin. Set ! scalar character arguments
    !>      allMember = val(:) .allin. Set(:)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \see
    !>  [operator(.in.)](@ref pm_arrayMembership_in)<br>
    !>  [operator(.allin.)](@ref pm_arrayMembership_allin)<br>
    !>  [operator(.anyin.)](@ref pm_arrayMembership_anyin)<br>
    !>  [operator(.inrange.)](@ref pm_arrayMembership_inrange)<br>
    !>  [operator(.anyinrange.)](@ref pm_arrayMembership_anyinrange)<br>
    !>  [operator(.allinrange.)](@ref pm_arrayMembership_allinrange)<br>
    !>  [operator(.divmul.)](@ref pm_mathDivMul_divmul)<br>
    !>  [operator(.subadd.)](@ref pm_mathSubAdd_subadd)<br>
    !>  [getMinMax](@ref pm_mathMinMax::getMinMax)<br>
    !>  [setMinMax](@ref pm_mathMinMax::setMinMax)<br>
    !>  [pm_arraySearch](@ref pm_arraySearch)<br>
    !>  [pm_arrayFind](@ref pm_arrayFind)<br>
    !>
    !>  \example{allin}
    !>  \include{lineno} example/pm_arrayMembership/allin/main.F90
    !>  \compilef{allin}
    !>  \output{allin}
    !>  \include{lineno} example/pm_arrayMembership/allin/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayMembership](@ref test_pm_arrayMembership)
    !>
    !>  \finmain{allin}
    !>
    !>  \author
    !>  \FatemehBagheri, Wednesday 1:35 PM, August 11, 2021, Dallas, TX
    interface operator(.allin.)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function allin_D0_D0_SK5(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D0_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: allMember
    end function
#endif

#if SK4_ENABLED
    pure module function allin_D0_D0_SK4(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D0_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: allMember
    end function
#endif

#if SK3_ENABLED
    pure module function allin_D0_D0_SK3(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D0_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: allMember
    end function
#endif

#if SK2_ENABLED
    pure module function allin_D0_D0_SK2(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D0_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: allMember
    end function
#endif

#if SK1_ENABLED
    pure module function allin_D0_D0_SK1(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D0_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: allMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function allin_D1_D1_SK5(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if SK4_ENABLED
    pure module function allin_D1_D1_SK4(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if SK3_ENABLED
    pure module function allin_D1_D1_SK3(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if SK2_ENABLED
    pure module function allin_D1_D1_SK2(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if SK1_ENABLED
    pure module function allin_D1_D1_SK1(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function allin_D1_D1_IK5(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if IK4_ENABLED
    pure module function allin_D1_D1_IK4(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if IK3_ENABLED
    pure module function allin_D1_D1_IK3(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if IK2_ENABLED
    pure module function allin_D1_D1_IK2(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if IK1_ENABLED
    pure module function allin_D1_D1_IK1(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module function allin_D1_D1_LK5(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if LK4_ENABLED
    pure module function allin_D1_D1_LK4(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if LK3_ENABLED
    pure module function allin_D1_D1_LK3(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if LK2_ENABLED
    pure module function allin_D1_D1_LK2(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if LK1_ENABLED
    pure module function allin_D1_D1_LK1(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function allin_D1_D1_CK5(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if CK4_ENABLED
    pure module function allin_D1_D1_CK4(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if CK3_ENABLED
    pure module function allin_D1_D1_CK3(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if CK2_ENABLED
    pure module function allin_D1_D1_CK2(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if CK1_ENABLED
    pure module function allin_D1_D1_CK1(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function allin_D1_D1_RK5(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if RK4_ENABLED
    pure module function allin_D1_D1_RK4(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if RK3_ENABLED
    pure module function allin_D1_D1_RK3(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if RK2_ENABLED
    pure module function allin_D1_D1_RK2(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

#if RK1_ENABLED
    pure module function allin_D1_D1_RK1(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allin_D1_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: allMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_arrayMembership_allinrange
    !>  Generate and return `.true.` if **all elements** of the input array-like object `val` are within a range specified by the input vector `Set(1:2)`, otherwise, return `.false.`.
    !>
    !>  \details
    !>  The functionality of this interface is equivalent to the result of the set-theoretic operation \f$\ms{val}\subset\ms{Set}\f$.<br>
    !>
    !>  \param[in]  val :   The input scalar of
    !>                      <ul>
    !>                          <li> type `character` of kind \SKALL or arbitrary length type parameter,
    !>                      </ul>
    !>                      or vector of
    !>                      <ul>
    !>                          <li> type `character` of kind \SKALL or arbitrary length type parameter, or
    !>                          <li> type `integer` of kind \IKALL, or
    !>                          <li> type `logical` of kind \LKALL, or
    !>                          <li> type `complex` of kind \CKALL, or
    !>                          <li> type `real` of kind \RKALL, or
    !>                      </ul>
    !>                      whose elements memberships in the `Set` will be checked.
    !>  \param[in]  Set :   The input scalar of the same `character` type and kind as `val` of length-type-parameter `2` or vector of the same type and kind as `val` of size `2`,
    !>                      representing the range `[Set(1) : Set(2)]` with respect to which the (elemental) membership of `val` will be checked.<br>
    !>                      By convention,
    !>                      <ol>
    !>                          <li>    If the input arguments are of type `character`, `integer`, `real`, then the values are compared as defined by the standard.
    !>                          <li>    If the input arguments are of type `complex`, then the input values are elementally compared [lexically](@ref pm_arrayCompareLex) (similar to string comparison).
    !>                          <li>    If the input arguments are of type `logical`, then `.false. < .true.` is assumed.
    !>                      </ol>
    !>
    !>  \return
    !>  `allMember`     :   The output scalar or vector of the same size as the input `val`, of type `logical` of default kind \LK.<br>
    !>                      It is `.true.` <b>if and only if all elements</b> of the input `val` are within the input range specified by `Set`.<br>
    !>
    !>  \interface{allinrange}
    !>  \code{.F90}
    !>
    !>      use pm_arrayMembership, only: operator(.allinrange.)
    !>
    !>      allMember = val .allinrange. Set(1:2) ! scalar character.
    !>      allMember = val(:) .allinrange. Set(1:2)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \note
    !>  Normally, the condition `Set(1) < Set(2)` is expected to hold for the input `Set` argument, otherwise, the specified range corresponds to an empty set.
    !>
    !>  \note
    !>  This generic interface offers a neat way of checking a value against a range returned by either
    !>  [getMinMax()](@ref pm_mathMinMax::getMinMax) or [operator(.subadd.)](@ref pm_mathSubAdd_subadd) or [operator(.divmul.)](@ref pm_mathDivMul_divmul).<br>
    !>
    !>  \see
    !>  [operator(.in.)](@ref pm_arrayMembership_in)<br>
    !>  [operator(.allin.)](@ref pm_arrayMembership_allin)<br>
    !>  [operator(.anyin.)](@ref pm_arrayMembership_anyin)<br>
    !>  [operator(.inrange.)](@ref pm_arrayMembership_inrange)<br>
    !>  [operator(.anyinrange.)](@ref pm_arrayMembership_anyinrange)<br>
    !>  [operator(.allinrange.)](@ref pm_arrayMembership_allinrange)<br>
    !>  [operator(.divmul.)](@ref pm_mathDivMul_divmul)<br>
    !>  [operator(.subadd.)](@ref pm_mathSubAdd_subadd)<br>
    !>  [getMinMax](@ref pm_mathMinMax::getMinMax)<br>
    !>  [setMinMax](@ref pm_mathMinMax::setMinMax)<br>
    !>  [pm_arraySearch](@ref pm_arraySearch)<br>
    !>  [pm_arrayFind](@ref pm_arrayFind)<br>
    !>
    !>  \example{allinrange}
    !>  \include{lineno} example/pm_arrayMembership/allinrange/main.F90
    !>  \compilef{allinrange}
    !>  \output{allinrange}
    !>  \include{lineno} example/pm_arrayMembership/allinrange/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayMembership](@ref test_pm_arrayMembership)
    !>
    !>  \finmain{allinrange}
    !>
    !>  \author
    !>  \FatemehBagheri, Wednesday 1:35 PM, August 11, 2021, Dallas, TX
    interface operator(.allinrange.)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function allinrange_D0_D0_SK5(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D0_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in)                :: val
        character(2,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: allMember
    end function
#endif

#if SK4_ENABLED
    pure module function allinrange_D0_D0_SK4(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D0_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in)                :: val
        character(2,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: allMember
    end function
#endif

#if SK3_ENABLED
    pure module function allinrange_D0_D0_SK3(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D0_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in)                :: val
        character(2,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: allMember
    end function
#endif

#if SK2_ENABLED
    pure module function allinrange_D0_D0_SK2(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D0_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in)                :: val
        character(2,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: allMember
    end function
#endif

#if SK1_ENABLED
    pure module function allinrange_D0_D0_SK1(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D0_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in)                :: val
        character(2,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: allMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function allinrange_D1_D1_SK5(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if SK4_ENABLED
    pure module function allinrange_D1_D1_SK4(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if SK3_ENABLED
    pure module function allinrange_D1_D1_SK3(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if SK2_ENABLED
    pure module function allinrange_D1_D1_SK2(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if SK1_ENABLED
    pure module function allinrange_D1_D1_SK1(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function allinrange_D1_D1_IK5(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if IK4_ENABLED
    pure module function allinrange_D1_D1_IK4(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if IK3_ENABLED
    pure module function allinrange_D1_D1_IK3(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if IK2_ENABLED
    pure module function allinrange_D1_D1_IK2(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if IK1_ENABLED
    pure module function allinrange_D1_D1_IK1(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module function allinrange_D1_D1_LK5(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if LK4_ENABLED
    pure module function allinrange_D1_D1_LK4(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if LK3_ENABLED
    pure module function allinrange_D1_D1_LK3(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if LK2_ENABLED
    pure module function allinrange_D1_D1_LK2(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if LK1_ENABLED
    pure module function allinrange_D1_D1_LK1(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function allinrange_D1_D1_CK5(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if CK4_ENABLED
    pure module function allinrange_D1_D1_CK4(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if CK3_ENABLED
    pure module function allinrange_D1_D1_CK3(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if CK2_ENABLED
    pure module function allinrange_D1_D1_CK2(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if CK1_ENABLED
    pure module function allinrange_D1_D1_CK1(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function allinrange_D1_D1_RK5(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if RK4_ENABLED
    pure module function allinrange_D1_D1_RK4(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if RK3_ENABLED
    pure module function allinrange_D1_D1_RK3(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if RK2_ENABLED
    pure module function allinrange_D1_D1_RK2(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

#if RK1_ENABLED
    pure module function allinrange_D1_D1_RK1(val, Set) result(allMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: allinrange_D1_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in)                :: Set(2)
        logical(LK)                                         :: allMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_arrayMembership_anyin
    !>  Generate and return `.true.` if **any elements** of the input array-like `val` are members of the input array-like object `Set`, otherwise, return `.false.`.
    !>
    !>  \details
    !>  The functionality of this interface is equivalent to testing whether the output of set-theoretic **intersection** operation \f$\ms{val}\subset\ms{Set}\f$ is non-empty.<br>
    !>
    !>  \param[in]  val :   The input scalar of
    !>                      <ul>
    !>                          <li> type `character` of kind \SKALL or arbitrary length type parameter,
    !>                      </ul>
    !>                      or vector of
    !>                      <ul>
    !>                          <li> type `character` of kind \SKALL or arbitrary length type parameter, or
    !>                          <li> type `integer` of kind \IKALL, or
    !>                          <li> type `logical` of kind \LKALL, or
    !>                          <li> type `complex` of kind \CKALL, or
    !>                          <li> type `real` of kind \RKALL, or
    !>                      </ul>
    !>                      whose elements memberships in the `Set` will be checked.
    !>  \param[in]  Set :   The input object of the same type, kind, and rank as the input `val`,
    !>                      representing the set with respect to which **any elements** of `val` will be checked.
    !>
    !>  \return
    !>  `anyMember`     :   The output scalar `logical` of default kind \LK.<br>
    !>                      It is `.true.` <b>if and only if any elements</b> of the input `val` are members of the input `Set`.<br>
    !>
    !>  \interface{anyin}
    !>  \code{.F90}
    !>
    !>      use pm_arrayMembership, only: operator(.anyin.)
    !>
    !>      anyMember = val .anyin. Set ! scalar character arguments
    !>      anyMember = val(:) .anyin. Set(:)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \see
    !>  [operator(.in.)](@ref pm_arrayMembership_in)<br>
    !>  [operator(.allin.)](@ref pm_arrayMembership_allin)<br>
    !>  [operator(.anyin.)](@ref pm_arrayMembership_anyin)<br>
    !>  [operator(.inrange.)](@ref pm_arrayMembership_inrange)<br>
    !>  [operator(.anyinrange.)](@ref pm_arrayMembership_anyinrange)<br>
    !>  [operator(.allinrange.)](@ref pm_arrayMembership_allinrange)<br>
    !>  [operator(.divmul.)](@ref pm_mathDivMul_divmul)<br>
    !>  [operator(.subadd.)](@ref pm_mathSubAdd_subadd)<br>
    !>  [getMinMax](@ref pm_mathMinMax::getMinMax)<br>
    !>  [setMinMax](@ref pm_mathMinMax::setMinMax)<br>
    !>  [pm_arraySearch](@ref pm_arraySearch)<br>
    !>  [pm_arrayFind](@ref pm_arrayFind)<br>
    !>
    !>  \example{anyin}
    !>  \include{lineno} example/pm_arrayMembership/anyin/main.F90
    !>  \compilef{anyin}
    !>  \output{anyin}
    !>  \include{lineno} example/pm_arrayMembership/anyin/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayMembership](@ref test_pm_arrayMembership)
    !>
    !>  \finmain{anyin}
    !>
    !>  \author
    !>  \FatemehBagheri, Wednesday 1:35 PM, August 11, 2021, Dallas, TX
    interface operator(.anyin.)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function anyin_D0_D0_SK5(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D0_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: anyMember
    end function
#endif

#if SK4_ENABLED
    pure module function anyin_D0_D0_SK4(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D0_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: anyMember
    end function
#endif

#if SK3_ENABLED
    pure module function anyin_D0_D0_SK3(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D0_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: anyMember
    end function
#endif

#if SK2_ENABLED
    pure module function anyin_D0_D0_SK2(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D0_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: anyMember
    end function
#endif

#if SK1_ENABLED
    pure module function anyin_D0_D0_SK1(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D0_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in)                :: val
        character(*,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: anyMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function anyin_D1_D1_SK5(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if SK4_ENABLED
    pure module function anyin_D1_D1_SK4(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if SK3_ENABLED
    pure module function anyin_D1_D1_SK3(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if SK2_ENABLED
    pure module function anyin_D1_D1_SK2(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if SK1_ENABLED
    pure module function anyin_D1_D1_SK1(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function anyin_D1_D1_IK5(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if IK4_ENABLED
    pure module function anyin_D1_D1_IK4(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if IK3_ENABLED
    pure module function anyin_D1_D1_IK3(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if IK2_ENABLED
    pure module function anyin_D1_D1_IK2(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if IK1_ENABLED
    pure module function anyin_D1_D1_IK1(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module function anyin_D1_D1_LK5(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if LK4_ENABLED
    pure module function anyin_D1_D1_LK4(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if LK3_ENABLED
    pure module function anyin_D1_D1_LK3(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if LK2_ENABLED
    pure module function anyin_D1_D1_LK2(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if LK1_ENABLED
    pure module function anyin_D1_D1_LK1(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function anyin_D1_D1_CK5(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if CK4_ENABLED
    pure module function anyin_D1_D1_CK4(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if CK3_ENABLED
    pure module function anyin_D1_D1_CK3(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if CK2_ENABLED
    pure module function anyin_D1_D1_CK2(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if CK1_ENABLED
    pure module function anyin_D1_D1_CK1(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function anyin_D1_D1_RK5(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if RK4_ENABLED
    pure module function anyin_D1_D1_RK4(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if RK3_ENABLED
    pure module function anyin_D1_D1_RK3(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if RK2_ENABLED
    pure module function anyin_D1_D1_RK2(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

#if RK1_ENABLED
    pure module function anyin_D1_D1_RK1(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyin_D1_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in), contiguous    :: Set(:)
        logical(LK)                                         :: anyMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_arrayMembership_anyinrange
    !>  Generate and return `.true.` if **any elements** of the input array-like object `val` are within a range specified by the input vector `Set(1:2)`, otherwise, return `.false.`.
    !>
    !>  \details
    !>  The functionality of this interface is equivalent to testing whether the output of set-theoretic **intersection** operation \f$\ms{val}\subset\ms{Set}\f$ is non-empty.<br>
    !>
    !>  \param[in]  val :   The input scalar of
    !>                      <ul>
    !>                          <li> type `character` of kind \SKALL or arbitrary length type parameter,
    !>                      </ul>
    !>                      or vector of
    !>                      <ul>
    !>                          <li> type `character` of kind \SKALL or arbitrary length type parameter, or
    !>                          <li> type `integer` of kind \IKALL, or
    !>                          <li> type `logical` of kind \LKALL, or
    !>                          <li> type `complex` of kind \CKALL, or
    !>                          <li> type `real` of kind \RKALL, or
    !>                      </ul>
    !>                      whose elements memberships in the `Set` will be checked.
    !>  \param[in]  Set :   The input scalar of the same `character` type and kind as `val` of length-type-parameter `2` or vector of the same type and kind as `val` of size `2`,
    !>                      representing the range `[Set(1) : Set(2)]` with respect to which the (elemental) membership of `val` will be checked.<br>
    !>                      By convention,
    !>                      <ol>
    !>                          <li>    If the input arguments are of type `character`, `integer`, `real`, then the values are compared as defined by the standard.
    !>                          <li>    If the input arguments are of type `complex`, then the input values are elementally compared [lexically](@ref pm_arrayCompareLex) (similar to string comparison).
    !>                          <li>    If the input arguments are of type `logical`, then `.false. < .true.` is assumed.
    !>                      </ol>
    !>
    !>  \return
    !>  `anyMember`     :   The output scalar or vector of the same size as the input `val`, of type `logical` of default kind \LK.<br>
    !>                      It is `.true.` <b>if and only if any elements</b> of the input `val` are within the input range specified by `Set`.<br>
    !>
    !>  \interface{anyinrange}
    !>  \code{.F90}
    !>
    !>      use pm_arrayMembership, only: operator(.anyinrange.)
    !>
    !>      anyMember = val .anyinrange. Set(1:2) ! scalar character.
    !>      anyMember = val(:) .anyinrange. Set(1:2)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \note
    !>  Normally, the condition `Set(1) < Set(2)` is expected to hold for the input `Set` argument, otherwise, the specified range corresponds to an empty set.
    !>
    !>  \note
    !>  This generic interface offers a neat way of checking a value against a range returned by either
    !>  [getMinMax()](@ref pm_mathMinMax::getMinMax) or [operator(.subadd.)](@ref pm_mathSubAdd_subadd) or [operator(.divmul.)](@ref pm_mathDivMul_divmul).<br>
    !>
    !>  \see
    !>  [operator(.in.)](@ref pm_arrayMembership_in)<br>
    !>  [operator(.allin.)](@ref pm_arrayMembership_allin)<br>
    !>  [operator(.anyin.)](@ref pm_arrayMembership_anyin)<br>
    !>  [operator(.inrange.)](@ref pm_arrayMembership_inrange)<br>
    !>  [operator(.anyinrange.)](@ref pm_arrayMembership_anyinrange)<br>
    !>  [operator(.allinrange.)](@ref pm_arrayMembership_allinrange)<br>
    !>  [operator(.divmul.)](@ref pm_mathDivMul_divmul)<br>
    !>  [operator(.subadd.)](@ref pm_mathSubAdd_subadd)<br>
    !>  [getMinMax](@ref pm_mathMinMax::getMinMax)<br>
    !>  [setMinMax](@ref pm_mathMinMax::setMinMax)<br>
    !>  [pm_arraySearch](@ref pm_arraySearch)<br>
    !>  [pm_arrayFind](@ref pm_arrayFind)<br>
    !>
    !>  \example{anyinrange}
    !>  \include{lineno} example/pm_arrayMembership/anyinrange/main.F90
    !>  \compilef{anyinrange}
    !>  \output{anyinrange}
    !>  \include{lineno} example/pm_arrayMembership/anyinrange/main.out.F90
    !>
    !>  \test
    !>  [test_pm_arrayMembership](@ref test_pm_arrayMembership)
    !>
    !>  \finmain{anyinrange}
    !>
    !>  \author
    !>  \FatemehBagheri, Wednesday 1:35 PM, August 11, 2021, Dallas, TX
    interface operator(.anyinrange.)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function anyinrange_D0_D0_SK5(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D0_D0_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in)                :: val
        character(2,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: anyMember
    end function
#endif

#if SK4_ENABLED
    pure module function anyinrange_D0_D0_SK4(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D0_D0_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in)                :: val
        character(2,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: anyMember
    end function
#endif

#if SK3_ENABLED
    pure module function anyinrange_D0_D0_SK3(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D0_D0_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in)                :: val
        character(2,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: anyMember
    end function
#endif

#if SK2_ENABLED
    pure module function anyinrange_D0_D0_SK2(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D0_D0_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in)                :: val
        character(2,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: anyMember
    end function
#endif

#if SK1_ENABLED
    pure module function anyinrange_D0_D0_SK1(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D0_D0_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in)                :: val
        character(2,SKC)        , intent(in)                :: Set
        logical(LK)                                         :: anyMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    pure module function anyinrange_D1_D1_SK5(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_SK5
#endif
        use pm_kind, only: SKC => SK5
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if SK4_ENABLED
    pure module function anyinrange_D1_D1_SK4(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_SK4
#endif
        use pm_kind, only: SKC => SK4
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if SK3_ENABLED
    pure module function anyinrange_D1_D1_SK3(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_SK3
#endif
        use pm_kind, only: SKC => SK3
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if SK2_ENABLED
    pure module function anyinrange_D1_D1_SK2(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_SK2
#endif
        use pm_kind, only: SKC => SK2
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if SK1_ENABLED
    pure module function anyinrange_D1_D1_SK1(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_SK1
#endif
        use pm_kind, only: SKC => SK1
        character(*,SKC)        , intent(in), contiguous    :: val(:)
        character(*,SKC)        , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function anyinrange_D1_D1_IK5(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if IK4_ENABLED
    pure module function anyinrange_D1_D1_IK4(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if IK3_ENABLED
    pure module function anyinrange_D1_D1_IK3(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if IK2_ENABLED
    pure module function anyinrange_D1_D1_IK2(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if IK1_ENABLED
    pure module function anyinrange_D1_D1_IK1(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in), contiguous    :: val(:)
        integer(IKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    pure module function anyinrange_D1_D1_LK5(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_LK5
#endif
        use pm_kind, only: LKC => LK5
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if LK4_ENABLED
    pure module function anyinrange_D1_D1_LK4(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_LK4
#endif
        use pm_kind, only: LKC => LK4
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if LK3_ENABLED
    pure module function anyinrange_D1_D1_LK3(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_LK3
#endif
        use pm_kind, only: LKC => LK3
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if LK2_ENABLED
    pure module function anyinrange_D1_D1_LK2(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_LK2
#endif
        use pm_kind, only: LKC => LK2
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if LK1_ENABLED
    pure module function anyinrange_D1_D1_LK1(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_LK1
#endif
        use pm_kind, only: LKC => LK1
        logical(LKC)            , intent(in), contiguous    :: val(:)
        logical(LKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function anyinrange_D1_D1_CK5(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if CK4_ENABLED
    pure module function anyinrange_D1_D1_CK4(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if CK3_ENABLED
    pure module function anyinrange_D1_D1_CK3(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if CK2_ENABLED
    pure module function anyinrange_D1_D1_CK2(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if CK1_ENABLED
    pure module function anyinrange_D1_D1_CK1(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in), contiguous    :: val(:)
        complex(CKC)            , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function anyinrange_D1_D1_RK5(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if RK4_ENABLED
    pure module function anyinrange_D1_D1_RK4(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if RK3_ENABLED
    pure module function anyinrange_D1_D1_RK3(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if RK2_ENABLED
    pure module function anyinrange_D1_D1_RK2(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

#if RK1_ENABLED
    pure module function anyinrange_D1_D1_RK1(val, Set) result(anyMember)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: anyinrange_D1_D1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in), contiguous    :: val(:)
        real(RKC)               , intent(in)                :: Set(2)
        logical(LK)                                         :: anyMember
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_arrayMembership