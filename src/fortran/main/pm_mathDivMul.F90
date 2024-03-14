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
!>  This module contains procedures and generic interfaces for evaluating the mathematical
!>  division and multiplication operators acting on integer, complex, or real values.<br>
!>
!>  \note
!>  The procedures of this module offer a handy and flexible way of membership checks.
!>
!>  \see
!>  [operator(.divmul.)](@ref pm_mathDivMul_divmul)<br>
!>  [operator(.subadd.)](@ref pm_mathSubAdd_subadd)<br>
!>  [operator(.allin.)](@ref pm_arrayMembership_allin)<br>
!>  [operator(.anyin.)](@ref pm_arrayMembership_anyin)<br>
!>  [operator(.inrange.)](@ref pm_arrayMembership_inrange)<br>
!>  [operator(.anyinrange.)](@ref pm_arrayMembership_anyinrange)<br>
!>  [operator(.allinrange.)](@ref pm_arrayMembership_allinrange)<br>
!>  [operator(.in.)](@ref pm_arrayMembership_in)<br>
!>  [getMinMax](@ref pm_mathMinMax::getMinMax)<br>
!>  [setMinMax](@ref pm_mathMinMax::setMinMax)<br>
!>
!>  \test
!>  [test_pm_mathDivMul](@ref test_pm_mathDivMul)
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathDivMul

    use pm_kind, only: IK, RK, SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathDivMul"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_mathDivMul_divmul
    !>  Generate and return the result of applying the mathematical unary or binary operators \f$\times/\f$ to the input argument(s).<br>
    !>
    !>  \details
    !>  Given two input arguments `ref` and `val`, the procedures of this generic interface return an array of size `2` whose elements are `[ref / val, ref * val]`.<br>
    !>  If `ref` is missing, then an appropriate default value is used.<br>
    !>
    !>  \param[in]  ref :   The input scalar of the same type and kind as the input `val`, representing the reference value in the operations.<br>
    !>                      (**optional**, default = `1.`)
    !>  \param[in]  val :   The input scalar of
    !>                      <ul>
    !>                          <li>    type `integer` of kind \IKALL, or
    !>                          <li>    type `complex` of kind \CKALL, or
    !>                          <li>    type `real` of kind \RKALL,
    !>                      </ul>
    !>                      representing the value upon which the operators \f$\times/\f$ will act.
    !>
    !>  \return
    !>  `binval`        :   The output vector of size `2` of the same type and kind as `val`, containing the result of the \f$\times/\f$ operations.
    !>
    !>  \interface{divmul}
    !>  \code{.F90}
    !>
    !>      use pm_mathDivMul, only: operator(.divmul.)
    !>
    !>      binval(1:2) = .divmul. val
    !>      binval(1:2) = ref .divmul. val
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `val /= 0` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [operator(.divmul.)](@ref pm_mathDivMul_divmul)<br>
    !>  [operator(.subadd.)](@ref pm_mathSubAdd_subadd)<br>
    !>  [operator(.allin.)](@ref pm_arrayMembership_allin)<br>
    !>  [operator(.anyin.)](@ref pm_arrayMembership_anyin)<br>
    !>  [operator(.inrange.)](@ref pm_arrayMembership_inrange)<br>
    !>  [operator(.anyinrange.)](@ref pm_arrayMembership_anyinrange)<br>
    !>  [operator(.allinrange.)](@ref pm_arrayMembership_allinrange)<br>
    !>  [operator(.in.)](@ref pm_arrayMembership_in)<br>
    !>  [getMinMax](@ref pm_mathMinMax::getMinMax)<br>
    !>  [setMinMax](@ref pm_mathMinMax::setMinMax)<br>
    !>
    !>  \example{divmul}
    !>  \include{lineno} example/pm_mathDivMul/divmul/main.F90
    !>  \compilef{divmul}
    !>  \output{divmul}
    !>  \include{lineno} example/pm_mathDivMul/divmul/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathDivMul](@ref test_pm_mathDivMul)
    !>
    !>  \finmain{divmul}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface operator(.divmul.)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getDivMulUnary_IK5(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulUnary_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC), intent(in)    :: val
        integer(IKC)                :: binval(2)
    end function
#endif

#if IK4_ENABLED
    PURE module function getDivMulUnary_IK4(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulUnary_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC), intent(in)    :: val
        integer(IKC)                :: binval(2)
    end function
#endif

#if IK3_ENABLED
    PURE module function getDivMulUnary_IK3(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulUnary_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC), intent(in)    :: val
        integer(IKC)                :: binval(2)
    end function
#endif

#if IK2_ENABLED
    PURE module function getDivMulUnary_IK2(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulUnary_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC), intent(in)    :: val
        integer(IKC)                :: binval(2)
    end function
#endif

#if IK1_ENABLED
    PURE module function getDivMulUnary_IK1(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulUnary_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC), intent(in)    :: val
        integer(IKC)                :: binval(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getDivMulUnary_CK5(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulUnary_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC), intent(in)    :: val
        complex(CKC)                :: binval(2)
    end function
#endif

#if CK4_ENABLED
    PURE module function getDivMulUnary_CK4(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulUnary_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC), intent(in)    :: val
        complex(CKC)                :: binval(2)
    end function
#endif

#if CK3_ENABLED
    PURE module function getDivMulUnary_CK3(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulUnary_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC), intent(in)    :: val
        complex(CKC)                :: binval(2)
    end function
#endif

#if CK2_ENABLED
    PURE module function getDivMulUnary_CK2(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulUnary_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC), intent(in)    :: val
        complex(CKC)                :: binval(2)
    end function
#endif

#if CK1_ENABLED
    PURE module function getDivMulUnary_CK1(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulUnary_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC), intent(in)    :: val
        complex(CKC)                :: binval(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDivMulUnary_RK5(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulUnary_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)    :: val
        real(RKC)                   :: binval(2)
    end function
#endif

#if RK4_ENABLED
    PURE module function getDivMulUnary_RK4(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulUnary_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)    :: val
        real(RKC)                   :: binval(2)
    end function
#endif

#if RK3_ENABLED
    PURE module function getDivMulUnary_RK3(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulUnary_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)    :: val
        real(RKC)                   :: binval(2)
    end function
#endif

#if RK2_ENABLED
    PURE module function getDivMulUnary_RK2(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulUnary_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)    :: val
        real(RKC)                   :: binval(2)
    end function
#endif

#if RK1_ENABLED
    PURE module function getDivMulUnary_RK1(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulUnary_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)    :: val
        real(RKC)                   :: binval(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getDivMulBinary_IK5(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulBinary_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC), intent(in)    :: ref, val
        integer(IKC)                :: binval(2)
    end function
#endif

#if IK4_ENABLED
    PURE module function getDivMulBinary_IK4(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulBinary_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC), intent(in)    :: ref, val
        integer(IKC)                :: binval(2)
    end function
#endif

#if IK3_ENABLED
    PURE module function getDivMulBinary_IK3(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulBinary_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC), intent(in)    :: ref, val
        integer(IKC)                :: binval(2)
    end function
#endif

#if IK2_ENABLED
    PURE module function getDivMulBinary_IK2(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulBinary_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC), intent(in)    :: ref, val
        integer(IKC)                :: binval(2)
    end function
#endif

#if IK1_ENABLED
    PURE module function getDivMulBinary_IK1(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulBinary_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC), intent(in)    :: ref, val
        integer(IKC)                :: binval(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getDivMulBinary_CK5(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulBinary_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC), intent(in)    :: ref, val
        complex(CKC)                :: binval(2)
    end function
#endif

#if CK4_ENABLED
    PURE module function getDivMulBinary_CK4(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulBinary_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC), intent(in)    :: ref, val
        complex(CKC)                :: binval(2)
    end function
#endif

#if CK3_ENABLED
    PURE module function getDivMulBinary_CK3(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulBinary_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC), intent(in)    :: ref, val
        complex(CKC)                :: binval(2)
    end function
#endif

#if CK2_ENABLED
    PURE module function getDivMulBinary_CK2(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulBinary_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC), intent(in)    :: ref, val
        complex(CKC)                :: binval(2)
    end function
#endif

#if CK1_ENABLED
    PURE module function getDivMulBinary_CK1(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulBinary_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC), intent(in)    :: ref, val
        complex(CKC)                :: binval(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getDivMulBinary_RK5(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulBinary_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)    :: ref, val
        real(RKC)                   :: binval(2)
    end function
#endif

#if RK4_ENABLED
    PURE module function getDivMulBinary_RK4(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulBinary_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)    :: ref, val
        real(RKC)                   :: binval(2)
    end function
#endif

#if RK3_ENABLED
    PURE module function getDivMulBinary_RK3(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulBinary_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)    :: ref, val
        real(RKC)                   :: binval(2)
    end function
#endif

#if RK2_ENABLED
    PURE module function getDivMulBinary_RK2(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulBinary_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)    :: ref, val
        real(RKC)                   :: binval(2)
    end function
#endif

#if RK1_ENABLED
    PURE module function getDivMulBinary_RK1(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDivMulBinary_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)    :: ref, val
        real(RKC)                   :: binval(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathDivMul ! LCOV_EXCL_LINE