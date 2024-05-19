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
!>  This module contains procedures and generic interfaces for evaluating the mathematical operator \f$\mp\f$ acting on integer, complex, or real values.<br>
!>
!>  \details
!>  The plus–minus sign, `±`, is a mathematical symbol with multiple meanings.<br>
!>  <ol>
!>      <li>    In mathematics, it generally indicates a choice of exactly two possible values,
!>              one of which is obtained through addition and the other through subtraction.
!>      <li>    In experimental sciences, the sign commonly indicates the confidence interval or uncertainty bounding
!>              a range of possible errors in a measurement, often the standard deviation or standard error.<br>
!>              The sign may also represent an inclusive range of values that a reading might have.
!>      <li>    In medicine, it means **with or without**.
!>      <li>    In engineering, the sign indicates the tolerance, which is the range of values that are
!>              considered to be acceptable, safe, or which comply with some standard or with a contract.
!>      <li>    In botany, it is used in morphological descriptions to notate **more or less**.
!>      <li>    In chemistry, the sign is used to indicate a racemic mixture.
!>      <li>    In electronics, this sign may indicate a dual voltage power supply,
!>              such as `±5` volts means `+5` volts and `-5` volts, when used with audio circuits and operational amplifiers.
!>  </ol>
!>  Given two input arguments `ref` and `val`, the procedures of this module return an array of size `2` whose elements are `[ref - val, ref + val]`.<br>
!>  If `ref` is missing, then an appropriate default value is used.<br>
!>
!>  \note
!>  The procedures of this module offer a handy and flexible way of membership checks via [operator(.inrange.)](@ref pm_arrayMembership_inrange).
!>
!>  \note
!>  The operator \f$\pm\f$ performs the opposite operation of \f$\mp\f$, that is, \f$\pm = -\mp\f$.<br>
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
!>  [test_pm_mathSubAdd](@ref test_pm_mathSubAdd)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathSubAdd

    use pm_kind, only: IK, RK, SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathSubAdd"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \anchor pm_mathSubAdd_subadd
    !>  Generate and return the result of applying the mathematical unary or binary operator \f$\mp\f$ to the input argument(s).
    !>
    !>  \details
    !>  Given two input arguments `ref` and `val`, the procedures of this generic interface return an array of size `2` whose elements are `[ref - val, ref + val]`.<br>
    !>  If `ref` is missing, then an appropriate default value is used.<br>
    !>
    !>  \param[in]  ref :   The input scalar of the same type and kind as the input `val`, representing the reference value in the operation.<br>
    !>                      (**optional**, default = `0.`)
    !>  \param[in]  val :   The input scalar of
    !>                      <ul>
    !>                          <li>    type `integer` of kind \IKALL, or
    !>                          <li>    type `complex` of kind \CKALL, or
    !>                          <li>    type `real` of kind \RKALL,
    !>                      </ul>
    !>                      representing the value upon which the operator \f$\mp\f$ will act.
    !>
    !>  \return
    !>  `binval`        :   The output vector of size `2` of the same type and kind as `val`, containing the result of the \f$\mp\f$ operation.
    !>
    !>  \interface{subadd}
    !>  \code{.F90}
    !>
    !>      use pm_mathSubAdd, only: operator(.subadd.)
    !>
    !>      binval(1:2) = .subadd. val
    !>      binval(1:2) = ref .subadd. val
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \note
    !>  The procedures of this module offer a handy and flexible way of membership checks via [operator(.inrange.)](@ref pm_arrayMembership_inrange).<br>
    !>
    !>  \note
    !>  The operator \f$\mp\f$ performs the opposite operation of \f$\mp\f$, that is, \f$\pm = -\mp\f$.<br>
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
    !>  \example{subadd}
    !>  \include{lineno} example/pm_mathSubAdd/subadd/main.F90
    !>  \compilef{subadd}
    !>  \output{subadd}
    !>  \include{lineno} example/pm_mathSubAdd/subadd/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathSubAdd](@ref test_pm_mathSubAdd)
    !>
    !>  \final{subadd}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface operator(.subadd.)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function getSubAddUnary_IK5(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddUnary_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG), intent(in)    :: val
        integer(IKG)                :: binval(2)
    end function
#endif

#if IK4_ENABLED
    pure module function getSubAddUnary_IK4(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddUnary_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG), intent(in)    :: val
        integer(IKG)                :: binval(2)
    end function
#endif

#if IK3_ENABLED
    pure module function getSubAddUnary_IK3(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddUnary_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG), intent(in)    :: val
        integer(IKG)                :: binval(2)
    end function
#endif

#if IK2_ENABLED
    pure module function getSubAddUnary_IK2(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddUnary_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG), intent(in)    :: val
        integer(IKG)                :: binval(2)
    end function
#endif

#if IK1_ENABLED
    pure module function getSubAddUnary_IK1(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddUnary_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG), intent(in)    :: val
        integer(IKG)                :: binval(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function getSubAddUnary_CK5(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddUnary_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    :: val
        complex(CKG)                :: binval(2)
    end function
#endif

#if CK4_ENABLED
    pure module function getSubAddUnary_CK4(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddUnary_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    :: val
        complex(CKG)                :: binval(2)
    end function
#endif

#if CK3_ENABLED
    pure module function getSubAddUnary_CK3(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddUnary_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    :: val
        complex(CKG)                :: binval(2)
    end function
#endif

#if CK2_ENABLED
    pure module function getSubAddUnary_CK2(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddUnary_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    :: val
        complex(CKG)                :: binval(2)
    end function
#endif

#if CK1_ENABLED
    pure module function getSubAddUnary_CK1(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddUnary_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    :: val
        complex(CKG)                :: binval(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function getSubAddUnary_RK5(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddUnary_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: val
        real(RKG)                   :: binval(2)
    end function
#endif

#if RK4_ENABLED
    pure module function getSubAddUnary_RK4(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddUnary_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: val
        real(RKG)                   :: binval(2)
    end function
#endif

#if RK3_ENABLED
    pure module function getSubAddUnary_RK3(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddUnary_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: val
        real(RKG)                   :: binval(2)
    end function
#endif

#if RK2_ENABLED
    pure module function getSubAddUnary_RK2(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddUnary_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: val
        real(RKG)                   :: binval(2)
    end function
#endif

#if RK1_ENABLED
    pure module function getSubAddUnary_RK1(val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddUnary_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: val
        real(RKG)                   :: binval(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure module function getSubAddBinary_IK5(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddBinary_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG), intent(in)    :: ref, val
        integer(IKG)                :: binval(2)
    end function
#endif

#if IK4_ENABLED
    pure module function getSubAddBinary_IK4(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddBinary_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG), intent(in)    :: ref, val
        integer(IKG)                :: binval(2)
    end function
#endif

#if IK3_ENABLED
    pure module function getSubAddBinary_IK3(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddBinary_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG), intent(in)    :: ref, val
        integer(IKG)                :: binval(2)
    end function
#endif

#if IK2_ENABLED
    pure module function getSubAddBinary_IK2(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddBinary_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG), intent(in)    :: ref, val
        integer(IKG)                :: binval(2)
    end function
#endif

#if IK1_ENABLED
    pure module function getSubAddBinary_IK1(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddBinary_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG), intent(in)    :: ref, val
        integer(IKG)                :: binval(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure module function getSubAddBinary_CK5(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddBinary_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    :: ref, val
        complex(CKG)                :: binval(2)
    end function
#endif

#if CK4_ENABLED
    pure module function getSubAddBinary_CK4(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddBinary_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    :: ref, val
        complex(CKG)                :: binval(2)
    end function
#endif

#if CK3_ENABLED
    pure module function getSubAddBinary_CK3(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddBinary_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    :: ref, val
        complex(CKG)                :: binval(2)
    end function
#endif

#if CK2_ENABLED
    pure module function getSubAddBinary_CK2(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddBinary_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    :: ref, val
        complex(CKG)                :: binval(2)
    end function
#endif

#if CK1_ENABLED
    pure module function getSubAddBinary_CK1(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddBinary_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    :: ref, val
        complex(CKG)                :: binval(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure module function getSubAddBinary_RK5(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddBinary_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: ref, val
        real(RKG)                   :: binval(2)
    end function
#endif

#if RK4_ENABLED
    pure module function getSubAddBinary_RK4(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddBinary_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: ref, val
        real(RKG)                   :: binval(2)
    end function
#endif

#if RK3_ENABLED
    pure module function getSubAddBinary_RK3(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddBinary_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: ref, val
        real(RKG)                   :: binval(2)
    end function
#endif

#if RK2_ENABLED
    pure module function getSubAddBinary_RK2(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddBinary_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: ref, val
        real(RKG)                   :: binval(2)
    end function
#endif

#if RK1_ENABLED
    pure module function getSubAddBinary_RK1(ref, val) result(binval)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSubAddBinary_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: ref, val
        real(RKG)                   :: binval(2)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathSubAdd ! LCOV_EXCL_LINE