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
!>  This module contains procedures and generic interfaces for testing for exceptional cases at runtime.<br>
!>
!>  \details
!>  Three exceptional cases are currently handled by the procedures of this module,
!>  <ol>
!>      <li>    **Addition Outflow**: Testing whether numerical negative or positive overflow can occur in an addition operation for arbitrary intrinsic numeric types and kinds,
!>              without the risk of facing runtime positive or negative outflow (overflow) errors.<br>
!>              There are two scenarios leading to an addition `a + b` outflow:
!>              <ol>
!>                  <li>    <b>Negative Overflow</b> occurs if and only if `a < 0 && b < 0 && a < min - b`.
!>                  <li>    <b>Positive Overflow</b> occurs if and only if `0 < a && 0 < b && max - b < a`.
!>              </ol>
!>              where `min` and `max` refer to the minimum and maximum possible values for the type and kind of the pair `(a,b)`.<br>
!>      <li>    **Infinity**: obtaining or testing for occurrences of `real` and `complex` IEEE-compliant infinity values.<br>
!>      <li>    **Not-A-Number (NAN)**: obtaining or testing for the occurrences of `real` and `complex` IEEE-compliant `NaN` (Not A Number) values.
!>  </ol>
!>
!>
!>  \test
!>  [test_pm_except](@ref test_pm_except)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_except

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_except"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if and only if the addition of the
    !>  two input `integer`, `complex` , or `real` values causes runtime negative or positive outflow (overflow).
    !>
    !>  \details
    !>  See the documentation of [pm_except](@ref pm_except) for more information.
    !>
    !>  \param[in]  a   :   The input scalar, or array of the same rank and shape as other array-like arguments, of
    !>                      <ol>
    !>                          <li>    type `integer` of kind \IKALL,
    !>                          <li>    type `complex` of kind \CKALL,
    !>                          <li>    type `real` of kind \RKALL,
    !>                      </ol>
    !>                      containing the number to be added to the other input argument `b`.
    !>  \param[in]  b   :   The input scalar, or array of the same rank and shape as other array-like arguments,
    !>                      of the same type and kind as the input argument `a`
    !>                      containing the number to be added to the other input argument `a`.
    !>
    !>  \return
    !>  `outflow`       :   The scalar, or array of the same rank and shape as other array-like arguments,
    !>                      of type `logical` of default kind \LK that is `.true.` if and only if `-huge(a) < a + b .or. huge(a) < a + b`.
    !>
    !>  \interface{isAddOutflow}
    !>  \code{.F90}
    !>
    !>      use pm_except, only: isAddOutflow
    !>      use pm_kind, only: LK
    !>      logical(LK) :: outflow
    !>
    !>      outflow = isAddOutflow(a, b)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  A complex addition is considered an outflow if either real or imaginary component addition causes an outflow.
    !>
    !>  \see
    !>  [getFactorial](@ref pm_mathFactorial::getFactorial)<br>
    !>
    !>  \example{isAddOutflow}
    !>  \include{lineno} example/pm_except/isAddOutflow/main.F90
    !>  \compilef{isAddOutflow}
    !>  \output{isAddOutflow}
    !>  \include{lineno} example/pm_except/isAddOutflow/main.out.F90
    !>
    !>  \test
    !>  [test_pm_except](@ref test_pm_except)
    !>
    !>  \final{isAddOutflow}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface isAddOutflow

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure elemental module function isAddOutflow_IK5(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflow_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if IK4_ENABLED
    pure elemental module function isAddOutflow_IK4(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflow_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if IK3_ENABLED
    pure elemental module function isAddOutflow_IK3(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflow_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if IK2_ENABLED
    pure elemental module function isAddOutflow_IK2(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflow_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if IK1_ENABLED
    pure elemental module function isAddOutflow_IK1(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflow_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isAddOutflow_CK5(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflow_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isAddOutflow_CK4(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflow_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isAddOutflow_CK3(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflow_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isAddOutflow_CK2(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflow_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isAddOutflow_CK1(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflow_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function isAddOutflow_RK5(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflow_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)       , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if RK4_ENABLED
    pure elemental module function isAddOutflow_RK4(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflow_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)       , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if RK3_ENABLED
    pure elemental module function isAddOutflow_RK3(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflow_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)       , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if RK2_ENABLED
    pure elemental module function isAddOutflow_RK2(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflow_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)       , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if RK1_ENABLED
    pure elemental module function isAddOutflow_RK1(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflow_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)       , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if and only if the addition of the
    !>  two input `integer`, `complex` , or `real` values causes runtime negative outflow (overflow).
    !>
    !>  \details
    !>  See the documentation of [pm_except](@ref pm_except) for more information.
    !>
    !>  \param[in]  a   :   The input scalar, or array of the same rank and shape as other array-like arguments, of
    !>                      <ol>
    !>                          <li>    type `integer` of kind \IKALL,
    !>                          <li>    type `complex` of kind \CKALL,
    !>                          <li>    type `real` of kind \RKALL,
    !>                      </ol>
    !>                      containing the number to be added to the other input argument `b`.
    !>  \param[in]  b   :   The input scalar, or array of the same rank and shape as other array-like arguments,
    !>                      of the same type and kind as the input argument `a`
    !>                      containing the number to be added to the other input argument `a`.
    !>
    !>  \return
    !>  `outflow`       :   The scalar, or array of the same rank and shape as other array-like arguments,
    !>                      of type `logical` of default kind \LK that is `.true.` if and only if `-huge(a) < a + b`.
    !>
    !>  \interface{isAddOutflowNeg}
    !>  \code{.F90}
    !>
    !>      use pm_except, only: isAddOutflowNeg
    !>      use pm_kind, only: LK
    !>      logical(LK) :: outflow
    !>
    !>      outflow = isAddOutflowNeg(a, b)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  A complex addition is considered an outflow if either real or imaginary component addition causes an outflow.
    !>
    !>  \see
    !>  [getFactorial](@ref pm_mathFactorial::getFactorial)<br>
    !>
    !>  \example{isAddOutflowNeg}
    !>  \include{lineno} example/pm_except/isAddOutflowNeg/main.F90
    !>  \compilef{isAddOutflowNeg}
    !>  \output{isAddOutflowNeg}
    !>  \include{lineno} example/pm_except/isAddOutflowNeg/main.out.F90
    !>
    !>  \test
    !>  [test_pm_except](@ref test_pm_except)
    !>
    !>  \final{isAddOutflowNeg}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface isAddOutflowNeg

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure elemental module function isAddOutflowNeg_IK5(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowNeg_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if IK4_ENABLED
    pure elemental module function isAddOutflowNeg_IK4(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowNeg_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if IK3_ENABLED
    pure elemental module function isAddOutflowNeg_IK3(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowNeg_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if IK2_ENABLED
    pure elemental module function isAddOutflowNeg_IK2(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowNeg_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if IK1_ENABLED
    pure elemental module function isAddOutflowNeg_IK1(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowNeg_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isAddOutflowNeg_CK5(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowNeg_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isAddOutflowNeg_CK4(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowNeg_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isAddOutflowNeg_CK3(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowNeg_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isAddOutflowNeg_CK2(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowNeg_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isAddOutflowNeg_CK1(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowNeg_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function isAddOutflowNeg_RK5(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowNeg_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)       , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if RK4_ENABLED
    pure elemental module function isAddOutflowNeg_RK4(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowNeg_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)       , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if RK3_ENABLED
    pure elemental module function isAddOutflowNeg_RK3(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowNeg_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)       , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if RK2_ENABLED
    pure elemental module function isAddOutflowNeg_RK2(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowNeg_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)       , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if RK1_ENABLED
    pure elemental module function isAddOutflowNeg_RK1(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowNeg_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)       , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if and only if the addition of the
    !>  two input `integer`, `complex` , or `real` values causes runtime positive outflow (overflow).
    !>
    !>  \details
    !>  See the documentation of [pm_except](@ref pm_except) for more information.
    !>
    !>  \param[in]  a   :   The input scalar, or array of the same rank and shape as other array-like arguments, of
    !>                      <ol>
    !>                          <li>    type `integer` of kind \IKALL,
    !>                          <li>    type `complex` of kind \CKALL,
    !>                          <li>    type `real` of kind \RKALL,
    !>                      </ol>
    !>                      containing the number to be added to the other input argument `b`.
    !>  \param[in]  b   :   The input scalar, or array of the same rank and shape as other array-like arguments,
    !>                      of the same type and kind as the input argument `a`
    !>                      containing the number to be added to the other input argument `a`.
    !>
    !>  \return
    !>  `outflow`       :   The scalar, or array of the same rank and shape as other array-like arguments,
    !>                      of type `logical` of default kind \LK that is `.true.` if and only if `huge(a) < a + b`.
    !>
    !>  \interface{isAddOutflowPos}
    !>  \code{.F90}
    !>
    !>      use pm_except, only: isAddOutflowPos
    !>      use pm_kind, only: LK
    !>      logical(LK) :: outflow
    !>
    !>      outflow = isAddOutflowPos(a, b)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  A complex addition is considered an outflow if either real or imaginary component addition causes an outflow.
    !>
    !>  \see
    !>  [getFactorial](@ref pm_mathFactorial::getFactorial)<br>
    !>
    !>  \example{isAddOutflowPos}
    !>  \include{lineno} example/pm_except/isAddOutflowPos/main.F90
    !>  \compilef{isAddOutflowPos}
    !>  \output{isAddOutflowPos}
    !>  \include{lineno} example/pm_except/isAddOutflowPos/main.out.F90
    !>
    !>  \test
    !>  [test_pm_except](@ref test_pm_except)
    !>
    !>  \final{isAddOutflowPos}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface isAddOutflowPos

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    pure elemental module function isAddOutflowPos_IK5(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowPos_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if IK4_ENABLED
    pure elemental module function isAddOutflowPos_IK4(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowPos_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if IK3_ENABLED
    pure elemental module function isAddOutflowPos_IK3(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowPos_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if IK2_ENABLED
    pure elemental module function isAddOutflowPos_IK2(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowPos_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if IK1_ENABLED
    pure elemental module function isAddOutflowPos_IK1(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowPos_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isAddOutflowPos_CK5(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowPos_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isAddOutflowPos_CK4(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowPos_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isAddOutflowPos_CK3(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowPos_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isAddOutflowPos_CK2(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowPos_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isAddOutflowPos_CK1(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowPos_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)    , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function isAddOutflowPos_RK5(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowPos_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)       , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if RK4_ENABLED
    pure elemental module function isAddOutflowPos_RK4(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowPos_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)       , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if RK3_ENABLED
    pure elemental module function isAddOutflowPos_RK3(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowPos_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)       , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if RK2_ENABLED
    pure elemental module function isAddOutflowPos_RK2(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowPos_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)       , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

#if RK1_ENABLED
    pure elemental module function isAddOutflowPos_RK1(a, b) result(outflow)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isAddOutflowPos_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)       , intent(in)        :: a, b
        logical(LK)                         :: outflow
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input value is an **IEEE-compliant** negative or positive infinity.
    !>  If the input value is a `complex` number, then the output is `.true.` if any of the two real or
    !>  imaginary components or both components are either negative or positive infinities.
    !>
    !>  \param[in]  x   :   The input scalar or array of arbitrary rank of either <br>
    !>                      &nbsp; type `complex` of kind \CKALL, or <br>
    !>                      &nbsp; type `real` of kind \RKALL, <br>
    !>                      &nbsp; whose value will be tested for being negative or positive infinity.
    !>
    !>  \return
    !>  `inf`           :   The output scalar or array of the same shape as the input `x` of type `logical` of default kind \LK
    !>                      whose value is `.true.` if the input `x` is a negative or positive infinity, otherwise it is `.false.`.
    !>
    !>  \interface{isInf}
    !>  \code{.F90}
    !>
    !>      use pm_except, only: isInf
    !>      use pm_kind, only: LK
    !>      logical(LK) :: inf
    !>
    !>      inf = isInf(x)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  A very simple test of negative or positive infinity of a variable `x` is the condition `abs(x) > huge(x)` that is `.true.` if `x` is a negative or positive infinity.
    !>
    !>  \note
    !>  The procedures under this generic interface are equivalent to `.not. (ieee_is_negative(x) .or. ieee_is_finite(x))`
    !>  from the `ieee_arithmetic` intrinsic module to detect negative or positive infinity.
    !>  This generic interface extends this function also to `complex` numbers.
    !>
    !>  \see
    !>  [isInf](@ref pm_except::isInf)<br>
    !>  [isNAN](@ref pm_except::isNAN)<br>
    !>  [getNAN](@ref pm_except::getNAN)<br>
    !>  [setNAN](@ref pm_except::setNAN)<br>
    !>  [isInfPos](@ref pm_except::isInfPos)<br>
    !>  [isInfNeg](@ref pm_except::isInfNeg)<br>
    !>  [getInfPos](@ref pm_except::getInfPos)<br>
    !>  [setInfPos](@ref pm_except::setInfPos)<br>
    !>  [getInfNeg](@ref pm_except::getInfNeg)<br>
    !>  [setInfNeg](@ref pm_except::setInfNeg)<br>
    !>
    !>  \example{isInf}
    !>  \include{lineno} example/pm_except/isInf/main.F90
    !>  \compilef{isInf}
    !>  \output{isInf}
    !>  \include{lineno} example/pm_except/isInf/main.out.F90
    !>
    !>  \test
    !>  [test_pm_except](@ref test_pm_except)
    !>
    !>  \final{isInf}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface isInf

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function isInf_RK5(x) result(inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInf_RK5
#endif
        use pm_kind, only: LK, RKG => RK5
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: inf
    end function
#endif

#if RK4_ENABLED
    pure elemental module function isInf_RK4(x) result(inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInf_RK4
#endif
        use pm_kind, only: LK, RKG => RK4
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: inf
    end function
#endif

#if RK3_ENABLED
    pure elemental module function isInf_RK3(x) result(inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInf_RK3
#endif
        use pm_kind, only: LK, RKG => RK3
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: inf
    end function
#endif

#if RK2_ENABLED
    pure elemental module function isInf_RK2(x) result(inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInf_RK2
#endif
        use pm_kind, only: LK, RKG => RK2
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: inf
    end function
#endif

#if RK1_ENABLED
    pure elemental module function isInf_RK1(x) result(inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInf_RK1
#endif
        use pm_kind, only: LK, RKG => RK1
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: inf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isInf_CK5(x) result(inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInf_CK5
#endif
        use pm_kind, only: LK, CKG => CK5
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: inf
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isInf_CK4(x) result(inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInf_CK4
#endif
        use pm_kind, only: LK, CKG => CK4
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: inf
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isInf_CK3(x) result(inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInf_CK3
#endif
        use pm_kind, only: LK, CKG => CK3
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: inf
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isInf_CK2(x) result(inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInf_CK2
#endif
        use pm_kind, only: LK, CKG => CK2
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: inf
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isInf_CK1(x) result(inf)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInf_CK1
#endif
        use pm_kind, only: LK, CKG => CK1
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: inf
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input value is an **IEEE-compliant** positive infinity.<br>
    !>
    !>  \details
    !>  If the input value is a `complex` number, then the output is `.true.` if any of the two real or
    !>  imaginary components or both are positive infinities.
    !>
    !>  \param[in]  x   :   The input scalar or array of arbitrary rank of either <br>
    !>                      <ol>
    !>                          <li>    type `complex` of kind \CKALL, or <br>
    !>                          <li>    type `real` of kind \RKALL, <br>
    !>                      </ol>
    !>                      whose value will be tested for being positive infinity.
    !>
    !>  \return
    !>  `infPos`        :   The output scalar or array of the same shape as the input `x` of type `logical` of default kind \LK
    !>                      whose value is `.true.` if the input `x` is a positive infinity, otherwise it is `.false.`.
    !>
    !>  \interface{isInfPos}
    !>  \code{.F90}
    !>
    !>      use pm_except, only: isInfPos
    !>      use pm_kind, only: LK
    !>      logical(LK) :: infPos
    !>
    !>      infPos = isInfPos(x)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  A very simple test of positive infinity of a variable `x` is the condition `x > huge(x)` that is `.true.` if `x` is a positive infinity.
    !>
    !>  \note
    !>  The procedures under this generic interface are equivalent to `.not. (ieee_is_negative(x) .or. ieee_is_finite(x))`
    !>  from the `ieee_arithmetic` intrinsic module to detect positive infinity.<br>
    !>  This generic interface extends this function also to `complex` numbers.<br>
    !>
    !>  \see
    !>  [isInf](@ref pm_except::isInf)<br>
    !>  [isNAN](@ref pm_except::isNAN)<br>
    !>  [getNAN](@ref pm_except::getNAN)<br>
    !>  [setNAN](@ref pm_except::setNAN)<br>
    !>  [isInfPos](@ref pm_except::isInfPos)<br>
    !>  [isInfNeg](@ref pm_except::isInfNeg)<br>
    !>  [getInfPos](@ref pm_except::getInfPos)<br>
    !>  [setInfPos](@ref pm_except::setInfPos)<br>
    !>  [getInfNeg](@ref pm_except::getInfNeg)<br>
    !>  [setInfNeg](@ref pm_except::setInfNeg)<br>
    !>
    !>  \example{isInfPos}
    !>  \include{lineno} example/pm_except/isInfPos/main.F90
    !>  \compilef{isInfPos}
    !>  \output{isInfPos}
    !>  \include{lineno} example/pm_except/isInfPos/main.out.F90
    !>
    !>  \test
    !>  [test_pm_except](@ref test_pm_except)
    !>
    !>  \final{isInfPos}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface isInfPos

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function isInfPos_RK5(x) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfPos_RK5
#endif
        use pm_kind, only: LK, RKG => RK5
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: infPos
    end function
#endif

#if RK4_ENABLED
    pure elemental module function isInfPos_RK4(x) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfPos_RK4
#endif
        use pm_kind, only: LK, RKG => RK4
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: infPos
    end function
#endif

#if RK3_ENABLED
    pure elemental module function isInfPos_RK3(x) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfPos_RK3
#endif
        use pm_kind, only: LK, RKG => RK3
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: infPos
    end function
#endif

#if RK2_ENABLED
    pure elemental module function isInfPos_RK2(x) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfPos_RK2
#endif
        use pm_kind, only: LK, RKG => RK2
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: infPos
    end function
#endif

#if RK1_ENABLED
    pure elemental module function isInfPos_RK1(x) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfPos_RK1
#endif
        use pm_kind, only: LK, RKG => RK1
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: infPos
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isInfPos_CK5(x) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfPos_CK5
#endif
        use pm_kind, only: LK, CKG => CK5
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: infPos
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isInfPos_CK4(x) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfPos_CK4
#endif
        use pm_kind, only: LK, CKG => CK4
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: infPos
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isInfPos_CK3(x) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfPos_CK3
#endif
        use pm_kind, only: LK, CKG => CK3
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: infPos
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isInfPos_CK2(x) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfPos_CK2
#endif
        use pm_kind, only: LK, CKG => CK2
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: infPos
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isInfPos_CK1(x) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfPos_CK1
#endif
        use pm_kind, only: LK, CKG => CK1
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: infPos
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an IEEE-compliant positive infinity.<br>
    !>
    !>  \details
    !>  This procedure is useful for deliberate initializations of `real` or `complex` numbers to positive infinity
    !>  without signaling the occurrence of an exception.
    !>
    !>  \param[out] mold    :   The input scalar or array of arbitrary rank of either<br>
    !>                          <ol>
    !>                              <li>    type `complex` of kind \CKALL, or <br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ol>
    !>                          whose type and kind will match the generated output `+Inf`.
    !>
    !>  \return
    !>  `infPos`            :   The output scalar or array of the same type and kind as
    !>                          the input `mold` containing the generated positive infinity value.
    !>
    !>  \interface{getInfPos}
    !>  \code{.F90}
    !>
    !>      use pm_except, only: getInfPos
    !>
    !>      infPos = getInfPos(mold = infPos)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isInf](@ref pm_except::isInf)<br>
    !>  [isNAN](@ref pm_except::isNAN)<br>
    !>  [getNAN](@ref pm_except::getNAN)<br>
    !>  [setNAN](@ref pm_except::setNAN)<br>
    !>  [isInfPos](@ref pm_except::isInfPos)<br>
    !>  [isInfNeg](@ref pm_except::isInfNeg)<br>
    !>  [getInfPos](@ref pm_except::getInfPos)<br>
    !>  [setInfPos](@ref pm_except::setInfPos)<br>
    !>  [getInfNeg](@ref pm_except::getInfNeg)<br>
    !>  [setInfNeg](@ref pm_except::setInfNeg)<br>
    !>
    !>  \example{getInfPos}
    !>  \include{lineno} example/pm_except/getInfPos/main.F90
    !>  \compilef{getInfPos}
    !>  \output{getInfPos}
    !>  \include{lineno} example/pm_except/getInfPos/main.out.F90
    !>
    !>  \test
    !>  [test_pm_except](@ref test_pm_except)
    !>
    !>  \final{getInfPos}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface getInfPos

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function getInfPos_CK5(mold) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfPos_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    :: mold
        complex(CKG)                :: infPos
    end function
#endif

#if CK4_ENABLED
    pure elemental module function getInfPos_CK4(mold) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfPos_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    :: mold
        complex(CKG)                :: infPos
    end function
#endif

#if CK3_ENABLED
    pure elemental module function getInfPos_CK3(mold) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfPos_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    :: mold
        complex(CKG)                :: infPos
    end function
#endif

#if CK2_ENABLED
    pure elemental module function getInfPos_CK2(mold) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfPos_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    :: mold
        complex(CKG)                :: infPos
    end function
#endif

#if CK1_ENABLED
    pure elemental module function getInfPos_CK1(mold) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfPos_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    :: mold
        complex(CKG)                :: infPos
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function getInfPos_RK5(mold) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfPos_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: mold
        real(RKG)                   :: infPos
    end function
#endif

#if RK4_ENABLED
    pure elemental module function getInfPos_RK4(mold) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfPos_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: mold
        real(RKG)                   :: infPos
    end function
#endif

#if RK3_ENABLED
    pure elemental module function getInfPos_RK3(mold) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfPos_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: mold
        real(RKG)                   :: infPos
    end function
#endif

#if RK2_ENABLED
    pure elemental module function getInfPos_RK2(mold) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfPos_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: mold
        real(RKG)                   :: infPos
    end function
#endif

#if RK1_ENABLED
    pure elemental module function getInfPos_RK1(mold) result(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfPos_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: mold
        real(RKG)                   :: infPos
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return an IEEE-compliant positive infinity.<br>
    !>
    !>  \details
    !>  This procedure is useful for deliberate initializations of `real` or `complex` numbers to positive infinity
    !>  without signaling the occurrence of an exception.
    !>
    !>  \param[out] infPos  :   The output scalar or array of arbitrary rank of either<br>
    !>                          <ol>
    !>                              <li>    type `complex` of kind \CKALL, or <br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ol>
    !>                          which will contain a positive infinity.
    !>
    !>  \interface{setInfPos}
    !>  \code{.F90}
    !>
    !>      use pm_except, only: setInfPos
    !>
    !>      call setInfPos(infPos)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isInf](@ref pm_except::isInf)<br>
    !>  [isNAN](@ref pm_except::isNAN)<br>
    !>  [getNAN](@ref pm_except::getNAN)<br>
    !>  [setNAN](@ref pm_except::setNAN)<br>
    !>  [isInfPos](@ref pm_except::isInfPos)<br>
    !>  [isInfNeg](@ref pm_except::isInfNeg)<br>
    !>  [getInfPos](@ref pm_except::getInfPos)<br>
    !>  [setInfPos](@ref pm_except::setInfPos)<br>
    !>  [getInfNeg](@ref pm_except::getInfNeg)<br>
    !>  [setInfNeg](@ref pm_except::setInfNeg)<br>
    !>
    !>  \example{setInfPos}
    !>  \include{lineno} example/pm_except/setInfPos/main.F90
    !>  \compilef{setInfPos}
    !>  \output{setInfPos}
    !>  \include{lineno} example/pm_except/setInfPos/main.out.F90
    !>
    !>  \test
    !>  [test_pm_except](@ref test_pm_except)
    !>
    !>  \final{setInfPos}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface setInfPos

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module subroutine setInfPos_CK5(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfPos_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(out)   :: infPos
    end subroutine
#endif

#if CK4_ENABLED
    pure elemental module subroutine setInfPos_CK4(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfPos_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(out)   :: infPos
    end subroutine
#endif

#if CK3_ENABLED
    pure elemental module subroutine setInfPos_CK3(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfPos_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(out)   :: infPos
    end subroutine
#endif

#if CK2_ENABLED
    pure elemental module subroutine setInfPos_CK2(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfPos_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(out)   :: infPos
    end subroutine
#endif

#if CK1_ENABLED
    pure elemental module subroutine setInfPos_CK1(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfPos_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(out)   :: infPos
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module subroutine setInfPos_RK5(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfPos_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)   :: infPos
    end subroutine
#endif

#if RK4_ENABLED
    pure elemental module subroutine setInfPos_RK4(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfPos_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)   :: infPos
    end subroutine
#endif

#if RK3_ENABLED
    pure elemental module subroutine setInfPos_RK3(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfPos_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)   :: infPos
    end subroutine
#endif

#if RK2_ENABLED
    pure elemental module subroutine setInfPos_RK2(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfPos_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)   :: infPos
    end subroutine
#endif

#if RK1_ENABLED
    pure elemental module subroutine setInfPos_RK1(infPos)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfPos_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)   :: infPos
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input value is an **IEEE-compliant** negative infinity.<br>
    !>
    !>  \details
    !>  If the input value is a `complex` number, then the output is `.true.` if any of the two real or
    !>  imaginary components or both are negative infinities.
    !>
    !>  \param[in]  x   :   The input scalar or array of arbitrary rank of either <br>
    !>                      <ol>
    !>                          <li>    type `complex` of kind \CKALL, or <br>
    !>                          <li>    type `real` of kind \RKALL, <br>
    !>                      </ol>
    !>                      whose value will be tested for being negative infinity.
    !>
    !>  \return
    !>  `infNeg`        :   The output scalar or array of the same shape as the input `x` of type `logical` of default kind \LK
    !>                      whose value is `.true.` if the input `x` is a negative infinity, otherwise it is `.false.`.
    !>
    !>  \interface{isInfNeg}
    !>  \code{.F90}
    !>
    !>      use pm_except, only: isInfNeg
    !>      use pm_kind, only: LK
    !>      logical(LK) :: infNeg
    !>
    !>      infNeg = isInfNeg(x)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \note
    !>  A very simple test of negative infinity of a variable `x` is the condition `x < -huge(x)` that is `.true.` if `x` is a negative infinity.
    !>
    !>  \note
    !>  The procedures under this generic interface are equivalent to `ieee_is_negative(x) .and. .not. ieee_is_finite(x)`
    !>  from the `ieee_arithmetic` intrinsic module to detect negative infinity.<br>
    !>  This generic interface extends this function also to `complex` numbers.<br>
    !>
    !>  \see
    !>  [isInf](@ref pm_except::isInf)<br>
    !>  [isNAN](@ref pm_except::isNAN)<br>
    !>  [getNAN](@ref pm_except::getNAN)<br>
    !>  [setNAN](@ref pm_except::setNAN)<br>
    !>  [isInfPos](@ref pm_except::isInfPos)<br>
    !>  [isInfNeg](@ref pm_except::isInfNeg)<br>
    !>  [getInfPos](@ref pm_except::getInfPos)<br>
    !>  [setInfPos](@ref pm_except::setInfPos)<br>
    !>  [getInfNeg](@ref pm_except::getInfNeg)<br>
    !>  [setInfNeg](@ref pm_except::setInfNeg)<br>
    !>
    !>  \example{isInfNeg}
    !>  \include{lineno} example/pm_except/isInfNeg/main.F90
    !>  \compilef{isInfNeg}
    !>  \output{isInfNeg}
    !>  \include{lineno} example/pm_except/isInfNeg/main.out.F90
    !>
    !>  \test
    !>  [test_pm_except](@ref test_pm_except)
    !>
    !>  \final{isInfNeg}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface isInfNeg

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isInfNeg_CK5(x) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfNeg_CK5
#endif
        use pm_kind, only: LK, CKG => CK5
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: infNeg
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isInfNeg_CK4(x) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfNeg_CK4
#endif
        use pm_kind, only: LK, CKG => CK4
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: infNeg
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isInfNeg_CK3(x) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfNeg_CK3
#endif
        use pm_kind, only: LK, CKG => CK3
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: infNeg
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isInfNeg_CK2(x) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfNeg_CK2
#endif
        use pm_kind, only: LK, CKG => CK2
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: infNeg
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isInfNeg_CK1(x) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfNeg_CK1
#endif
        use pm_kind, only: LK, CKG => CK1
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: infNeg
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function isInfNeg_RK5(x) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfNeg_RK5
#endif
        use pm_kind, only: LK, RKG => RK5
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: infNeg
    end function
#endif

#if RK4_ENABLED
    pure elemental module function isInfNeg_RK4(x) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfNeg_RK4
#endif
        use pm_kind, only: LK, RKG => RK4
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: infNeg
    end function
#endif

#if RK3_ENABLED
    pure elemental module function isInfNeg_RK3(x) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfNeg_RK3
#endif
        use pm_kind, only: LK, RKG => RK3
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: infNeg
    end function
#endif

#if RK2_ENABLED
    pure elemental module function isInfNeg_RK2(x) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfNeg_RK2
#endif
        use pm_kind, only: LK, RKG => RK2
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: infNeg
    end function
#endif

#if RK1_ENABLED
    pure elemental module function isInfNeg_RK1(x) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInfNeg_RK1
#endif
        use pm_kind, only: LK, RKG => RK1
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: infNeg
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an IEEE-compliant negative infinity.
    !>
    !>  \details
    !>  This procedure is useful for deliberate initializations of `real` or `complex` numbers to negative infinity
    !>  without signaling the occurrence of an exception.
    !>
    !>  \param[out] mold    :   The input scalar or array of arbitrary rank of either<br>
    !>                          <ol>
    !>                              <li>    type `complex` of kind \CKALL, or <br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ol>
    !>                          whose type and kind will match the generated output `-Inf`.
    !>
    !>  \return
    !>  `infNeg`            :   The output scalar or array of the same type and kind as
    !>                          the input `mold` containing the generated negative infinity value.
    !>
    !>  \interface{getInfNeg}
    !>  \code{.F90}
    !>
    !>      use pm_except, only: getInfNeg
    !>
    !>      infNeg = getInfNeg(mold = infNeg)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isInf](@ref pm_except::isInf)<br>
    !>  [isNAN](@ref pm_except::isNAN)<br>
    !>  [getNAN](@ref pm_except::getNAN)<br>
    !>  [setNAN](@ref pm_except::setNAN)<br>
    !>  [isInfPos](@ref pm_except::isInfPos)<br>
    !>  [isInfNeg](@ref pm_except::isInfNeg)<br>
    !>  [getInfPos](@ref pm_except::getInfPos)<br>
    !>  [setInfPos](@ref pm_except::setInfPos)<br>
    !>  [getInfNeg](@ref pm_except::getInfNeg)<br>
    !>  [setInfNeg](@ref pm_except::setInfNeg)<br>
    !>
    !>  \example{getInfNeg}
    !>  \include{lineno} example/pm_except/getInfNeg/main.F90
    !>  \compilef{getInfNeg}
    !>  \output{getInfNeg}
    !>  \include{lineno} example/pm_except/getInfNeg/main.out.F90
    !>
    !>  \test
    !>  [test_pm_except](@ref test_pm_except)
    !>
    !>  \final{getInfNeg}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface getInfNeg

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function getInfNeg_CK5(mold) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfNeg_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    :: mold
        complex(CKG)                :: infNeg
    end function
#endif

#if CK4_ENABLED
    pure elemental module function getInfNeg_CK4(mold) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfNeg_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    :: mold
        complex(CKG)                :: infNeg
    end function
#endif

#if CK3_ENABLED
    pure elemental module function getInfNeg_CK3(mold) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfNeg_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    :: mold
        complex(CKG)                :: infNeg
    end function
#endif

#if CK2_ENABLED
    pure elemental module function getInfNeg_CK2(mold) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfNeg_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    :: mold
        complex(CKG)                :: infNeg
    end function
#endif

#if CK1_ENABLED
    pure elemental module function getInfNeg_CK1(mold) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfNeg_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    :: mold
        complex(CKG)                :: infNeg
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function getInfNeg_RK5(mold) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfNeg_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: mold
        real(RKG)                   :: infNeg
    end function
#endif

#if RK4_ENABLED
    pure elemental module function getInfNeg_RK4(mold) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfNeg_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: mold
        real(RKG)                   :: infNeg
    end function
#endif

#if RK3_ENABLED
    pure elemental module function getInfNeg_RK3(mold) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfNeg_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: mold
        real(RKG)                   :: infNeg
    end function
#endif

#if RK2_ENABLED
    pure elemental module function getInfNeg_RK2(mold) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfNeg_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: mold
        real(RKG)                   :: infNeg
    end function
#endif

#if RK1_ENABLED
    pure elemental module function getInfNeg_RK1(mold) result(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getInfNeg_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: mold
        real(RKG)                   :: infNeg
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return an IEEE-compliant negative infinity.<br>
    !>
    !>  \brief
    !>  This procedure is useful for deliberate initializations of `real` or `complex` numbers to negative infinity
    !>  without signaling the occurrence of an exception.
    !>
    !>  \param[out] infNeg  :   The output scalar or array of arbitrary rank of either<br>
    !>                          <ol>
    !>                              <li>    type `complex` of kind \CKALL, or <br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ol>
    !>                          which will contain a negative infinity.
    !>
    !>  \interface{setInfNeg}
    !>  \code{.F90}
    !>
    !>      use pm_except, only: setInfNeg
    !>
    !>      call setInfNeg(infNeg)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isInf](@ref pm_except::isInf)<br>
    !>  [isNAN](@ref pm_except::isNAN)<br>
    !>  [getNAN](@ref pm_except::getNAN)<br>
    !>  [setNAN](@ref pm_except::setNAN)<br>
    !>  [isInfPos](@ref pm_except::isInfPos)<br>
    !>  [isInfNeg](@ref pm_except::isInfNeg)<br>
    !>  [getInfPos](@ref pm_except::getInfPos)<br>
    !>  [setInfPos](@ref pm_except::setInfPos)<br>
    !>  [getInfNeg](@ref pm_except::getInfNeg)<br>
    !>  [setInfNeg](@ref pm_except::setInfNeg)<br>
    !>
    !>  \example{setInfNeg}
    !>  \include{lineno} example/pm_except/setInfNeg/main.F90
    !>  \compilef{setInfNeg}
    !>  \output{setInfNeg}
    !>  \include{lineno} example/pm_except/setInfNeg/main.out.F90
    !>
    !>  \test
    !>  [test_pm_except](@ref test_pm_except)
    !>
    !>  \final{setInfNeg}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface setInfNeg

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module subroutine setInfNeg_CK5(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfNeg_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(out)   :: infNeg
    end subroutine
#endif

#if CK4_ENABLED
    pure elemental module subroutine setInfNeg_CK4(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfNeg_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(out)   :: infNeg
    end subroutine
#endif

#if CK3_ENABLED
    pure elemental module subroutine setInfNeg_CK3(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfNeg_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(out)   :: infNeg
    end subroutine
#endif

#if CK2_ENABLED
    pure elemental module subroutine setInfNeg_CK2(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfNeg_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(out)   :: infNeg
    end subroutine
#endif

#if CK1_ENABLED
    pure elemental module subroutine setInfNeg_CK1(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfNeg_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(out)   :: infNeg
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module subroutine setInfNeg_RK5(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfNeg_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)   :: infNeg
    end subroutine
#endif

#if RK4_ENABLED
    pure elemental module subroutine setInfNeg_RK4(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfNeg_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)   :: infNeg
    end subroutine
#endif

#if RK3_ENABLED
    pure elemental module subroutine setInfNeg_RK3(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfNeg_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)   :: infNeg
    end subroutine
#endif

#if RK2_ENABLED
    pure elemental module subroutine setInfNeg_RK2(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfNeg_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)   :: infNeg
    end subroutine
#endif

#if RK1_ENABLED
    pure elemental module subroutine setInfNeg_RK1(infNeg)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setInfNeg_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)   :: infNeg
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if the input value is an **IEEE-compliant** `NAN` (Not a Number)
    !>  or if the input value `x` is not equal to its input copy `xcopy`, a characteristic behavior of `NAN` values.<br>
    !>
    !>  \details
    !>  If the input value is a `complex` number, then the output is `.true.`
    !>  if any of the real or imaginary components or both components are `NaN`.
    !>
    !>  \param[in]  x       :   The input scalar or array of the same rank as other input array-like arguments of either<br>
    !>                          <ol>
    !>                              <li>    type `complex` of kind \CKALL, or <br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ol>
    !>                          whose value will be tested for being `NaN`.
    !>  \param[in]  xcopy   :   The input scalar or array of the same rank as other input array-like arguments, 
    !>                          of the same type and kind and **value** as the input argument `x`.<br>
    !>                          <ol>
    !>                              <li>    When `xcopy` is missing, the procedure relies on IEEE-compliant intrinsic routines for detecting `NaN`.
    !>                              <li>    When `xcopy` is present, its value is directly compared with the input `x` for non-equality to detect `NaN`.<br>
    !>                          </ol>
    !>                          (**optional**. If missing, IEEE-compliance is assumed for detecting `NaN`.)
    !>
    !>  \return
    !>  `isNotANumber`      :   The output scalar or array of the same shape as the input `x` of type `logical` of default kind \LK
    !>                          whose value is `.true.` if the input `x` is `NaN`, otherwise it is `.false.`.
    !>
    !>  \interface{isNAN}
    !>  \code{.F90}
    !>
    !>      use pm_kind, only: LK
    !>      use pm_except, only: isNAN
    !>      logical(LK) :: isNotANumber
    !>
    !>      isNotANumber = isNAN(x)
    !>      isNotANumber = isNAN(x, xcopy)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Keep in mind that `NaN` might not even be defined for non-IEEE-compliant compilers
    !>  or when certain aggressive IEEE-breaking optimization flags are specified like
    !>  `-O3 -ffast-math` with GNU `gfortran` or `-fast` with Intel `ifort` compilers.<br>
    !>  In such cases, this algorithm with only one input argument which relies on IEEE-compliance may spectacularly fail.<br>
    !>  After all, why would one need `NaN` in fast mode?<br>
    !>  The two-arguments interface relies on the observation that `NaN` is the only value that is not equal to itself.<br>
    !>  This approach to detecting `NaN` works as long as the compiler does not perform aggressive optimizations
    !>  to inline the `NaN`-detection procedure, thereby allowing to eliminate the direct
    !>  comparison of `x` with itself through aggressive optimizations.<br>
    !>  As such, the procedures under the two-arguments interface are declared `impure` to present inlining of the procedures.<br>
    !>
    !>  \remark
    !>  The procedures of this generic interface are `pure` with one input argument and `impure` with two input arguments.<br>
    !>
    !>  \elemental
    !>
    !>  \note
    !>  According to the IEEE definition, if a `real` variable `x` is `NaN` then `x /= x` always yields `.true.`.<br>
    !>  This can be simple fast test of `NaN` when the processor is IEEE-compliant.<br>
    !>  However, keep in mind that some compilers might optimize such comparison away.<br>
    !>
    !>  \note
    !>  The procedures under this generic interface with single input argument `x` 
    !>  use `ieee_is_nan(x)` from the `ieee_arithmetic` intrinsic module to detect `NaN` values.<br>
    !>  This generic interface extends this function also to `complex` numbers.<br>
    !>
    !>  \note
    !>  Note that the Intel and GNU Fortran compilers also have the `isnan(x)` extension function.
    !>
    !>  \see
    !>  [isInf](@ref pm_except::isInf)<br>
    !>  [isNAN](@ref pm_except::isNAN)<br>
    !>  [getNAN](@ref pm_except::getNAN)<br>
    !>  [setNAN](@ref pm_except::setNAN)<br>
    !>  [isInfPos](@ref pm_except::isInfPos)<br>
    !>  [isInfNeg](@ref pm_except::isInfNeg)<br>
    !>  [getInfPos](@ref pm_except::getInfPos)<br>
    !>  [setInfPos](@ref pm_except::setInfPos)<br>
    !>  [getInfNeg](@ref pm_except::getInfNeg)<br>
    !>  [setInfNeg](@ref pm_except::setInfNeg)<br>
    !>
    !>  \example{isNAN}
    !>  \include{lineno} example/pm_except/isNAN/main.F90
    !>  \compilef{isNAN}
    !>  \output{isNAN}
    !>  \include{lineno} example/pm_except/isNAN/main.out.F90
    !>
    !>  \test
    !>  [test_pm_except](@ref test_pm_except)
    !>
    !>  \todo
    !>  \pvhigh
    !>  The implementation of the test for `NaN` should be improved to
    !>  a bitwise comparison that is also valid with non-IEEE-compliant processors.
    !>
    !>  \final{isNAN}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface isNAN

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function isNANIEEE_CK5(x) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANIEEE_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: isNotANumber
    end function
#endif

#if CK4_ENABLED
    pure elemental module function isNANIEEE_CK4(x) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANIEEE_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: isNotANumber
    end function
#endif

#if CK3_ENABLED
    pure elemental module function isNANIEEE_CK3(x) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANIEEE_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: isNotANumber
    end function
#endif

#if CK2_ENABLED
    pure elemental module function isNANIEEE_CK2(x) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANIEEE_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: isNotANumber
    end function
#endif

#if CK1_ENABLED
    pure elemental module function isNANIEEE_CK1(x) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANIEEE_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    :: x
        logical(LK)                 :: isNotANumber
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function isNANIEEE_RK5(x) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANIEEE_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: isNotANumber
    end function
#endif

#if RK4_ENABLED
    pure elemental module function isNANIEEE_RK4(x) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANIEEE_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: isNotANumber
    end function
#endif

#if RK3_ENABLED
    pure elemental module function isNANIEEE_RK3(x) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANIEEE_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: isNotANumber
    end function
#endif

#if RK2_ENABLED
    pure elemental module function isNANIEEE_RK2(x) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANIEEE_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: isNotANumber
    end function
#endif

#if RK1_ENABLED
    pure elemental module function isNANIEEE_RK1(x) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANIEEE_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: x
        logical(LK)                 :: isNotANumber
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure elemental module function isNANXNEQ_CK5(x, xcopy) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANXNEQ_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    :: x, xcopy
        logical(LK)                 :: isNotANumber
    end function
#endif

#if CK4_ENABLED
    impure elemental module function isNANXNEQ_CK4(x, xcopy) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANXNEQ_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    :: x, xcopy
        logical(LK)                 :: isNotANumber
    end function
#endif

#if CK3_ENABLED
    impure elemental module function isNANXNEQ_CK3(x, xcopy) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANXNEQ_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    :: x, xcopy
        logical(LK)                 :: isNotANumber
    end function
#endif

#if CK2_ENABLED
    impure elemental module function isNANXNEQ_CK2(x, xcopy) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANXNEQ_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    :: x, xcopy
        logical(LK)                 :: isNotANumber
    end function
#endif

#if CK1_ENABLED
    impure elemental module function isNANXNEQ_CK1(x, xcopy) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANXNEQ_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    :: x, xcopy
        logical(LK)                 :: isNotANumber
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function isNANXNEQ_RK5(x, xcopy) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANXNEQ_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: x, xcopy
        logical(LK)                 :: isNotANumber
    end function
#endif

#if RK4_ENABLED
    impure elemental module function isNANXNEQ_RK4(x, xcopy) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANXNEQ_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: x, xcopy
        logical(LK)                 :: isNotANumber
    end function
#endif

#if RK3_ENABLED
    impure elemental module function isNANXNEQ_RK3(x, xcopy) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANXNEQ_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: x, xcopy
        logical(LK)                 :: isNotANumber
    end function
#endif

#if RK2_ENABLED
    impure elemental module function isNANXNEQ_RK2(x, xcopy) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANXNEQ_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: x, xcopy
        logical(LK)                 :: isNotANumber
    end function
#endif

#if RK1_ENABLED
    impure elemental module function isNANXNEQ_RK1(x, xcopy) result(isNotANumber)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isNANXNEQ_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: x, xcopy
        logical(LK)                 :: isNotANumber
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return an IEEE-compliant quiet `NAN` (Not a Number).<br>
    !>
    !>  \details
    !>  This procedure is useful for deliberate initializations of `real` or `complex` numbers to `NaN`
    !>  without signaling the occurrence of an exception.
    !>
    !>  \param[out] mold    :   The input scalar or array of arbitrary rank of either<br>
    !>                          <ol>
    !>                              <li>    type `complex` of kind \CKALL, or <br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ol>
    !>                          whose type and kind will match the generated output `NaN`.
    !>
    !>  \return
    !>  `nan`               :   The output scalar or array of the same type and kind as
    !>                          the input `mold` containing the generated `NaN` value.
    !>
    !>  \interface{getNAN}
    !>  \code{.F90}
    !>
    !>      use pm_except, only: getNAN
    !>
    !>      nan = getNAN(mold = nan)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Keep in mind that `NaN` might not even be defined for non-IEEE-compliant compilers
    !>  or when certain aggressive IEEE-breaking optimization flags are specified like
    !>  `-O3 -ffast-math` with GNU `gfortran` or `-fast` with Intel `ifort` compilers.
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isInf](@ref pm_except::isInf)<br>
    !>  [isNAN](@ref pm_except::isNAN)<br>
    !>  [getNAN](@ref pm_except::getNAN)<br>
    !>  [setNAN](@ref pm_except::setNAN)<br>
    !>  [isInfPos](@ref pm_except::isInfPos)<br>
    !>  [isInfNeg](@ref pm_except::isInfNeg)<br>
    !>  [getInfPos](@ref pm_except::getInfPos)<br>
    !>  [setInfPos](@ref pm_except::setInfPos)<br>
    !>  [getInfNeg](@ref pm_except::getInfNeg)<br>
    !>  [setInfNeg](@ref pm_except::setInfNeg)<br>
    !>
    !>  \example{getNAN}
    !>  \include{lineno} example/pm_except/getNAN/main.F90
    !>  \compilef{getNAN}
    !>  \output{getNAN}
    !>  \include{lineno} example/pm_except/getNAN/main.out.F90
    !>
    !>  \test
    !>  [test_pm_except](@ref test_pm_except)
    !>
    !>  \final{getNAN}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface getNAN

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module function getNAN_CK5(mold) result(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNAN_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(in)    :: mold
        complex(CKG)                :: nan
    end function
#endif

#if CK4_ENABLED
    pure elemental module function getNAN_CK4(mold) result(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNAN_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(in)    :: mold
        complex(CKG)                :: nan
    end function
#endif

#if CK3_ENABLED
    pure elemental module function getNAN_CK3(mold) result(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNAN_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(in)    :: mold
        complex(CKG)                :: nan
    end function
#endif

#if CK2_ENABLED
    pure elemental module function getNAN_CK2(mold) result(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNAN_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(in)    :: mold
        complex(CKG)                :: nan
    end function
#endif

#if CK1_ENABLED
    pure elemental module function getNAN_CK1(mold) result(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNAN_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(in)    :: mold
        complex(CKG)                :: nan
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module function getNAN_RK5(mold) result(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNAN_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(in)    :: mold
        real(RKG)                   :: nan
    end function
#endif

#if RK4_ENABLED
    pure elemental module function getNAN_RK4(mold) result(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNAN_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(in)    :: mold
        real(RKG)                   :: nan
    end function
#endif

#if RK3_ENABLED
    pure elemental module function getNAN_RK3(mold) result(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNAN_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(in)    :: mold
        real(RKG)                   :: nan
    end function
#endif

#if RK2_ENABLED
    pure elemental module function getNAN_RK2(mold) result(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNAN_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(in)    :: mold
        real(RKG)                   :: nan
    end function
#endif

#if RK1_ENABLED
    pure elemental module function getNAN_RK1(mold) result(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNAN_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(in)    :: mold
        real(RKG)                   :: nan
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return an IEEE-compliant quiet `NAN` (Not a Number).<br>
    !>
    !>  \details
    !>  This procedure is useful for deliberate initializations of `real` or `complex` numbers to `NaN`
    !>  without signaling the occurrence of an exception.
    !>
    !>  \param[out] nan :   The output scalar or array of arbitrary rank of either<br>
    !>                      <ol>
    !>                          <li>    type `complex` of kind \CKALL, or <br>
    !>                          <li>    type `real` of kind \RKALL, <br>
    !>                      </ol>
    !>                      which will contain a `NaN`.
    !>
    !>  \interface{setNAN}
    !>  \code{.F90}
    !>
    !>      use pm_except, only: setNAN
    !>
    !>      call setNAN(nan)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  Keep in mind that `NaN` might not even be defined for non-IEEE-compliant compilers
    !>  or when certain aggressive IEEE-breaking optimization flags are specified like
    !>  `-O3 -ffast-math` with GNU `gfortran` or `-fast` with Intel `ifort` compilers.
    !>
    !>  \pure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isInf](@ref pm_except::isInf)<br>
    !>  [isNAN](@ref pm_except::isNAN)<br>
    !>  [getNAN](@ref pm_except::getNAN)<br>
    !>  [setNAN](@ref pm_except::setNAN)<br>
    !>  [isInfPos](@ref pm_except::isInfPos)<br>
    !>  [isInfNeg](@ref pm_except::isInfNeg)<br>
    !>  [getInfPos](@ref pm_except::getInfPos)<br>
    !>  [setInfPos](@ref pm_except::setInfPos)<br>
    !>  [getInfNeg](@ref pm_except::getInfNeg)<br>
    !>  [setInfNeg](@ref pm_except::setInfNeg)<br>
    !>
    !>  \example{setNAN}
    !>  \include{lineno} example/pm_except/setNAN/main.F90
    !>  \compilef{setNAN}
    !>  \output{setNAN}
    !>  \include{lineno} example/pm_except/setNAN/main.out.F90
    !>
    !>  \test
    !>  [test_pm_except](@ref test_pm_except)
    !>
    !>  \final{setNAN}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface setNAN

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    pure elemental module subroutine setNAN_CK5(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNAN_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG), intent(out)   :: nan
    end subroutine
#endif

#if CK4_ENABLED
    pure elemental module subroutine setNAN_CK4(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNAN_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG), intent(out)   :: nan
    end subroutine
#endif

#if CK3_ENABLED
    pure elemental module subroutine setNAN_CK3(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNAN_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG), intent(out)   :: nan
    end subroutine
#endif

#if CK2_ENABLED
    pure elemental module subroutine setNAN_CK2(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNAN_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG), intent(out)   :: nan
    end subroutine
#endif

#if CK1_ENABLED
    pure elemental module subroutine setNAN_CK1(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNAN_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG), intent(out)   :: nan
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    pure elemental module subroutine setNAN_RK5(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNAN_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)   , intent(out)   :: nan
    end subroutine
#endif

#if RK4_ENABLED
    pure elemental module subroutine setNAN_RK4(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNAN_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)   , intent(out)   :: nan
    end subroutine
#endif

#if RK3_ENABLED
    pure elemental module subroutine setNAN_RK3(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNAN_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)   , intent(out)   :: nan
    end subroutine
#endif

#if RK2_ENABLED
    pure elemental module subroutine setNAN_RK2(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNAN_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)   , intent(out)   :: nan
    end subroutine
#endif

#if RK1_ENABLED
    pure elemental module subroutine setNAN_RK1(nan)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setNAN_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)   , intent(out)   :: nan
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_except