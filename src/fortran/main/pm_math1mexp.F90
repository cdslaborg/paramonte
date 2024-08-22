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
!>  This module contains procedures and generic interfaces for computing `1 - exp(x)` <b>more precisely for tiny `x`</b>.
!>
!>  \details
!>  This module is particularly useful for computing `1 - exp(x)` when `x` is near or smaller than the computer precision of the desired `real` kind.<br>
!>  This precision is the value returned by the Fortran intrinsic `epsilon()` procedure.<br>
!>
!>  \see
!>  [Beebe, 2002, Computation of expm1(x) = exp(x) − 1](https://scholar.google.com/scholar?cluster=14303983886882776474&hl=en&as_sdt=0,44&scioq=Beebe,+2002,+Computation+of+expm1(x)+%3D+exp(x)+%E2%88%92+1)<br>
!>  [getLog1p](@ref pm_mathLog1p::getLog1p)<br>
!>  [get1mexp](@ref pm_math1mexp::get1mexp)<br>
!>  [getLogAddExp](@ref pm_mathLogAddExp::getLogAddExp)<br>
!>  [getLogSubExp](@ref pm_mathLogSubExp::getLogSubExp)<br>
!>  [getLogSumExp](@ref pm_mathLogSumExp::getLogSumExp)<br>
!>
!>  \final
!>
!>  \test
!>  [test_pm_math1mexp](@ref test_pm_math1mexp)
!>
!>  \author
!>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_math1mexp

    use pm_kind, only: SK, IK, LK
    use pm_control, only: selection_type, selection

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_math1mexp"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the expression `1 - exp(x)` robustly (without numerical underflow).
    !>
    !>  \param[in]  x           :   The input scalar, or array of arbitrary rank, shape, and size, of either <br>
    !>                              <ul>
    !>                                  <li>    type `complex` of kind \CKALL, or<br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              representing the `x` value whose `1 - exp(x)` is to be returned.<br>
    !>  \param[in]  control     :   The input scalar object that can be,
    !>                              <ol>
    !>                                  <li>    the constant [selection](@ref pm_control::selection) or equivalently, 
    !>                                          an object of type [selection_type](@ref pm_control::selection_type).<br>
    !>                                          Specifying this option enables the runtime checks for underflow occurrence via branching and dynamic dispatch.<br>
    !>                                          **Enabling this option can aid runtime efficiency when the `x` ratio is expected to cause underflow**.<br>
    !>                                          This option avoids the computation of a logarithm term, leading to better runtime efficiency.<br>
    !>                                          Note that logarithm is highly expensive (on the order of ~200 CPU cycles).<br>
    !>                                          See [the relevant benchmark here](#benchmark-get1mexp).<br>
    !>                              </ol>
    !>                              The presence of this argument is merely for compile-time resolution of the procedures of this generic interface.<br>
    !>                              (**optional**. If missing, a sequence control flow will be assumed without dynamic dispatch.)
    !>
    !>  \return
    !>  `onemexp`               :   The output scalar (or array) of the same type and kind (and shape) as the input `x` representing `1 - exp(x)` without underflow.
    !>
    !>  \interface{get1mexp}
    !>  \code{.F90}
    !>
    !>      use pm_math1mexp, only: get1mexp, selection
    !>
    !>      onemexp = get1mexp(x)
    !>      onemexp = get1mexp(x, control = selection)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `x < log(huge(x))` must hold.<br>
    !>  When the input arguments are of type `complex`, this condition must hold for the `real` components of the numbers.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getLog1p](@ref pm_mathLog1p::getLog1p)<br>
    !>  [get1mexp](@ref pm_math1mexp::get1mexp)<br>
    !>  [getLogAddExp](@ref pm_mathLogAddExp::getLogAddExp)<br>
    !>  [getLogSubExp](@ref pm_mathLogSubExp::getLogSubExp)<br>
    !>  [getLogSumExp](@ref pm_mathLogSumExp::getLogSumExp)<br>
    !>
    !>  \example{get1mexp}
    !>  \include{lineno} example/pm_math1mexp/get1mexp/main.F90
    !>  \compilef{get1mexp}
    !>  \output{get1mexp}
    !>  \include{lineno} example/pm_math1mexp/get1mexp/main.out.F90
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{get1mexp, The effects of `control` on runtime efficiency}
    !>      The following program compares the runtime performance of [get1mexp](@ref pm_math1mexp::get1mexp)
    !>      algorithm with and without checking for underflows.
    !>  \include{lineno} benchmark/pm_math1mexp/get1mexp/main.F90
    !>  \compilefb{get1mexp}
    !>  \postprocb{get1mexp}
    !>  \include{lineno} benchmark/pm_math1mexp/get1mexp/main.py
    !>  \visb{get1mexp}
    !>  \image html benchmark/pm_math1mexp/get1mexp/benchmark.get1mexp.runtime.png width=1000
    !>  \image html benchmark/pm_math1mexp/get1mexp/benchmark.get1mexp.runtime.ratio.png width=1000
    !>  \moralb{get1mexp}
    !>      -#  If the input value `x` is known to be smaller than `-log(huge)` then it is generally beneficial 
    !>          to call the interface with `control` argument set to [selection](@ref pm_control::selection).<br>
    !>
    !>  \test
    !>  [test_pm_math1mexp](@ref test_pm_math1mexp)
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX
    interface get1mexp

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module function get1mexpSeq_CK5(x) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSeq_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in)    :: x
        complex(CKG)                        :: onemexp
    end function
#endif

#if CK4_ENABLED
    PURE elemental module function get1mexpSeq_CK4(x) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSeq_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in)    :: x
        complex(CKG)                        :: onemexp
    end function
#endif

#if CK3_ENABLED
    PURE elemental module function get1mexpSeq_CK3(x) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSeq_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in)    :: x
        complex(CKG)                        :: onemexp
    end function
#endif

#if CK2_ENABLED
    PURE elemental module function get1mexpSeq_CK2(x) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSeq_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in)    :: x
        complex(CKG)                        :: onemexp
    end function
#endif

#if CK1_ENABLED
    PURE elemental module function get1mexpSeq_CK1(x) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSeq_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in)    :: x
        complex(CKG)                        :: onemexp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function get1mexpSeq_RK5(x) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSeq_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    :: x
        real(RKG)                           :: onemexp
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function get1mexpSeq_RK4(x) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSeq_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    :: x
        real(RKG)                           :: onemexp
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function get1mexpSeq_RK3(x) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSeq_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    :: x
        real(RKG)                           :: onemexp
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function get1mexpSeq_RK2(x) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSeq_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    :: x
        real(RKG)                           :: onemexp
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function get1mexpSeq_RK1(x) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSeq_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    :: x
        real(RKG)                           :: onemexp
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module function get1mexpSel_CK5(x, control) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSel_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in)    :: x
        complex(CKG)                        :: onemexp
        type(selection_type), intent(in)    :: control
    end function
#endif

#if CK4_ENABLED
    PURE elemental module function get1mexpSel_CK4(x, control) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSel_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in)    :: x
        complex(CKG)                        :: onemexp
        type(selection_type), intent(in)    :: control
    end function
#endif

#if CK3_ENABLED
    PURE elemental module function get1mexpSel_CK3(x, control) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSel_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in)    :: x
        complex(CKG)                        :: onemexp
        type(selection_type), intent(in)    :: control
    end function
#endif

#if CK2_ENABLED
    PURE elemental module function get1mexpSel_CK2(x, control) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSel_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in)    :: x
        complex(CKG)                        :: onemexp
        type(selection_type), intent(in)    :: control
    end function
#endif

#if CK1_ENABLED
    PURE elemental module function get1mexpSel_CK1(x, control) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSel_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in)    :: x
        complex(CKG)                        :: onemexp
        type(selection_type), intent(in)    :: control
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function get1mexpSel_RK5(x, control) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSel_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    :: x
        real(RKG)                           :: onemexp
        type(selection_type), intent(in)    :: control
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function get1mexpSel_RK4(x, control) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSel_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    :: x
        real(RKG)                           :: onemexp
        type(selection_type), intent(in)    :: control
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function get1mexpSel_RK3(x, control) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSel_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    :: x
        real(RKG)                           :: onemexp
        type(selection_type), intent(in)    :: control
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function get1mexpSel_RK2(x, control) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSel_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    :: x
        real(RKG)                           :: onemexp
        type(selection_type), intent(in)    :: control
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function get1mexpSel_RK1(x, control) result(onemexp)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: get1mexpSel_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    :: x
        real(RKG)                           :: onemexp
        type(selection_type), intent(in)    :: control
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_math1mexp ! LCOV_EXCL_LINE