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
!>  This module contains procedures and generic interfaces for computing `log(1 + x)` <b>more precisely for tiny `x`</b>.
!>
!>  \details
!>  This module is particularly useful for computing `log(1 + x)` when `x` is smaller than the computer precision of the desired `real` kind.<br>
!>  This precision is the value returned by the Fortran intrinsic `epsilon()` procedure.<br>
!>
!>  \see
!>  [getLog1p](@ref pm_mathLog1p::getLog1p)<br>
!>  [get1mexp](@ref pm_math1mexp::get1mexp)<br>
!>  [getLogAddExp](@ref pm_mathLogAddExp::getLogAddExp)<br>
!>  [getLogSubExp](@ref pm_mathLogSubExp::getLogSubExp)<br>
!>  [getLogSumExp](@ref pm_mathLogSumExp::getLogSumExp)<br>
!>
!>  \finmain
!>
!>  \test
!>  [test_pm_mathLog1p](@ref test_pm_mathLog1p)<br>
!>
!>  \author
!>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathLog1p

    use pm_kind, only: SK, IK, LK
    use pm_control, only: selection_type, selection

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathLog1p"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the expression `log(1 + x)` robustly (without numerical underflow).
    !>
    !>  \param[in]  x           :   The input scalar, or array of arbitrary rank, shape, and size, of either <br>
    !>                              <ul>
    !>                                  <li>    type `complex` of kind \CKALL, or<br>
    !>                                  <li>    type `real` of kind \RKALL, <br>
    !>                              </ul>
    !>                              representing the `x` value whose `log(1 + x)` is to be returned.<br>
    !                               INTERNAL NOTE: Based on the benchmark results,
    !                               there does not appear to exist any performance benefits in having a control optional argument.
    !                               As such, it is removed from the interface.
    !   \param      control     :   The input scalar object that can be,
    !                               <ol>
    !                                   <li>    the constant [selection](@ref pm_control::selection) or equivalently,
    !                                           an object of type [selection_type](@ref pm_control::selection_type).<br>
    !                                           Specifying this option enables the runtime checks for underflow occurrence via branching and dynamic dispatch.<br>
    !                                           Enabling this option can aid runtime efficiency when the `x` ratio is expected to cause underflow.<br>
    !                                           In such cases, the logarithm can be avoided if `control =` [selection](@ref pm_control::selection), leading
    !                                           to better runtime efficiency since logarithm is highly expensive (on the order of ~200 CPU cycles).<br>
    !                                           See [the relevant benchmark here](#benchmark-getLog1p).<br>
    !                               </ol>
    !                               The presence of this argument is merely for compile-time resolution of the procedures of this generic interface.<br>
    !                               (**optional**. If missing, a sequence control flow will be assumed without dynamic dispatch.)
    !>
    !>  \return
    !>  `log1p`                 :   The output scalar (or array) of the same type and kind (and shape) as the input `x` representing `log(1 + x)` without underflow.
    !>
    !>  \interface
    !>  \code{.F90}
    !>
    !>      use pm_mathLog1p, only: getLog1p, selection
    !>
    !>      log1p = getLog1p(x)
    !>      log1p = getLog1p(x, control)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `x > -1.` must hold.<br>
    !>  The condition `x < huge(x)` must hold.<br>
    !>  When the input arguments are of type `complex`, these conditions must hold for the `real` components of the numbers.<br>
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
    !>  \example
    !>  \include{lineno} example/pm_mathLog1p/getLog1p/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_mathLog1p/getLog1p/main.out.F90
    !>
    !  \benchmarks
    !
    !  \benchmark{getLog1p, The effects of `control` on runtime efficiency}
    !      The following program compares the runtime performance of [getLog1p](@ref pm_mathLog1p::getLog1p)
    !      algorithm with and without checking for underflows.
    !  \include{lineno} benchmark/pm_mathLog1p/getLog1p/main.F90
    !  \compilefb{getLog1p}
    !  \postprocb{getLog1p}
    !  \include{lineno} benchmark/pm_mathLog1p/getLog1p/main.py
    !  \visb{getLog1p}
    !  \image html benchmark/pm_mathLog1p/getLog1p/benchmark.getLog1p.runtime.png width=1000
    !  \image html benchmark/pm_mathLog1p/getLog1p/benchmark.getLog1p.runtime.ratio.png width=1000
    !  \moralb{getLog1p}
    !      -#  If the input value `x` is known to be frequently `< tiny(x)` or `> epsilon(x)` then it is generally beneficial
    !          to call the interface with `control` argument set to [selection](@ref pm_control::selection).<br>
    !>
    !>  \test
    !>  [test_pm_mathLog1p](@ref test_pm_mathLog1p)<br>
    !>
    !>  \finmain
    !>
    !>  \author
    !>  \AmirShahmoradi, Thursday 1:45 AM, August 22, 2019, Dallas, TX
    interface getLog1p

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE elemental module function getLog1pSeq_CK5(x) result(log1p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSeq_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)        , intent(in)    :: x
        complex(CKC)                        :: log1p
    end function
#endif

#if CK4_ENABLED
    PURE elemental module function getLog1pSeq_CK4(x) result(log1p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSeq_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC), intent(in)    :: x
        complex(CKC)                :: log1p
    end function
#endif

#if CK3_ENABLED
    PURE elemental module function getLog1pSeq_CK3(x) result(log1p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSeq_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC), intent(in)    :: x
        complex(CKC)                :: log1p
    end function
#endif

#if CK2_ENABLED
    PURE elemental module function getLog1pSeq_CK2(x) result(log1p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSeq_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC), intent(in)    :: x
        complex(CKC)                :: log1p
    end function
#endif

#if CK1_ENABLED
    PURE elemental module function getLog1pSeq_CK1(x) result(log1p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSeq_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC), intent(in)    :: x
        complex(CKC)                :: log1p
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLog1pSeq_RK5(x) result(log1p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSeq_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)    :: x
        real(RKC)                   :: log1p
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLog1pSeq_RK4(x) result(log1p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSeq_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)    :: x
        real(RKC)                   :: log1p
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLog1pSeq_RK3(x) result(log1p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSeq_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)    :: x
        real(RKC)                   :: log1p
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLog1pSeq_RK2(x) result(log1p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSeq_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)    :: x
        real(RKC)                   :: log1p
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLog1pSeq_RK1(x) result(log1p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSeq_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)    :: x
        real(RKC)                   :: log1p
    end function
#endif

!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if CK5_ENABLED
!    PURE elemental module function getLog1pSel_CK5(x, control) result(log1p)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSel_CK5
!#endif
!        use pm_kind, only: CKC => CK5
!        type(selection_type), intent(in)    :: control
!        complex(CKC)        , intent(in)    :: x
!        complex(CKC)                        :: log1p
!    end function
!#endif
!
!#if CK4_ENABLED
!    PURE elemental module function getLog1pSel_CK4(x, control) result(log1p)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSel_CK4
!#endif
!        use pm_kind, only: CKC => CK4
!        type(selection_type), intent(in)    :: control
!        complex(CKC)        , intent(in)    :: x
!        complex(CKC)                        :: log1p
!    end function
!#endif
!
!#if CK3_ENABLED
!    PURE elemental module function getLog1pSel_CK3(x, control) result(log1p)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSel_CK3
!#endif
!        use pm_kind, only: CKC => CK3
!        type(selection_type), intent(in)    :: control
!        complex(CKC)        , intent(in)    :: x
!        complex(CKC)                        :: log1p
!    end function
!#endif
!
!#if CK2_ENABLED
!    PURE elemental module function getLog1pSel_CK2(x, control) result(log1p)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSel_CK2
!#endif
!        use pm_kind, only: CKC => CK2
!        type(selection_type), intent(in)    :: control
!        complex(CKC)        , intent(in)    :: x
!        complex(CKC)                        :: log1p
!    end function
!#endif
!
!#if CK1_ENABLED
!    PURE elemental module function getLog1pSel_CK1(x, control) result(log1p)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSel_CK1
!#endif
!        use pm_kind, only: CKC => CK1
!        type(selection_type), intent(in)    :: control
!        complex(CKC)        , intent(in)    :: x
!        complex(CKC)                        :: log1p
!    end function
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE elemental module function getLog1pSel_RK5(x, control) result(log1p)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSel_RK5
!#endif
!        use pm_kind, only: RKC => RK5
!        type(selection_type), intent(in)    :: control
!        real(RKC)           , intent(in)    :: x
!        real(RKC)                           :: log1p
!    end function
!#endif
!
!#if RK4_ENABLED
!    PURE elemental module function getLog1pSel_RK4(x, control) result(log1p)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSel_RK4
!#endif
!        use pm_kind, only: RKC => RK4
!        type(selection_type), intent(in)    :: control
!        real(RKC)           , intent(in)    :: x
!        real(RKC)                           :: log1p
!    end function
!#endif
!
!#if RK3_ENABLED
!    PURE elemental module function getLog1pSel_RK3(x, control) result(log1p)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSel_RK3
!#endif
!        use pm_kind, only: RKC => RK3
!        type(selection_type), intent(in)    :: control
!        real(RKC)           , intent(in)    :: x
!        real(RKC)                           :: log1p
!    end function
!#endif
!
!#if RK2_ENABLED
!    PURE elemental module function getLog1pSel_RK2(x, control) result(log1p)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSel_RK2
!#endif
!        use pm_kind, only: RKC => RK2
!        type(selection_type), intent(in)    :: control
!        real(RKC)           , intent(in)    :: x
!        real(RKC)                           :: log1p
!    end function
!#endif
!
!#if RK1_ENABLED
!    PURE elemental module function getLog1pSel_RK1(x, control) result(log1p)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getLog1pSel_RK1
!#endif
!        use pm_kind, only: RKC => RK1
!        type(selection_type), intent(in)    :: control
!        real(RKC)           , intent(in)    :: x
!        real(RKC)                           :: log1p
!    end function
!#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathLog1p ! LCOV_EXCL_LINE