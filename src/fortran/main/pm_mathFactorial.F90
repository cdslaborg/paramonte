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
!>  This module contains procedures and generic interfaces for the Factorial function.
!>
!>  \details
!>  This module provides a method [getFactorial](@ref pm_mathFactorial::getFactorial) to calculate the factorial of a number.
!>  However, note that factorial can easily overflow the integer representations of computers and even real numbers.
!>  As such, a better safer method of computing the factorial is to compute its natural logarithm via
!>  [getLogFactorial](@ref pm_mathFactorial::getLogFactorial).
!>
!>  \benchmarks
!>
!>  \benchmark{getFactorial_vs_getLogFactorial, The runtime performance of [getFactorial](@ref pm_mathFactorial::getFactorial) vs. [getLogFactorial](@ref pm_mathFactorial::getLogFactorial)}
!>  \include{lineno} benchmark/pm_mathFactorial/getFactorial_vs_getLogFactorial/main.F90
!>  \compilefb{getFactorial_vs_getLogFactorial}
!>  \postprocb{getFactorial_vs_getLogFactorial}
!>  \include{lineno} benchmark/pm_mathFactorial/getFactorial_vs_getLogFactorial/main.py
!>  \visb{getFactorial_vs_getLogFactorial}
!>  \image html benchmark/pm_mathFactorial/getFactorial_vs_getLogFactorial/benchmark.getFactorial_vs_getLogFactorial.runtime.png width=1000
!>  \image html benchmark/pm_mathFactorial/getFactorial_vs_getLogFactorial/benchmark.getFactorial_vs_getLogFactorial.runtime.ratio.png width=1000
!>  \moralb{getFactorial_vs_getLogFactorial}
!>      -#  The procedures under the generic interface [getFactorial](@ref pm_mathFactorial::getFactorial) compute the factorial
!>          using its default definition while the procedures under the generic interface [getLogFactorial](@ref pm_mathFactorial::getLogFactorial)
!>          use the Fortran intrinsic `log_gamma()` to compute the `log(factorial)` which is then converted to `factorial` in the benchmark code.<br>
!>          Based on the benchmark results, the safe method of computing the factorial as a real number (thus avoiding the potential numerical overflow)
!>          is about 3-5 times slower than the direct definition of the factorial, although the performance gap closes at large input whole numbers.<br>
!>
!>  \test
!>  [test_pm_mathFactorial](@ref test_pm_mathFactorial)
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

module pm_mathFactorial

    use pm_kind, only: IK, RK, SK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathFactorial"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the factorial of the input positive integer.
    !>
    !>  \details
    !>  The factorial of an integer number \f$n\f$ is defined as,
    !>
    !>  \f{equation}{
    !>      \large
    !>      n! = \prod_{i=1}^{n} i = \Gamma(n+1) ~.
    !>  \f}
    !>
    !>  Note that the factorial of a number can readily overflow the maximum integer values representable by computers.
    !>  As such, [getLogFactorial](@ref pm_mathFactorial::getLogFactorial) is a safer alternative to use.
    !>
    !>  \param[in]  n   :   The input scalar or array of arbitrary rank of type `integer` of kind \IKALL
    !>                      containing the non-negative integer whose \f$\log(n!)\f$ is to be computed on return.
    !>
    !>  \return
    !>  `factorial`     :   The output scalar or array of the same shape as the input `n` representing the factorial of `n`.
    !>
    !>  \interface{getFactorial}
    !>  \code{.F90}
    !>
    !>      use pm_mathFactorial, only: getFactorial
    !>
    !>      factorial = getFactorial(n)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0_IKC <= n` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getFactorial](@ref pm_mathFactorial::getFactorial)<br>
    !>
    !>  \example{getFactorial}
    !>  \include{lineno} example/pm_mathFactorial/getFactorial/main.F90
    !>  \compilef{getFactorial}
    !>  \output{getFactorial}
    !>  \include{lineno} example/pm_mathFactorial/getFactorial/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathFactorial](@ref test_pm_mathFactorial)
    !>
    !>  \finmain{getFactorial}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface getFactorial

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE elemental module function getFactorial_IK5(n) result(factorial)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFactorial_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC), intent(in)    :: n
        integer(IKC)                :: factorial
    end function
#endif

#if IK4_ENABLED
    PURE elemental module function getFactorial_IK4(n) result(factorial)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFactorial_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC), intent(in)    :: n
        integer(IKC)                :: factorial
    end function
#endif

#if IK3_ENABLED
    PURE elemental module function getFactorial_IK3(n) result(factorial)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFactorial_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC), intent(in)    :: n
        integer(IKC)                :: factorial
    end function
#endif

#if IK2_ENABLED
    PURE elemental module function getFactorial_IK2(n) result(factorial)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFactorial_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC), intent(in)    :: n
        integer(IKC)                :: factorial
    end function
#endif

#if IK1_ENABLED
    PURE elemental module function getFactorial_IK1(n) result(factorial)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFactorial_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC), intent(in)    :: n
        integer(IKC)                :: factorial
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the **natural logarithm** of the factorial of the input positive **whole** real number.
    !>
    !>  \details
    !>  The factorial of an integer number \f$n\f$ is defined as,
    !>
    !>  \f{equation}{
    !>      \large
    !>      n! = \prod_{i=1}^{n} i = \Gamma(n+1) ~.
    !>  \f}
    !>
    !>  Therefore,
    !>
    !>  \f{equation}{
    !>      \large
    !>      \log(n!) = \sum_{i=1}^{n} i = \log\bigg(\Gamma(n+1)\bigg) ~.
    !>  \f}
    !>
    !>  Note that the factorial of a number can readily overflow the maximum integer or even real values representable by computers.
    !>  This is the primary reason for computing the natural logarithm of the factorial instead of the factorial.
    !>
    !>  \param[in]  x   :   The input scalar or array of arbitrary rank of type `real` of kind \RKALL containing
    !>                      the whole number (integer) whose \f$\log(x!)\f$ is to be computed on return.
    !>
    !>  \return
    !>  `logFactorial`  :   The output scalar or array of the same shape as the input `x`
    !>                      representing the natural logarithm of the factorial of `x`.
    !>
    !>  \interface{getLogFactorial}
    !>  \code{.F90}
    !>
    !>      use pm_mathFactorial, only: getLogFactorial
    !>
    !>      logFactorial = getLogFactorial(x)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The input `x` must be a positive whole number.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [getFactorial](@ref pm_mathFactorial::getFactorial)<br>
    !>
    !>  \example{getLogFactorial}
    !>  \include{lineno} example/pm_mathFactorial/getLogFactorial/main.F90
    !>  \compilef{getLogFactorial}
    !>  \output{getLogFactorial}
    !>  \include{lineno} example/pm_mathFactorial/getLogFactorial/main.out.F90
    !>  \postproc{getLogFactorial}
    !>  \include{lineno} example/pm_mathFactorial/getLogFactorial/main.py
    !>  \vis{getLogFactorial}
    !>  \image html pm_mathFactorial/getLogFactorial/getLogFactorial.RK.png width=700
    !>
    !>  \test
    !>  [test_pm_mathFactorial](@ref test_pm_mathFactorial)
    !>
    !>  \finmain{getLogFactorial}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface getLogFactorial

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module function getLogFactorial_RK5(x) result(logFactorial)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogFactorial_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)    :: x
        real(RKC)                   :: logFactorial
    end function
#endif

#if RK4_ENABLED
    PURE elemental module function getLogFactorial_RK4(x) result(logFactorial)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogFactorial_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)    :: x
        real(RKC)                   :: logFactorial
    end function
#endif

#if RK3_ENABLED
    PURE elemental module function getLogFactorial_RK3(x) result(logFactorial)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogFactorial_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)    :: x
        real(RKC)                   :: logFactorial
    end function
#endif

#if RK2_ENABLED
    PURE elemental module function getLogFactorial_RK2(x) result(logFactorial)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogFactorial_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)    :: x
        real(RKC)                   :: logFactorial
    end function
#endif

#if RK1_ENABLED
    PURE elemental module function getLogFactorial_RK1(x) result(logFactorial)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogFactorial_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)    :: x
        real(RKC)                   :: logFactorial
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathFactorial ! LCOV_EXCL_LINE