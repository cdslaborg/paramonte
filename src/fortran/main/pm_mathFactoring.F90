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
!>  This module contains procedures and generic interfaces for computing the **prime factors** of integers.
!>
!>  \details
!>  In mathematics, **factorization** or **factoring** consists of writing a number or another mathematical object as a product of several factors,
!>  usually smaller or simpler objects of the same kind.<br>
!>  For example, \f$3 \times 5\f$ is a factorization of the integer \f$15\f$.<br>
!>
!>  The factoring of an integer number \f$n\f$ is defined such that,
!>
!>  \f{equation}{
!>      \large
!>      n = \prod_{i=1}^{m} f_i ~,
!>  \f}
!>
!>  where \f$m\f$ is the number of factors of \f$n\f$ and \f$f_i\f$ are the prime factors.<br>
!>
!>  \note
!>  Factorization of the lengths of data is an important prerequisite for computing
!>  the [Discrete Fourier Transform](@ref pm_fftpack) of `real` or `complex` data.<br>
!>
!>  [pm_fftpack](@ref pm_fftpack)<br>
!>  [pm_mathFactorial](@ref pm_mathFactorial)<br>
!>
!>  \test
!>  [test_pm_mathFactoring](@ref test_pm_mathFactoring)
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_mathFactoring

    use pm_kind, only: SK, IK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_mathFactoring"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the factoring of the input positive integer.
    !>
    !>  \details
    !>  The factoring of an integer number \f$n\f$ is defined such that,
    !>
    !>  \f{equation}{
    !>      \large
    !>      n = \prod_{i=1}^{m} f_i ~,
    !>  \f}
    !>
    !>  where \f$m\f$ is the number of factors of \f$n\f$ and \f$f_i\f$ are the factors.<br>
    !>
    !>  \param[in]  posint  :   The input scalar of type `integer` of kind \IKALL
    !>                          containing the positive integer (larger than `1`) whose factoring is to be computed on return.
    !>
    !>  \return
    !>  `Factoring`         :   The output `allocatable` vector of the same type and kind as the input `posint`,
    !>                          containing the **prime factors** of the input integer `posint`.<br>
    !>
    !>  \interface{getFactoring}
    !>  \code{.F90}
    !>
    !>      use pm_mathFactoring, only: getFactoring
    !>
    !>      Factoring = getFactoring(posint)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `1 < posint` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warning
    !>  The input argument `posint` has the `value` attribute.<br>
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getFactorial](@ref pm_mathFactorial::getFactorial)<br>
    !>
    !>  \example{getFactoring}
    !>  \include{lineno} example/pm_mathFactoring/getFactoring/main.F90
    !>  \compilef{getFactoring}
    !>  \output{getFactoring}
    !>  \include{lineno} example/pm_mathFactoring/getFactoring/main.out.F90
    !>
    !>  \test
    !>  [test_pm_mathFactoring](@ref test_pm_mathFactoring)
    !>
    !>  \todo
    !>  \phigh
    !>  The performance of the procedures under this generic interface can be improved by using more efficient
    !>  algorithms and lookup tables for various ranges of input values. The following schemes could be implemented,<br>
    !>  <ul>
    !>      <li>    \f$posint \leq 2^{16}\f$: Lookup table.<br>
    !>      <li>    \f$posint \leq 2^{70}\f$: Richard Brent modification of Pollard rho algorithm.
    !>      <li>    \f$posint \leq 10^{50}\f$: Lenstra elliptic curve factorization.
    !>      <li>    \f$posint \leq 10^{100}\f$: Quadratic Sieve.
    !>      <li>    \f$posint \leq 10^{100}\f$: General Number Field Sieve.
    !>  </ul>
    !>
    !>  \final{getFactoring}
    !>
    !>  \author
    !>  \AmirShahmoradi, April 23, 2017, 1:36 AM, Institute for Computational Engineering and Sciences (ICES), University of Texas at Austin
    interface getFactoring

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function getFactoring_IK5(posint) result(Factoring)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFactoring_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG), value                 :: posint
        integer(IKG), allocatable           :: Factoring(:)
    end function
#endif

#if IK4_ENABLED
    PURE module function getFactoring_IK4(posint) result(Factoring)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFactoring_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG), value                 :: posint
        integer(IKG), allocatable           :: Factoring(:)
    end function
#endif

#if IK3_ENABLED
    PURE module function getFactoring_IK3(posint) result(Factoring)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFactoring_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG), value                 :: posint
        integer(IKG), allocatable           :: Factoring(:)
    end function
#endif

#if IK2_ENABLED
    PURE module function getFactoring_IK2(posint) result(Factoring)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFactoring_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG), value                 :: posint
        integer(IKG), allocatable           :: Factoring(:)
    end function
#endif

#if IK1_ENABLED
    PURE module function getFactoring_IK1(posint) result(Factoring)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFactoring_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG), value                 :: posint
        integer(IKG), allocatable           :: Factoring(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_mathFactoring
