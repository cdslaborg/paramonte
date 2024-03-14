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
!>  This module contains classes and procedures for generating **Bernoulli**-distributed random numbers.<br>
!>
!>  \details
!>  The **Bernoulli** distribution, named after Swiss mathematician Jacob Bernoulli,
!>  is the **discrete probability distribution** of a random variable which takes the value \f$1\f$ with probability \f$p\f$
!>  and the value \f$0\f$ with probability \f$q = 1-p\f$.<br>
!>
!>  Less formally, it can be thought of as a model for the set of possible outcomes of any single experiment that asks a yes–no question.<br>
!>  The Bernoulli distribution is a special case of the **Binomial distribution** with \f$n = 1\f$ (with \f$n\f$ the number of trials).<br>
!>  It is also sometimes called the **two-point distribution** because of having only two possible outcomes.<br>
!>
!>  \test
!>  [test_pm_distBern](@ref test_pm_distBern)
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distBern

    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = SK_"@pm_distBern"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a scalar (or a vector of length `size` or an array of the same shape as the input `p`)
    !>  containing the odds of getting a `head` in a single (or a series) of coin-flipping experiment(s).<br>
    !>
    !>  \details
    !>  The coin can be optionally biased toward getting a head (`.true.`) with the input user-specified probability `p`.<br>
    !>
    !>  \param[in]  p       :   The input scalar or array of arbitrary rank of type `real` of kind \RKALL within the range `[0, 1]`,
    !>                          representing the probability of success/yes/true (or getting *head* in a coin-flip).<br>
    !>                          A value other than `p = 0.5` corresponds to a biased coin-flip experiment.<br>
    !>                          <ul>
    !>                              <li>    If the input argument `size` is present, then `p` must be a scalar.<br>
    !>                              <li>    If the input argument `size` is missing, then `p` can be a scalar or array of arbitrary rank and shape.<br>
    !>                          </ul>
    !>                          (**optional**, default = `0.5`)
    !>  \param[in]  size    :   The input scalar of type `integer` of default kind \IK, representing the size of the output `rand` vector.<br>
    !>                          (**optional**, if present, then `p` must be a scalar and the output will be a vector.)
    !>
    !>  \return
    !>  `rand`              :   The output scalar or vector of length `size`, or array of the same shape as the input `p` of type `logical` of default kind \LK,
    !>                          containing the Bernoulli-distributed random output value.
    !>
    !>  \interface{isHead}
    !>  \code{.F90}
    !>
    !>      use pm_distBern, only: isHead
    !>
    !>      rand = isHead() ! single fair-coin flipping.
    !>      rand = isHead(p) ! elemental coin-flipping with `p` as the odds of success.
    !>      rand(1:size) = isHead(p, size) ! non-elemental coin-flipping of `size` times with `p` as the odds of success.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The input value for `p` must be a number between `0` and `1`.<br>
    !>  \vericon
    !>
    !>  \warning
    !>  Keep in mind that the interface `isHead()` uses a default single-precision value `p = .5`.<br>
    !>  On the current hardware technology, this corresponds to a 32-bit `real` precision.<br>
    !>  This implies a repeat cycle of \f$\sim10^6\f$ for the generated random numbers.<br>
    !>  If more random numbers are needed, specify the value `p = .5_RKC`,
    !>  where `RKC` refers to a higher real precision kind capable of
    !>  generating the desired number of random coin-flipping.<br>
    !>
    !>  \impure
    !>
    !>  \remark
    !>  The procedures under this generic interface are `elemental` when `p` is present and the `size` argument is missing.<br>
    !>
    !>  \note
    !>  The procedure [isHead()](@ref pm_distBern::isHead) without any input arguments is
    !>  equivalent to [getUnifRand()](@ref pm_distUnif::getUnifRand) without any input arguments.
    !>
    !>  \see
    !>  [getBernRand](@ref pm_distBern::getBernRand)<br>
    !>  [setBernRand](@ref pm_distBern::setBernRand)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>
    !>  \example{isHead}
    !>  \include{lineno} example/pm_distBern/isHead/main.F90
    !>  \compilef{isHead}
    !>  \output{isHead}
    !>  \include{lineno} example/pm_distBern/isHead/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distBern](@ref test_pm_distBern)
    !>
    !>  \finmain{isHead}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface isHead

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function isHeadDD() result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isHeadDD
#endif
        logical(LK)                                 :: rand
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    impure module function isHeadDS(size) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isHeadDS
#endif
        integer(IK)     , intent(in)                :: size
        logical(LK)                                 :: rand(size)
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function isHeadPD_RK5(p) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isHeadPD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                    :: p
        logical(LK)                                 :: rand
    end function
#endif

#if RK4_ENABLED
    impure elemental module function isHeadPD_RK4(p) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isHeadPD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                    :: p
        logical(LK)                                 :: rand
    end function
#endif

#if RK3_ENABLED
    impure elemental module function isHeadPD_RK3(p) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isHeadPD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                    :: p
        logical(LK)                                 :: rand
    end function
#endif

#if RK2_ENABLED
    impure elemental module function isHeadPD_RK2(p) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isHeadPD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                    :: p
        logical(LK)                                 :: rand
    end function
#endif

#if RK1_ENABLED
    impure elemental module function isHeadPD_RK1(p) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isHeadPD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                    :: p
        logical(LK)                                 :: rand
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function isHeadPS_RK5(p, size) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isHeadPS_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)       , intent(in)                :: p
        integer(IK)     , intent(in)                :: size
        logical(LK)                                 :: rand(size)
    end function
#endif

#if RK4_ENABLED
    impure module function isHeadPS_RK4(p, size) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isHeadPS_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)       , intent(in)                :: p
        integer(IK)     , intent(in)                :: size
        logical(LK)                                 :: rand(size)
    end function
#endif

#if RK3_ENABLED
    impure module function isHeadPS_RK3(p, size) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isHeadPS_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)       , intent(in)                :: p
        integer(IK)     , intent(in)                :: size
        logical(LK)                                 :: rand(size)
    end function
#endif

#if RK2_ENABLED
    impure module function isHeadPS_RK2(p, size) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isHeadPS_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)       , intent(in)                :: p
        integer(IK)     , intent(in)                :: size
        logical(LK)                                 :: rand(size)
    end function
#endif

#if RK1_ENABLED
    impure module function isHeadPS_RK1(p, size) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isHeadPS_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)       , intent(in)                :: p
        integer(IK)     , intent(in)                :: size
        logical(LK)                                 :: rand(size)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a scalar or array of rank `1` of length `size` or the same rank and size as `p` of
    !>  Bernoulli-distributed random values (`0` or `1`) with probability of getting `1` set by the input success probability `p`.<br>
    !>
    !>  \param[in]  p       :   The input scalar or array of arbitrary rank of type `real` of kind \RKALL, representing the probability of success/yes/true.
    !>                          It must be a number between `0` and `1`. It must be scalar if the input argument `size` is present.
    !>  \param[in]  size    :   The input scalar of type `integer` of default kind \IK, representing the size of the output `rand` vector.<br>
    !>                          (**optional**, if present, then `p` must be a scalar.)
    !>
    !>  \return
    !>  `rand`              :   The output scalar or vector of length `size`, or array of the same shape as the input `p` of type `integer` of default kind \IK,
    !>                          containing the Bernoulli-distributed random output value.
    !>
    !>  \interface{getBernRand}
    !>  \code{.F90}
    !>
    !>      use pm_distBern, only: getBernRand
    !>
    !>      rand = getBernRand(p) ! elemental
    !>      rand(1:size) = getBernRand(p, size) ! non-elemental
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The input value for `p` must be a number between `0` and `1`.<br>
    !>  \vericon
    !>
    !>  \impure
    !>
    !>  \elemental
    !>  The procedures under this generic interface are not `elemental` when the `size` argument is present.
    !>
    !>  \see
    !>  [isHead](@ref pm_distBern::isHead)<br>
    !>  [setBernRand](@ref pm_distBern::setBernRand)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>
    !>  \example{getBernRand}
    !>  \include{lineno} example/pm_distBern/getBernRand/main.F90
    !>  \compilef{getBernRand}
    !>  \output{getBernRand}
    !>  \include{lineno} example/pm_distBern/getBernRand/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distBern](@ref test_pm_distBern)
    !>
    !>  \finmain{getBernRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface getBernRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure elemental module function getBernRandPD_RK5(p) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBernRandPD_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(in)                    :: p
        integer(IK)                                 :: rand
    end function
#endif

#if RK4_ENABLED
    impure elemental module function getBernRandPD_RK4(p) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBernRandPD_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(in)                    :: p
        integer(IK)                                 :: rand
    end function
#endif

#if RK3_ENABLED
    impure elemental module function getBernRandPD_RK3(p) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBernRandPD_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(in)                    :: p
        integer(IK)                                 :: rand
    end function
#endif

#if RK2_ENABLED
    impure elemental module function getBernRandPD_RK2(p) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBernRandPD_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(in)                    :: p
        integer(IK)                                 :: rand
    end function
#endif

#if RK1_ENABLED
    impure elemental module function getBernRandPD_RK1(p) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBernRandPD_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(in)                    :: p
        integer(IK)                                 :: rand
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getBernRandPS_RK5(p, size) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBernRandPS_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)       , intent(in)                :: p
        integer(IK)     , intent(in)                :: size
        integer(IK)                                 :: rand(size)
    end function
#endif

#if RK4_ENABLED
    impure module function getBernRandPS_RK4(p, size) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBernRandPS_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)       , intent(in)                :: p
        integer(IK)     , intent(in)                :: size
        integer(IK)                                 :: rand(size)
    end function
#endif

#if RK3_ENABLED
    impure module function getBernRandPS_RK3(p, size) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBernRandPS_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)       , intent(in)                :: p
        integer(IK)     , intent(in)                :: size
        integer(IK)                                 :: rand(size)
    end function
#endif

#if RK2_ENABLED
    impure module function getBernRandPS_RK2(p, size) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBernRandPS_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)       , intent(in)                :: p
        integer(IK)     , intent(in)                :: size
        integer(IK)                                 :: rand(size)
    end function
#endif

#if RK1_ENABLED
    impure module function getBernRandPS_RK1(p, size) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBernRandPS_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)       , intent(in)                :: p
        integer(IK)     , intent(in)                :: size
        integer(IK)                                 :: rand(size)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a scalar or array of arbitrary rank of Bernoulli-distributed random values (`0` or `1`),
    !>  with the specified input success probability `p`.
    !>
    !>  \param[out] rand    :   The output scalar (or array of the same rank, shape, and size as other array-like arguments), of<br>
    !>                          <ul>
    !>                              <li>    type `integer` of kind \IKALL, <br>
    !>                              <li>    type `logical` of kind \LKALL, <br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ul>
    !>                          containing the Bernoulli-distributed random output value.<br>
    !>                          <ul>
    !>                              <li>    If `rand` is `real`, then an output value of `1.` means success/yes/true and `0.` means failure/no/false.<br>
    !>                              <li>    If `rand` is `integer`, then an output value of `1` means success/yes/true and `0` means failure/no/false.<br>
    !>                              <li>    If `rand` is `logical`, then an output value of `.true.` means success/yes/true and `.false.` means failure/no/false.<br>
    !>                          </ul>
    !>  \param[in]  urand   :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of<br>
    !>                          <ul>
    !>                              <li>    type `real` of the same kind as `rand` if `rand` is `real`, or<br>
    !>                              <li>    type `real` of the kind \RKALL, if `rand` is of type `logical` or `integer`,<br>
    !>                          </ul>
    !>                          containing uniformly-distributed random value(s) with the range `0 <= urand < 1`.<br>
    !>                          Such random value(s) can be readily obtained via the Fortran intrinsic procedure `random_number()` or via [getUnifRand()](@ref pm_distUnif::getUnifRand).<br>
    !>                          Supplying this argument ensures the purity of the procedures, allowing further compiler optimizations.<br>
    !>  \param[in]  p       :   The input scalar (or array of the same rank, shape, and size as other array-like arguments), of<br>
    !>                          <ul>
    !>                              <li>    type `real` of the same kind as `rand` if `rand` is `real`, or<br>
    !>                              <li>    type `real` of the kind \RKALL, if `rand` is of type `logical` or `integer`,<br>
    !>                          </ul>
    !>                          representing the success probability of in the Bernoulli trial.<br>
    !>                          Note that `p` must be a number in the range `[0,1]`.
    !>
    !>  \interface{setBernRand}
    !>  \code{.F90}
    !>
    !>      use pm_distBern, only: setBernRand
    !>
    !>      call setBernRand(rand, urand)
    !>      call setBernRand(rand, urand, p)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The conditions `0 <= p <= 1` must hold for the corresponding input arguments.<br>
    !>  The conditions `0 <= urand < 1` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \elemental
    !>
    !>  \see
    !>  [isHead](@ref pm_distBern::isHead)<br>
    !>  [getBernRand](@ref pm_distBern::getBernRand)<br>
    !>  [getUnifRand](@ref pm_distUnif::getUnifRand)<br>
    !>  [setUnifRand](@ref pm_distUnif::setUnifRand)<br>
    !>
    !>  \example{setBernRand}
    !>  \include{lineno} example/pm_distBern/setBernRand/main.F90
    !>  \compilef{setBernRand}
    !>  \output{setBernRand}
    !>  \include{lineno} example/pm_distBern/setBernRand/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distBern](@ref test_pm_distBern)
    !>
    !>  \finmain{setBernRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    interface setBernRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED && RK5_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK5_RK5(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK5_RK5
#endif
        use pm_kind, only: IKC => IK5, RKC => RK5
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK5_ENABLED && RK4_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK5_RK4(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK5_RK4
#endif
        use pm_kind, only: IKC => IK5, RKC => RK4
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK5_ENABLED && RK3_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK5_RK3(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK5_RK3
#endif
        use pm_kind, only: IKC => IK5, RKC => RK3
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK5_ENABLED && RK2_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK5_RK2(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK5_RK2
#endif
        use pm_kind, only: IKC => IK5, RKC => RK2
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK5_ENABLED && RK1_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK5_RK1(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK5_RK1
#endif
        use pm_kind, only: IKC => IK5, RKC => RK1
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK4_ENABLED && RK5_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK4_RK5(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK4_RK5
#endif
        use pm_kind, only: IKC => IK4, RKC => RK5
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK4_ENABLED && RK4_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK4_RK4(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK4_RK4
#endif
        use pm_kind, only: IKC => IK4, RKC => RK4
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK4_ENABLED && RK3_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK4_RK3(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK4_RK3
#endif
        use pm_kind, only: IKC => IK4, RKC => RK3
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK4_ENABLED && RK2_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK4_RK2(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK4_RK2
#endif
        use pm_kind, only: IKC => IK4, RKC => RK2
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK4_ENABLED && RK1_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK4_RK1(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK4_RK1
#endif
        use pm_kind, only: IKC => IK4, RKC => RK1
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK3_ENABLED && RK5_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK3_RK5(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK3_RK5
#endif
        use pm_kind, only: IKC => IK3, RKC => RK5
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK3_ENABLED && RK4_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK3_RK4(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK3_RK4
#endif
        use pm_kind, only: IKC => IK3, RKC => RK4
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK3_ENABLED && RK3_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK3_RK3(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK3_RK3
#endif
        use pm_kind, only: IKC => IK3, RKC => RK3
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK3_ENABLED && RK2_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK3_RK2(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK3_RK2
#endif
        use pm_kind, only: IKC => IK3, RKC => RK2
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK3_ENABLED && RK1_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK3_RK1(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK3_RK1
#endif
        use pm_kind, only: IKC => IK3, RKC => RK1
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK2_ENABLED && RK5_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK2_RK5(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK2_RK5
#endif
        use pm_kind, only: IKC => IK2, RKC => RK5
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK2_ENABLED && RK4_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK2_RK4(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK2_RK4
#endif
        use pm_kind, only: IKC => IK2, RKC => RK4
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK2_ENABLED && RK3_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK2_RK3(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK2_RK3
#endif
        use pm_kind, only: IKC => IK2, RKC => RK3
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK2_ENABLED && RK2_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK2_RK2(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK2_RK2
#endif
        use pm_kind, only: IKC => IK2, RKC => RK2
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK2_ENABLED && RK1_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK2_RK1(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK2_RK1
#endif
        use pm_kind, only: IKC => IK2, RKC => RK1
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK1_ENABLED && RK5_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK1_RK5(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK1_RK5
#endif
        use pm_kind, only: IKC => IK1, RKC => RK5
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK1_ENABLED && RK4_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK1_RK4(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK1_RK4
#endif
        use pm_kind, only: IKC => IK1, RKC => RK4
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK1_ENABLED && RK3_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK1_RK3(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK1_RK3
#endif
        use pm_kind, only: IKC => IK1, RKC => RK3
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK1_ENABLED && RK2_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK1_RK2(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK1_RK2
#endif
        use pm_kind, only: IKC => IK1, RKC => RK2
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if IK1_ENABLED && RK1_ENABLED
    PURE elemental module subroutine setBernRandRUP_IK1_RK1(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_IK1_RK1
#endif
        use pm_kind, only: IKC => IK1, RKC => RK1
        integer(IKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED && RK5_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK5_RK5(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK5_RK5
#endif
        use pm_kind, only: LKC => LK5, RKC => RK5
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK5_ENABLED && RK4_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK5_RK4(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK5_RK4
#endif
        use pm_kind, only: LKC => LK5, RKC => RK4
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK5_ENABLED && RK3_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK5_RK3(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK5_RK3
#endif
        use pm_kind, only: LKC => LK5, RKC => RK3
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK5_ENABLED && RK2_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK5_RK2(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK5_RK2
#endif
        use pm_kind, only: LKC => LK5, RKC => RK2
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK5_ENABLED && RK1_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK5_RK1(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK5_RK1
#endif
        use pm_kind, only: LKC => LK5, RKC => RK1
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK4_ENABLED && RK5_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK4_RK5(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK4_RK5
#endif
        use pm_kind, only: LKC => LK4, RKC => RK5
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK4_ENABLED && RK4_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK4_RK4(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK4_RK4
#endif
        use pm_kind, only: LKC => LK4, RKC => RK4
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK4_ENABLED && RK3_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK4_RK3(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK4_RK3
#endif
        use pm_kind, only: LKC => LK4, RKC => RK3
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK4_ENABLED && RK2_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK4_RK2(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK4_RK2
#endif
        use pm_kind, only: LKC => LK4, RKC => RK2
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK4_ENABLED && RK1_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK4_RK1(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK4_RK1
#endif
        use pm_kind, only: LKC => LK4, RKC => RK1
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK3_ENABLED && RK5_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK3_RK5(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK3_RK5
#endif
        use pm_kind, only: LKC => LK3, RKC => RK5
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK3_ENABLED && RK4_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK3_RK4(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK3_RK4
#endif
        use pm_kind, only: LKC => LK3, RKC => RK4
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK3_ENABLED && RK3_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK3_RK3(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK3_RK3
#endif
        use pm_kind, only: LKC => LK3, RKC => RK3
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK3_ENABLED && RK2_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK3_RK2(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK3_RK2
#endif
        use pm_kind, only: LKC => LK3, RKC => RK2
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK3_ENABLED && RK1_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK3_RK1(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK3_RK1
#endif
        use pm_kind, only: LKC => LK3, RKC => RK1
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK2_ENABLED && RK5_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK2_RK5(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK2_RK5
#endif
        use pm_kind, only: LKC => LK2, RKC => RK5
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK2_ENABLED && RK4_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK2_RK4(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK2_RK4
#endif
        use pm_kind, only: LKC => LK2, RKC => RK4
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK2_ENABLED && RK3_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK2_RK3(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK2_RK3
#endif
        use pm_kind, only: LKC => LK2, RKC => RK3
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK2_ENABLED && RK2_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK2_RK2(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK2_RK2
#endif
        use pm_kind, only: LKC => LK2, RKC => RK2
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK2_ENABLED && RK1_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK2_RK1(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK2_RK1
#endif
        use pm_kind, only: LKC => LK2, RKC => RK1
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK1_ENABLED && RK5_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK1_RK5(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK1_RK5
#endif
        use pm_kind, only: LKC => LK1, RKC => RK5
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK1_ENABLED && RK4_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK1_RK4(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK1_RK4
#endif
        use pm_kind, only: LKC => LK1, RKC => RK4
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK1_ENABLED && RK3_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK1_RK3(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK1_RK3
#endif
        use pm_kind, only: LKC => LK1, RKC => RK3
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK1_ENABLED && RK2_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK1_RK2(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK1_RK2
#endif
        use pm_kind, only: LKC => LK1, RKC => RK2
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if LK1_ENABLED && RK1_ENABLED
    PURE elemental module subroutine setBernRandRUP_LK1_RK1(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_LK1_RK1
#endif
        use pm_kind, only: LKC => LK1, RKC => RK1
        logical(LKC), intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE elemental module subroutine setBernRandRUP_RK5_RK5(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_RK5_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)   , intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if RK4_ENABLED
    PURE elemental module subroutine setBernRandRUP_RK4_RK4(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_RK4_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)   , intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if RK3_ENABLED
    PURE elemental module subroutine setBernRandRUP_RK3_RK3(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_RK3_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)   , intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if RK2_ENABLED
    PURE elemental module subroutine setBernRandRUP_RK2_RK2(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_RK2_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)   , intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

#if RK1_ENABLED
    PURE elemental module subroutine setBernRandRUP_RK1_RK1(rand, urand, p)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setBernRandRUP_RK1_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)   , intent(out)                   :: rand
        real(RKC)   , intent(in)                    :: urand, p
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distBern