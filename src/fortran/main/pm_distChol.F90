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
!>  This module contains classes and procedures for generating random upper or lower Cholesky factor triangular matrices.
!>
!>  \details
!>  The output random Cholesky factors can be subsequently used to generate random positive-definite matrices.<br>
!>  Note that every real positive definite matrix \f$M\f$ has a Cholesky factorization
!>  \f{equation}{
!>      M = LL*
!>  \f}
!>  where \f$L\f$ is a uniquely defined lower triangular matrix with positive diagonal entries.<br>
!>  Therefore, \f$M\f$ can be constructed from a given random \f$L\f$.<br>
!>  This approach, called **Gram method** is fast, however, the resulting matrix \f$M\f$ does not possess any particular structure.<br>
!>
!>  \see
!>  [pm_distCov](@ref pm_distCov)<br>
!>
!>  \test
!>  [test_pm_distChol](@ref test_pm_distChol)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_distChol

    use pm_matrixSubset, only: subset_type, uppDia_type, lowDia_type, uppDia, lowDia
    use pm_distUnif, only: rngf_type, xoshiro256ssw_type
    use pm_kind, only: SK, IK, LK

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_distChol"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return a random upper and lower Cholesky factorization.<br>
    !>
    !>  \brief
    !>  See the documentation of [pm_distChol](@ref pm_distChol) for details.<br>
    !>  See also [setCholRand](@ref pm_distChol::setCholRand) for generating only the upper or lower subset of the Cholesky factor.<br>
    !>
    !>  \param[in]  mold        :   The input scalar of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              whose type and kind determines the type and kind of the output `rand`.<br>
    !>                              The value of `mold` is ignored entirely within the algorithm.<br>
    !>  \param[in]  ndim        :   The input positive scalar of type `integer` of default kind \IK,
    !>                              representing the rank of the matrix (the number of dimensions) of shape `(ndim, ndim)`.<br>
    !>  \param[in]  subset      :   The input scalar constant that can be any of the following:<br>
    !>                              <ol>
    !>                                  <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) or an object of type [uppDia_type](@ref pm_matrixSubset::uppDia_type)
    !>                                          implying that the upper-diagonal triangular block of the argument `rand` must be used while the lower subset is not referenced.<br>
    !>                                  <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) or an object of type [lowDia_type](@ref pm_matrixSubset::lowDia_type)
    !>                                          implying that the lower-diagonal triangular block of the argument `rand` must be used while the upper subset is not referenced.<br>
    !>                              </ol>
    !>                              This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.<br>
    !>                              Beware that the oppsite subset will not be set to zero explicitly and it is the user's responsibility to ensure it.<br>
    !>                              The generic interface [getCholRand](@ref pm_distChol::getCholRand) can be used to generate
    !>                              random Cholesky factors whose complementary subset elements are all set to zero.<br>
    !>
    !>  \return
    !>  `rand`                  :   The output matrix of shape `(1:ndim, 1:ndim)` of the same type and kind as the input argument `mold`,
    !>                              containing an upper or lower Cholesky factorization as specified by the input `subset`.<br>
    !>
    !>  \interface{getCholRand}
    !>  \code{.F90}
    !>
    !>      use pm_distChol, only: getCholRand
    !>
    !>      rand(1:ndim, 1:ndim) = getCholRand(mold, ndim, subset)
    !>
    !>  \endcode
    !>
    !>  \see
    !>  [getCovRand](@ref pm_distCov::getCovRand)<br>
    !>  [setCovRand](@ref pm_distCov::setCovRand)<br>
    !>  [getCholRand](@ref pm_distChol::getCholRand)<br>
    !>  [setCholRand](@ref pm_distChol::setCholRand)<br>
    !>
    !>  \example{getCholRand}
    !>  \include{lineno} example/pm_distChol/getCholRand/main.F90
    !>  \compilef{getCholRand}
    !>  \output{getCholRand}
    !>  \include{lineno} example/pm_distChol/getCholRand/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distChol](@ref test_pm_distChol)
    !>
    !>  \final{getCholRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    interface getCholRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getCholRandRNGD_CK5(mold, ndim, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCholRandRNGD_CK5
#endif
        use pm_kind, only: TKG => CK5
        integer(IK)             , intent(in)                    :: ndim
        complex(TKG)            , intent(in)                    :: mold
        class(*)                , intent(in)                    :: subset
        complex(TKG)                                            :: rand(ndim, ndim)
    end function
#endif

#if CK4_ENABLED
    impure module function getCholRandRNGD_CK4(mold, ndim, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCholRandRNGD_CK4
#endif
        use pm_kind, only: TKG => CK4
        integer(IK)             , intent(in)                    :: ndim
        complex(TKG)            , intent(in)                    :: mold
        class(*)                , intent(in)                    :: subset
        complex(TKG)                                            :: rand(ndim, ndim)
    end function
#endif

#if CK3_ENABLED
    impure module function getCholRandRNGD_CK3(mold, ndim, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCholRandRNGD_CK3
#endif
        use pm_kind, only: TKG => CK3
        integer(IK)             , intent(in)                    :: ndim
        complex(TKG)            , intent(in)                    :: mold
        class(*)                , intent(in)                    :: subset
        complex(TKG)                                            :: rand(ndim, ndim)
    end function
#endif

#if CK2_ENABLED
    impure module function getCholRandRNGD_CK2(mold, ndim, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCholRandRNGD_CK2
#endif
        use pm_kind, only: TKG => CK2
        integer(IK)             , intent(in)                    :: ndim
        complex(TKG)            , intent(in)                    :: mold
        class(*)                , intent(in)                    :: subset
        complex(TKG)                                            :: rand(ndim, ndim)
    end function
#endif

#if CK1_ENABLED
    impure module function getCholRandRNGD_CK1(mold, ndim, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCholRandRNGD_CK1
#endif
        use pm_kind, only: TKG => CK1
        integer(IK)             , intent(in)                    :: ndim
        complex(TKG)            , intent(in)                    :: mold
        class(*)                , intent(in)                    :: subset
        complex(TKG)                                            :: rand(ndim, ndim)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getCholRandRNGD_RK5(mold, ndim, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCholRandRNGD_RK5
#endif
        use pm_kind, only: TKG => RK5
        integer(IK)             , intent(in)                    :: ndim
        real(TKG)               , intent(in)                    :: mold
        class(*)                , intent(in)                    :: subset
        real(TKG)                                               :: rand(ndim, ndim)
    end function
#endif

#if RK4_ENABLED
    impure module function getCholRandRNGD_RK4(mold, ndim, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCholRandRNGD_RK4
#endif
        use pm_kind, only: TKG => RK4
        integer(IK)             , intent(in)                    :: ndim
        real(TKG)               , intent(in)                    :: mold
        class(*)                , intent(in)                    :: subset
        real(TKG)                                               :: rand(ndim, ndim)
    end function
#endif

#if RK3_ENABLED
    impure module function getCholRandRNGD_RK3(mold, ndim, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCholRandRNGD_RK3
#endif
        use pm_kind, only: TKG => RK3
        integer(IK)             , intent(in)                    :: ndim
        real(TKG)               , intent(in)                    :: mold
        class(*)                , intent(in)                    :: subset
        real(TKG)                                               :: rand(ndim, ndim)
    end function
#endif

#if RK2_ENABLED
    impure module function getCholRandRNGD_RK2(mold, ndim, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCholRandRNGD_RK2
#endif
        use pm_kind, only: TKG => RK2
        integer(IK)             , intent(in)                    :: ndim
        real(TKG)               , intent(in)                    :: mold
        class(*)                , intent(in)                    :: subset
        real(TKG)                                               :: rand(ndim, ndim)
    end function
#endif

#if RK1_ENABLED
    impure module function getCholRandRNGD_RK1(mold, ndim, subset) result(rand)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCholRandRNGD_RK1
#endif
        use pm_kind, only: TKG => RK1
        integer(IK)             , intent(in)                    :: ndim
        real(TKG)               , intent(in)                    :: mold
        class(*)                , intent(in)                    :: subset
        real(TKG)                                               :: rand(ndim, ndim)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return a random upper or lower Cholesky factorization.<br>
    !>
    !>  \brief
    !>  See the documentation of [pm_distChol](@ref pm_distChol) for details.
    !>
    !>  \param[inout]   rng     :   The input/output scalar that can be an object of,
    !>                              <ol>
    !>                                  <li>    type [rngf_type](@ref pm_distUnif::rngf_type),
    !>                                          implying the use of intrinsic Fortran uniform RNG.<br>
    !>                                  <li>    type [xoshiro256ssw_type](@ref pm_distUnif::xoshiro256ssw_type),
    !>                                          implying the use of [xoshiro256**](https://prng.di.unimi.it/) uniform RNG.<br>
    !>                              </ol>
    !>  \param[out]     rand    :   The output matrix of shape `(1:ndim, 1:ndim)` of,<br>
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing a random (optionally power-law-distributed determinant) positive-definite matrix.<br>
    !>                              The output `rand` can of `complex` type **if and only if** the optional input argument `method` is missing.<br>
    !>  \param[in]      subset  :   The input scalar constant that can be any of the following:<br>
    !>                              <ol>
    !>                                  <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) or an object of type [uppDia_type](@ref pm_matrixSubset::uppDia_type)
    !>                                          implying that the upper-diagonal triangular block of the argument `rand` must be used while the lower subset is not referenced.<br>
    !>                                  <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) or an object of type [lowDia_type](@ref pm_matrixSubset::lowDia_type)
    !>                                          implying that the lower-diagonal triangular block of the argument `rand` must be used while the upper subset is not referenced.<br>
    !>                              </ol>
    !>                              This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.<br>
    !>                              Beware that the oppsite subset will not be set to zero explicitly and it is the user responsibility to ensure it.<br>
    !>                              If needed, the generic interface [setMatInit](@ref pm_matrixInit::setMatInit) can be used to set the complementary subset to zero.<br>
    !>                              The generic interface [getCholRand](@ref pm_distChol::getCholRand) can be used to generate
    !>                              random Cholesky factors whose complementary subset elements are all set to zero.<br>
    !>
    !>  \interface{setCholRand}
    !>  \code{.F90}
    !>
    !>      use pm_distChol, only: setCholRand
    !>
    !>      call setCholRand(rng, rand(1:ndim, 1:ndim), subset)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(rand, 1) == size(rand, 2)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \see
    !>  [getCovRand](@ref pm_distCov::getCovRand)<br>
    !>  [setCovRand](@ref pm_distCov::setCovRand)<br>
    !>  [getCholRand](@ref pm_distChol::getCholRand)<br>
    !>  [setCholRand](@ref pm_distChol::setCholRand)<br>
    !>
    !>  \example{setCholRand}
    !>  \include{lineno} example/pm_distChol/setCholRand/main.F90
    !>  \compilef{setCholRand}
    !>  \output{setCholRand}
    !>  \include{lineno} example/pm_distChol/setCholRand/main.out.F90
    !>
    !>  \test
    !>  [test_pm_distChol](@ref test_pm_distChol)
    !>
    !>  \final{setCholRand}
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>

    ! setCholRandG RK

    interface setCholRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    module subroutine setCholRandRNGF_RK5(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGF_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(rngf_type)         , intent(in)                    :: rng
        real(TKG)               , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    module subroutine setCholRandRNGF_RK4(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGF_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(rngf_type)         , intent(in)                    :: rng
        real(TKG)               , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    module subroutine setCholRandRNGF_RK3(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGF_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(rngf_type)         , intent(in)                    :: rng
        real(TKG)               , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    module subroutine setCholRandRNGF_RK2(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGF_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(rngf_type)         , intent(in)                    :: rng
        real(TKG)               , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    module subroutine setCholRandRNGF_RK1(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGF_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(rngf_type)         , intent(in)                    :: rng
        real(TKG)               , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setCholRandRNGX_RK5(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGX_RK5
#endif
        use pm_kind, only: TKG => RK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(TKG)               , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setCholRandRNGX_RK4(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGX_RK4
#endif
        use pm_kind, only: TKG => RK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(TKG)               , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setCholRandRNGX_RK3(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGX_RK3
#endif
        use pm_kind, only: TKG => RK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(TKG)               , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setCholRandRNGX_RK2(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGX_RK2
#endif
        use pm_kind, only: TKG => RK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(TKG)               , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setCholRandRNGX_RK1(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGX_RK1
#endif
        use pm_kind, only: TKG => RK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        real(TKG)               , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! setCholRandG CK

    interface setCholRand

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    module subroutine setCholRandRNGF_CK5(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGF_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(rngf_type)         , intent(in)                    :: rng
        complex(TKG)            , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    module subroutine setCholRandRNGF_CK4(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGF_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(rngf_type)         , intent(in)                    :: rng
        complex(TKG)            , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    module subroutine setCholRandRNGF_CK3(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGF_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(rngf_type)         , intent(in)                    :: rng
        complex(TKG)            , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    module subroutine setCholRandRNGF_CK2(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGF_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(rngf_type)         , intent(in)                    :: rng
        complex(TKG)            , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    module subroutine setCholRandRNGF_CK1(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGF_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(rngf_type)         , intent(in)                    :: rng
        complex(TKG)            , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setCholRandRNGX_CK5(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGX_CK5
#endif
        use pm_kind, only: TKG => CK5
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        complex(TKG)            , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setCholRandRNGX_CK4(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGX_CK4
#endif
        use pm_kind, only: TKG => CK4
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        complex(TKG)            , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setCholRandRNGX_CK3(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGX_CK3
#endif
        use pm_kind, only: TKG => CK3
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        complex(TKG)            , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setCholRandRNGX_CK2(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGX_CK2
#endif
        use pm_kind, only: TKG => CK2
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        complex(TKG)            , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setCholRandRNGX_CK1(rng, rand, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setCholRandRNGX_CK1
#endif
        use pm_kind, only: TKG => CK1
        type(xoshiro256ssw_type), intent(inout)                 :: rng
        complex(TKG)            , intent(out)   , contiguous    :: rand(:,:)
        class(*)                , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_distChol