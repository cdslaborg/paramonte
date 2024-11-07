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
!>  This module contains procedures and generic interfaces for computing the Cholesky factorization of positive definite matrices.<br>
!>
!>  \details
!>
!>  Quick start
!>  -----------
!>
!>  The recommended routines for computing the Cholesky factorization are the procedures under
!>  the generic interface [setChoLow()](@ref setChoLow) with assumed-shape dummy arguments.<br>
!>  The explicit-shape procedures are not recommended as the array bounds cannot be checked and the return of error
!>  status for the Cholesky factorization failure is implicit (via the first element of the output argument `dia`).<br>
!>
!>  Cholesky factorization
!>  ----------------------
!>
!>  In linear algebra, the Cholesky decomposition or Cholesky factorization is a decomposition of a Hermitian,
!>  positive-definite matrix into the product of a lower triangular matrix and its conjugate transpose.<br>
!>  It was discovered by André-Louis Cholesky for real matrices, and posthumously published in 1924.<br>
!>
!>  Algorithmic efficiency
!>  ----------------------
!>
!>  When it is applicable, the Cholesky decomposition is roughly twice as efficient as
!>  the [LU decomposition](@ref pm_matrixLUP) for solving systems of linear equations.<br>
!>
!>  Definition
!>  ----------
!>
!>  The Cholesky decomposition of a Hermitian positive-definite matrix \f$A\f$ is a decomposition of the form,
!>  \f{equation}{
!>      \mathbf{A} = \mathbf{LL}^{*} ~,
!>  \f}
!>  where \f$L\f$ is a lower triangular matrix with real and positive diagonal entries, and \f$L*\f$ denotes the conjugate transpose of \f$L\f$.<br>
!>  Every Hermitian positive-definite matrix (and thus also every real-valued symmetric positive-definite matrix) has a unique Cholesky decomposition.<br>
!>  The converse holds trivially: if A can be written as \f$LL*\f$ for some invertible \f$L\f$, lower triangular or otherwise, then \f$A\f$ is Hermitian and positive definite.<br>
!>  When \f$A\f$ is a real matrix (hence symmetric positive-definite), the factorization may be written as,
!>  \f{equation}{
!>      \mathbf{A} = \mathbf{LL}^{\mathsf{T}} ~,
!>  \f}
!>  where \f$L\f$ is a real lower triangular matrix with positive diagonal entries.<br>
!>
!>  Algorithm
!>  ---------
!>
!>  For complex Hermitian matrix, the following formula applies<br>
!>  \f{eqnarray}{
!>      L_{j,j} &=& \sqrt{A_{j,j} - \sum_{k=1}^{j-1}L_{j,k}^{*}L_{j,k}} ~, \\
!>      L_{i,j} &=& \frac{1}{L_{j,j}} \left(A_{i,j} - \sum_{k=1}^{j-1}L_{j,k}^{*}L_{i,k}\right) \quad {\text{for }} i > j ~.
!>  \f}
!>  Therefore, the computation of the entry \f$(i, j)\f$ depends on the entries to the left and above.<br>
!>  The computation is usually arranged in either of the following orders:<br>
!>  <ol>
!>      <li>    The **Cholesky–Banachiewicz algorithm** (row-major) starts from the upper left corner of the matrix L and proceeds to calculate the matrix row by row:<br>
!>              \code{.F90}
!>                  do i = 1, size(A,1)
!>                      L(i,i) = sqrt(A(i,i) - dot_product(L(i,1:i-1), L(i,1:i-1)))
!>                      L(i+1:,i) = (A(i+1:,i) - matmul(conjg(L(i,1:i-1)), L(i+1:,1:i-1))) / L(i,i)
!>                  end do
!>              \endcode
!>      <li>    The **Cholesky–Crout algorithm** starts from the upper left corner of the matrix L and proceeds to calculate the matrix column by column:<br>
!>              \code{.F90}
!>                  do i = 1, size(A,1)
!>                      L(i,i) = sqrt(A(i,i) - dot_product(L(1:i-1,i), L(1:i-1,i)))
!>                      L(i,i+1:) = (A(i,i+1:) - matmul(conjg(L(1:i-1,i)), L(1:i-1,i+1:))) / L(i,i)
!>                  end do
!>              \endcode
!>  </ol>
!>  Either pattern of access allows the entire computation to be performed **in-place** if desired.<br>
!>  Given that Fortran is a column-major language, the **Cholesky–Crout algorithm** algorithm can potentially be faster on the contemporary hardware.<br>
!>  This means that the computations can be done more efficiently if one uses the upper-diagonal triangle of a
!>  Hermitian matrix to compute the corresponding upper-diagonal triangular Cholesky factorization.<br>
!>
!>  \benchmarks
!>
!>  \benchmark{cholesky_assumed_vs_explicit, Cholesky factorization - assumed-shape vs. explicit-shape dummy arguments}
!>  Here is a code snippet to compare the performances of the Cholesky factorization routine [setChoLow()](@ref setChoLow)
!>  via two different interfaces:
!>  <ol>
!>      <li>
!>          Assumed-shape dummy arguments with explicit return of failure error flag (**recommended**).
!>      </li>
!>      <li>
!>          Explicit-shape dummy arguments with implicit return of failure error flag (**not recommended**).
!>      </li>
!>  </ol>
!>  Overall, no significant difference between the two approaches is observed,
!>  meaning that the safe interface with assumed-shape arrays is much better to use.<br>
!>
!>  \include{lineno} benchmark/pm_matrixChol/cholesky_assumed_vs_explicit/main.F90
!>  \compilefb{cholesky_assumed_vs_explicit}
!>  \postprocb{cholesky_assumed_vs_explicit}
!>  \include{lineno} benchmark/pm_matrixChol/cholesky_assumed_vs_explicit/main.py
!>  \visb{cholesky_assumed_vs_explicit}
!>  \image html benchmark/pm_matrixChol/cholesky_assumed_vs_explicit/benchmark.cholesky_assumed_vs_explicit.runtime.png width=1000
!>  \image html benchmark/pm_matrixChol/cholesky_assumed_vs_explicit/benchmark.cholesky_assumed_vs_explicit.runtime.ratio.png width=1000
!>
!>  \todo
!>  \phigh
!>  A benchmark comparing the performance of the two computational algorithms above
!>  should be implemented and gauge the impact of row vs. column major access pattern.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!>  \cond
#if 1
#define _CONTIGUOUS_
#else
#define _CONTIGUOUS_, contiguous
#endif
!>  \endcond

module pm_matrixChol

    use pm_kind, only: SK, IK, LK
    use pm_control, only: iteration, iteration_type
    use pm_control, only: recursion, recursion_type
    use pm_matrixSubset, only: lowDia, lowDia_type
    use pm_matrixSubset, only: uppDia, uppDia_type
    use pm_matrixPack, only: rdpack, rdpack_type
    use pm_matrixPack, only: rfpack, rfpack_type
    use pm_matrixPack, only: lfpack, lfpack_type
    use pm_array, only: nothing, nothing_type
    use pm_matrixTrans, only: transHerm, transHerm_type
    use pm_matrixTrans, only: trans_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_matrixChol"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  \legacy
    !>  Return the lower-triangle of the Cholesky factorization \f$L\f$ of the symmetric positive-definite
    !>  real-valued matrix represented by the upper-triangle of the input argument \f$\ms{mat} = L.L^T\f$.<br>
    !>
    !>  \details
    !>  On input, the upper triangle and diagonal of `mat` must be specified, which remains intact on output.<br>
    !>
    !>  \param[inout]   mat     :   The input array of shape `(ndim, ndim)` of type `real` of kind \RKALL.<br>
    !>                              <ul>
    !>                                  <li>    On input, the upper triangle and diagonals of `mat` contain the square symmetric
    !>                                          positive-definite (covariance) matrix whose Cholesky factorization is to be computed.<br>
    !>                                  <li>    On output, the lower triangle of `mat` contains the lower triangle of the Cholesky
    !>                                          factorization of the matrix, while the upper triangle and diagonals remain intact.<br>
    !>                              </ul>
    !>  \param[out]     dia     :   The output vector of the same type and kind as `mat` containing
    !>                              the diagonal elements of the lower triangle of the Cholesky factorization.<br>
    !>                              If the Cholesky factorization fails, `dia(1) = -idim` will be set, where `idim` is the column index causing the singularity.<br>
    !>  \param[in]      ndim    :   The input integer of default kind \IK representing the rank of the input square matrix `(ndim,ndim)`.<br>
    !>
    !>  \interface{setChoLow}
    !>  \code{.F90}
    !>
    !>      use pm_matrixChol, only: setChoLow
    !>
    !>      call setChoLow(mat, dia, ndim) ! Explicit-shape dummy argument unsafe interface (for benchmarking purposes)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  This generic interface with explicit-shape dummy arguments is not recommended for general usage.<br>
    !>  It exists merely for benchmarking purposes against the alternative modern interface [setMatChol](@ref pm_matrixChol::setMatChol).<br>
    !>  However, if used, and the Cholesky factorization fails, `dia(1)` will be set to `-idim` where `idim` is the index of the column causing the singularity.<br>
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  If `ndim = 1`, `mat` will not be touched, only `sqrt(mat)` will be written to the output `dia`.<br>
    !>
    !>  \see
    !>  [getMatChol](@ref pm_matrixChol::getMatChol)<br>
    !>  [setMatChol](@ref pm_matrixChol::setMatChol)<br>
    !>
    !>  \example{setChoLow}
    !>  \include{lineno} example/pm_matrixChol/setChoLow/main.F90
    !>  \compilef{setChoLow}
    !>  \output{setChoLow}
    !>  \include{lineno} example/pm_matrixChol/setChoLow/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixChol](@ref test_pm_matrixChol)
    !>
    !>  \final{setChoLow}
    !>
    !>  \author
    !>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setChoLow

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setChoLow_RK5(mat, dia, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoLow_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK), intent(in)    :: ndim
        real(RKG)  , intent(inout) :: mat(ndim,ndim)
        real(RKG)  , intent(out)   :: dia(ndim)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setChoLow_RK4(mat, dia, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoLow_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK), intent(in)    :: ndim
        real(RKG)  , intent(inout) :: mat(ndim,ndim)
        real(RKG)  , intent(out)   :: dia(ndim)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setChoLow_RK3(mat, dia, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoLow_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK), intent(in)    :: ndim
        real(RKG)  , intent(inout) :: mat(ndim,ndim)
        real(RKG)  , intent(out)   :: dia(ndim)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setChoLow_RK2(mat, dia, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoLow_RK2
#endif
        use pm_kind, only: IK, RKG => RK2
        integer(IK), intent(in)    :: ndim
        real(RKG)  , intent(inout) :: mat(ndim,ndim)
        real(RKG)  , intent(out)   :: dia(ndim)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setChoLow_RK1(mat, dia, ndim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setChoLow_RK1
#endif
        use pm_kind, only: IK, RKG => RK1
        integer(IK), intent(in)    :: ndim
        real(RKG)  , intent(inout) :: mat(ndim,ndim)
        real(RKG)  , intent(out)   :: dia(ndim)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface setChoLow

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the upper or the lower Cholesky factorization of the input symmetric positive-definite
    !>  matrix represented by the upper-triangle of the input argument \f$\ms{mat} = L.L^T\f$.<br>
    !>
    !>  \details
    !>  This generic interface is merely a wrapper for the lower-level, potentially faster,
    !>  more comprehensive generic interface [setMatChol](@ref pm_matrixChol::setMatChol).<br>
    !>
    !>  \param[in]      mat         :   The input square positive-definite matrix of shape `(ndim, ndim)` of
    !>                                  <ol>
    !>                                      <li>     type `complex` of kind \CKALL,
    !>                                      <li>     type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the specified `subset` of the positive definite matrix whose Cholesky factorization is to be computed.<br>
    !>  \param[in]      subset      :   The input scalar that can be either,
    !>                                  <ol>
    !>                                      <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) implying that
    !>                                              the upper-diagonal triangular block of `mat` should be used for building the Cholesky factor
    !>                                              while the lower triangular part of `mat` is not referenced within the algorithm.
    !>                                      <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) implying that
    !>                                              the lower-diagonal triangular block of `mat` should be used for building the Cholesky factor
    !>                                              while the upper triangular part of `mat` is not referenced within the algorithm.
    !>                                  </ol>
    !>                                  This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.<br>
    !>  \param[in]      operation   :   The input scalar that can be either,
    !>                                  <ol>
    !>                                      <li>    the constant [nothing](@ref pm_array::nothing) implying that
    !>                                              the same subset of the output `chol` as that of `mat` specified by the input `subset` must be populated with the Cholesky factorization.<br>
    !>                                      <li>    the constant [transHerm](@ref pm_matrixTrans::transHerm) implying that
    !>                                              the Hermitian transpose of the computed Cholesky factorization must be written to the output `chol`.<br>
    !>                                              For example, specifying `subset = uppDia, operation = transHerm` will lead to using
    !>                                              the upper-diagonal subset of the input `mat` to compute the lower-diagonal triangle of `chol`.<br>
    !>                                  </ol>
    !>                                  (**optional**, default = [nothing](@ref pm_array::nothing))
    !>
    !>  \return
    !>  `chol`                      :   The output matrix of the same type, kind, and shape as the input `mat`
    !>                                  containing the requested subset of the Cholesky factorization of `mat`.<br>
    !>                                  <ol>
    !>                                      <li>    If `subset = uppDia, operation = nothing`, then the upper-diagonal Cholesky factorization of `mat` is computed.<br>
    !>                                      <li>    If `subset = lowDia, operation = nothing`, then the lower-diagonal Cholesky factorization of `mat` is computed.<br>
    !>                                      <li>    If `subset = uppDia, operation = transHerm`, then the lower-diagonal Cholesky factorization of `mat` is computed.<br>
    !>                                      <li>    If `subset = lowDia, operation = transHerm`, then the upper-diagonal Cholesky factorization of `mat` is computed.<br>
    !>                                  </ol>
    !>                                  By definition, the opposite-triangular elements of `chol` that do not contain the Cholesky factorization will be set to zero.<br>
    !>                                  The program will call `error stop` is the Cholesky factorization fails.
    !>
    !>  \interface{getMatChol}
    !>  \code{.F90}
    !>
    !>      use pm_matrixChol, only: getMatChol
    !>
    !>      chol(1:ndim, 1:ndim) = getMatChol(mat(1:ndim, 1:ndim), subset, operation = operation)
    !>
    !>  \endcode
    !>
    !>  \note
    !>  This generic interface is a convenience wrapper around the more efficient
    !>  generic subroutine interface [setMatChol](@ref pm_matrixChol::setMatChol).<br>
    !>
    !>  \pure
    !>  The procedures of this generic interface are `impure` when the output argument `failed` is present.
    !>
    !>  \see
    !>  [pm_matrixDet](@ref pm_matrixDet)<br>
    !>  [pm_matrixInv](@ref pm_matrixInv)<br>
    !>  [pm_matrixLUP](@ref pm_matrixLUP)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [iteration](@ref pm_control::iteration)<br>
    !>  [recursion](@ref pm_control::recursion)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>
    !>  \example{getMatChol}
    !>  \include{lineno} example/pm_matrixChol/getMatChol/main.F90
    !>  \compilef{getMatChol}
    !>  \output{getMatChol}
    !>  \include{lineno} example/pm_matrixChol/getMatChol/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixChol](@ref test_pm_matrixChol)
    !>
    !>  \final{getMatChol}
    !>
    !>  \author
    !>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! subset UXD

    interface getMatChol

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatChol_UXD_CK5(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_UXD_CK5
#endif
        use pm_kind, only: CKG => CK5
        class(*)                , intent(in)    , optional      :: operation
        type(uppDia_type)       , intent(in)                    :: subset
        complex(CKG)            , intent(in)                    :: mat(:,:)
        complex(CKG)                                            :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatChol_UXD_CK4(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_UXD_CK4
#endif
        use pm_kind, only: CKG => CK4
        class(*)                , intent(in)    , optional      :: operation
        type(uppDia_type)       , intent(in)                    :: subset
        complex(CKG)            , intent(in)                    :: mat(:,:)
        complex(CKG)                                            :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatChol_UXD_CK3(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_UXD_CK3
#endif
        use pm_kind, only: CKG => CK3
        class(*)                , intent(in)    , optional      :: operation
        type(uppDia_type)       , intent(in)                    :: subset
        complex(CKG)            , intent(in)                    :: mat(:,:)
        complex(CKG)                                            :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatChol_UXD_CK2(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_UXD_CK2
#endif
        use pm_kind, only: CKG => CK2
        class(*)                , intent(in)    , optional      :: operation
        type(uppDia_type)       , intent(in)                    :: subset
        complex(CKG)            , intent(in)                    :: mat(:,:)
        complex(CKG)                                            :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatChol_UXD_CK1(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_UXD_CK1
#endif
        use pm_kind, only: CKG => CK1
        class(*)                , intent(in)    , optional      :: operation
        type(uppDia_type)       , intent(in)                    :: subset
        complex(CKG)            , intent(in)                    :: mat(:,:)
        complex(CKG)                                            :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatChol_UXD_RK5(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_UXD_RK5
#endif
        use pm_kind, only: RKG => RK5
        class(*)                , intent(in)    , optional      :: operation
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)                    :: mat(:,:)
        real(RKG)                                               :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatChol_UXD_RK4(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_UXD_RK4
#endif
        use pm_kind, only: RKG => RK4
        class(*)                , intent(in)    , optional      :: operation
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)                    :: mat(:,:)
        real(RKG)                                               :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatChol_UXD_RK3(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_UXD_RK3
#endif
        use pm_kind, only: RKG => RK3
        class(*)                , intent(in)    , optional      :: operation
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)                    :: mat(:,:)
        real(RKG)                                               :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatChol_UXD_RK2(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_UXD_RK2
#endif
        use pm_kind, only: RKG => RK2
        class(*)                , intent(in)    , optional      :: operation
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)                    :: mat(:,:)
        real(RKG)                                               :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatChol_UXD_RK1(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_UXD_RK1
#endif
        use pm_kind, only: RKG => RK1
        class(*)                , intent(in)    , optional      :: operation
        type(uppDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)                    :: mat(:,:)
        real(RKG)                                               :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! subset XLD

    interface getMatChol

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatChol_XLD_CK5(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_XLD_CK5
#endif
        use pm_kind, only: CKG => CK5
        class(*)                , intent(in)    , optional      :: operation
        type(lowDia_type)       , intent(in)                    :: subset
        complex(CKG)            , intent(in)                    :: mat(:,:)
        complex(CKG)                                            :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatChol_XLD_CK4(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_XLD_CK4
#endif
        use pm_kind, only: CKG => CK4
        class(*)                , intent(in)    , optional      :: operation
        type(lowDia_type)       , intent(in)                    :: subset
        complex(CKG)            , intent(in)                    :: mat(:,:)
        complex(CKG)                                            :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatChol_XLD_CK3(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_XLD_CK3
#endif
        use pm_kind, only: CKG => CK3
        class(*)                , intent(in)    , optional      :: operation
        type(lowDia_type)       , intent(in)                    :: subset
        complex(CKG)            , intent(in)                    :: mat(:,:)
        complex(CKG)                                            :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatChol_XLD_CK2(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_XLD_CK2
#endif
        use pm_kind, only: CKG => CK2
        class(*)                , intent(in)    , optional      :: operation
        type(lowDia_type)       , intent(in)                    :: subset
        complex(CKG)            , intent(in)                    :: mat(:,:)
        complex(CKG)                                            :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatChol_XLD_CK1(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_XLD_CK1
#endif
        use pm_kind, only: CKG => CK1
        class(*)                , intent(in)    , optional      :: operation
        type(lowDia_type)       , intent(in)                    :: subset
        complex(CKG)            , intent(in)                    :: mat(:,:)
        complex(CKG)                                            :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatChol_XLD_RK5(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_XLD_RK5
#endif
        use pm_kind, only: RKG => RK5
        class(*)                , intent(in)    , optional      :: operation
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)                    :: mat(:,:)
        real(RKG)                                               :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatChol_XLD_RK4(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_XLD_RK4
#endif
        use pm_kind, only: RKG => RK4
        class(*)                , intent(in)    , optional      :: operation
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)                    :: mat(:,:)
        real(RKG)                                               :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatChol_XLD_RK3(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_XLD_RK3
#endif
        use pm_kind, only: RKG => RK3
        class(*)                , intent(in)    , optional      :: operation
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)                    :: mat(:,:)
        real(RKG)                                               :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatChol_XLD_RK2(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_XLD_RK2
#endif
        use pm_kind, only: RKG => RK2
        class(*)                , intent(in)    , optional      :: operation
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)                    :: mat(:,:)
        real(RKG)                                               :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatChol_XLD_RK1(mat, subset, operation) result(chol)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatChol_XLD_RK1
#endif
        use pm_kind, only: RKG => RK1
        class(*)                , intent(in)    , optional      :: operation
        type(lowDia_type)       , intent(in)                    :: subset
        real(RKG)               , intent(in)                    :: mat(:,:)
        real(RKG)                                               :: chol(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Compute and return the lower/upper-triangle of the Cholesky factorization \f$L\f$ of the input Symmetric/Hermitian positive-definite triangular matrix.<br>
    !>
    !>  \param[inout]   mat         :   The input or input/output array of rank `2` of square shape `(1 : ndim, 1 : ndim)` of,
    !>                                  <ol>
    !>                                      <li>     type `complex` of kind \CKALL, or
    !>                                      <li>     type `real` of kind \RKALL.
    !>                                  </ol>
    !>                                  The intent of the argument `mat` depends on the presence or settings of the rest of the arguments of the generic interface.<br>
    !>                                  <ol>
    !>                                      <li>    On input, the specified upper/lower-diagonal `subset` of `mat` must contain the square symmetric positive-definite matrix whose Cholesky factorization is to be computed.<br>
    !>                                      <li>    On output,
    !>                                              <ol>
    !>                                                  <li>    If the optional output argument `chol` is missing, then `mat` has `intent(inout)` and the same specified triangular `subset` of `mat`
    !>                                                          will contain the Cholesky factorization of the matrix, while the other triangle of `mat` remains intact.<br>
    !>                                                          In this case, the blocking algorithms of BLAS/LAPACK are used for the factorization.<br>
    !>                                                  <li>    If the optional output argument `chol` is present, then the input argument `mat` has `intent(in)` and remains intact on return.<br>
    !>                                                          The specified triangular `subset` of `mat` will be used to construct the same (or the opposite) triangular subset of the Cholesky factorization of `mat`.<br>
    !>                                                          The exact storage format of the output `chol` is determined by the input arguments `subset` and `operation`.<br>
    !>                                                          In all cases, the other side of `chol` remains intact upon return.<br>
    !>                                                          For example,
    !>                                                          <ol>
    !>                                                              <li>    setting `subset = uppDia, operation = nothing` will lead to reading the matrix from the upper-diagonal triangular subset of `mat`
    !>                                                                      and writing the upper-triangular Cholesky factorization to the upper-diagonal triangular subset of `chol` while its upper triangle remains intact.<br>
    !>                                                              <li>    setting `subset = uppDia, operation = transHerm` will lead to reading the matrix from the upper-diagonal triangular subset of `mat`
    !>                                                                      and writing the lower-triangular Cholesky factorization to the lower-diagonal triangular subset of `chol` while its upper triangle remains intact.<br>
    !>                                                              <li>    setting `subset = lowDia, operation = nothing` will lead to reading the matrix from the lower-diagonal triangular subset of `mat`
    !>                                                                      and writing the corresponding Cholesky factorization to the upper-diagonal triangular subset of `chol` while its upper triangle remains intact.<br>
    !>                                                              <li>    setting `subset = lowDia, operation = transHerm` will lead to reading the matrix from the lower-diagonal triangular subset of `mat`
    !>                                                                      and writing the lower-triangular Cholesky factorization to the lower-diagonal triangular subset of `chol` while its upper triangle remains intact.<br>
    !>                                                          </ol>
    !>                                                          The input `mat` and `chol` can be two opposite subsets of the same matrix.<br>
    !>                                                          This format is particularly useful for efficient packing of `mat` with its Cholesky factorization with losing the diagonals of either matrix.<br>
    !>                                                          The default algorithms with this format are the non-blocking basic (Cholesky–Crout / Cholesky–Banachiewicz) algorithms.<br>
    !>                                                          These algorithms are very fast for low rank input `mat` but slower than the blocked LAPACK algorithms at higher matrix ranks.<br>
    !>                                              </ol>
    !>                                  </ol>
    !>  \param[in]      subset      :   The input scalar that can be either,
    !>                                  <ol>
    !>                                      <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) implying that
    !>                                              the upper-diagonal triangular block of `mat` should be used for building
    !>                                              the Cholesky factor while the lower triangular part of `mat` remains intact on output.
    !>                                      <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) implying that
    !>                                              the lower-diagonal triangular block of `mat` should be used for building
    !>                                              the Cholesky factor while the upper triangular part of `mat` remains intact on output.
    !>                                  </ol>
    !>                                  This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.
    !>  \param[out]     info        :   The output scalar `integer` of default kind \IK that is `0` <b>if and only if</b> the Cholesky factorization succeeds.<br>
    !>                                  Otherwise, it is set to the order of the leading minor of the specified input subset of `mat` that is not positive definite,
    !>                                  indicating the occurrence of an error and that the factorization could not be completed.<br>
    !>                                  **A non-zero value implies failure of the Cholesky factorization.**<br>
    !>  \param[inout]   chol        :   The input/output matrix of the same type, kind, and shape as the input `mat` containing
    !>                                  the Cholesky factorization in its triangular subset dictated by the input arguments `subset` and `operation`.<br>
    !>                                  See the documentation of the input argument `mat` above for more information.<br>
    !>                                  (**optional**. It must be present **if and only if** the input argument `operation` is present and `control` is missing.)
    !>  \param[in]      operation   :   The input scalar that can be either,
    !>                                  <ol>
    !>                                      <li>    the constant [nothing](@ref pm_array::nothing) implying that
    !>                                              the same subset of the output `chol` as the input `subset` must be populated with the Cholesky factorization.<br>
    !>                                      <li>    the constant [transHerm](@ref pm_matrixTrans::transHerm) implying that
    !>                                              the complex transpose of the computed Cholesky factorization must be written to the output `chol`.<br>
    !>                                  </ol>
    !>                                  (**optional**, default = [nothing](@ref pm_array::nothing). It must be present **if and only if** the output argument `chol` is also present.)
    !>  \param[in]      control     :   The input scalar that can be,<br>
    !>                                  <ol>
    !>                                      <li>    the constant [recursion](@ref pm_control::recursion)
    !>                                              implying that the `recursive` version of the algorithm must be used.
    !>                                  </ol>
    !>                                  (**optional**. It can be present **only if** the input arguments `chol` and `operation` are both missing.)
    !>  \param[in]      ndim        :   The input scalar `integer` of default kind \IK representing the order of the square subset `(roff : roff + ndim, coff : coff + ndim)` of the input matrix of arbitrary shape `(:,:)`.<br>
    !>                                  (**optional**, default = `size(mat, 1)`.)
    !>  \param[in]      roff        :   The input non-negative scalar of type `integer` of default kind \IK containing the <b>off</b>set from the first <b>r</b>ow of the input matrix `mat`,
    !>                                  such that `mat(1 + roff, 1 + roff)` marks the top-left corner of the block of `mat` used in the algorithm.<br>
    !>                                  (**optional**, default = `0`.)
    !>  \param[in]      coff        :   The input non-negative scalar of type `integer` of default kind \IK containing the <b>off</b>set from the first <b>c</b>olumn of the input matrix `mat`,
    !>                                  such that `mat(1 + roff, 1 + roff)` marks the top-left corner of the block of `mat` used in the algorithm.<br>
    !>                                  (**optional**, default = `0`.)
    !>  \param[in]      roffc       :   The input non-negative scalar of type `integer` of default kind \IK containing the <b>off</b>set from the first <b>r</b>ow of the output matrix `chol`,
    !>                                  such that `chol(1 + roffc, 1 + roffc)` marks the top-left corner of the block of `chol` used in the algorithm.<br>
    !>                                  (**optional**, default = `0`.)
    !>  \param[in]      coffc       :   The input non-negative scalar of type `integer` of default kind \IK containing the <b>off</b>set from the first <b>c</b>olumn of the output matrix `chol`,
    !>                                  such that `chol(1 + roffc, 1 + roffc)` marks the top-left corner of the block of `chol` used in the algorithm.<br>
    !>                                  (**optional**, default = `0`.)
    !>  \param[in]      bdim        :   The input scalar `integer` of default kind \IK representing the optimal block size to be used in the **blocked algorithm**.<br>
    !>                                  For any input `ndim <= bdim`, the procedure uses the unblocked recursive algorithm, otherwise, the blocked algorithm.<br>
    !>                                  (**optional**. default = `64`. It can be present **only if** the input argument `control` is missing.)
    !>
    !>  \interface{setMatChol}
    !>  \code{.F90}
    !>
    !>      use pm_matrixChol, only: setMatChol
    !>
    !>      ! simplified interface.
    !>
    !>      call setMatChol(mat, subset, info)
    !>      call setMatChol(mat, subset, info, control) ! control = recursion
    !>      call setMatChol(mat, subset, info, control, bdim = bdim) ! control = iteration
    !>      call setMatChol(mat, subset, info, chol, operation)
    !>
    !>      ! contiguous interface.
    !>
    !>      call setMatChol(mat, subset, info, ndim, roff, coff)
    !>      call setMatChol(mat, subset, info, control, ndim, roff, coff) ! control = recursion
    !>      call setMatChol(mat, subset, info, control, ndim, roff, coff, bdim = bdim) ! control = iteration
    !>      call setMatChol(mat, subset, info, chol, operation, ndim, roff, coff, roffc, coffc)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= ndim` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= [roff, coff])` must hold for the corresponding input arguments.<br>
    !>  The condition `all([roff, coff] + ndim <= shape(mat))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(chol) == shape(mat))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [pm_matrixDet](@ref pm_matrixDet)<br>
    !>  [pm_matrixInv](@ref pm_matrixInv)<br>
    !>  [pm_matrixLUP](@ref pm_matrixLUP)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [iteration](@ref pm_control::iteration)<br>
    !>  [recursion](@ref pm_control::recursion)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>
    !>  \lapack{3.11}
    !>  `SPOTRF2`, `DPOTRF2`, `CPOTRF2`, and `ZPOTRF2`.<br>
    !>
    !>  \example{setMatChol}
    !>  \include{lineno} example/pm_matrixChol/setMatChol/main.F90
    !>  \compilef{setMatChol}
    !>  \output{setMatChol}
    !>  \include{lineno} example/pm_matrixChol/setMatChol/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixChol](@ref test_pm_matrixChol)<br>
    !>
    !>  \naming
    !>  \code{.F90}
    !>      setMC_IMP_ABI_UXD_ONO_CK5
    !>            ||| ||| ||| ||| |||
    !>            ||| ||| ||| ||| The output type and type kind parameter.
    !>            ||| ||| ||| The operation on the output Cholesky subset: ONO (NO operation, nothing) / OTH (Hermitian Transpose).
    !>            ||| ||| The subset of the input mat to be used: UXD (upper-diagonal), XLD (lower-diagonal).
    !>            ||| ||The algorithm control: I/R/X = Iteration/Recursion/unknown (default)
    !>            ||| |The algorithm access pattern: B/N/X = blocking / non-blocking / unknown (default)
    !>            ||| `A` stands for algorithm specified by the following letters.
    !>            The explicitness of the interface: IMP => implicit, EXP => explicit.
    !>  \endcode
    !>
    !>  \final{setMatChol}
    !>
    !>  \author
    !>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! implicit default

    interface setMatChol

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMC_IMP_AXX_UXD_ONO_CK5(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_UXD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMC_IMP_AXX_UXD_ONO_CK4(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_UXD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMC_IMP_AXX_UXD_ONO_CK3(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_UXD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMC_IMP_AXX_UXD_ONO_CK2(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_UXD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMC_IMP_AXX_UXD_ONO_CK1(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_UXD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMC_IMP_AXX_UXD_ONO_RK5(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_UXD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMC_IMP_AXX_UXD_ONO_RK4(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_UXD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMC_IMP_AXX_UXD_ONO_RK3(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_UXD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMC_IMP_AXX_UXD_ONO_RK2(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_UXD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMC_IMP_AXX_UXD_ONO_RK1(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_UXD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMC_IMP_AXX_XLD_ONO_CK5(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_XLD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMC_IMP_AXX_XLD_ONO_CK4(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_XLD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMC_IMP_AXX_XLD_ONO_CK3(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_XLD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMC_IMP_AXX_XLD_ONO_CK2(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_XLD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMC_IMP_AXX_XLD_ONO_CK1(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_XLD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMC_IMP_AXX_XLD_ONO_RK5(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_XLD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMC_IMP_AXX_XLD_ONO_RK4(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_XLD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMC_IMP_AXX_XLD_ONO_RK3(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_XLD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMC_IMP_AXX_XLD_ONO_RK2(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_XLD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMC_IMP_AXX_XLD_ONO_RK1(mat, subset, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_AXX_XLD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit, unblocked, nothing

    interface setMatChol

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_ONO_CK5(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_ONO_CK4(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_ONO_CK3(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_ONO_CK2(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_ONO_CK1(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_ONO_RK5(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_ONO_RK4(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_ONO_RK3(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_ONO_RK2(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_ONO_RK1(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_ONO_CK5(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_ONO_CK4(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_ONO_CK3(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_ONO_CK2(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_ONO_CK1(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_ONO_RK5(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_ONO_RK4(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_ONO_RK3(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_ONO_RK2(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_ONO_RK1(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit, unblocked, transHerm

    interface setMatChol

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_OTH_CK5(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_OTH_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_OTH_CK4(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_OTH_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_OTH_CK3(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_OTH_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_OTH_CK2(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_OTH_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_OTH_CK1(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_OTH_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_OTH_RK5(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_OTH_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_OTH_RK4(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_OTH_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_OTH_RK3(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_OTH_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_OTH_RK2(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_OTH_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMC_IMP_ANI_UXD_OTH_RK1(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_UXD_OTH_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_OTH_CK5(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_OTH_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_OTH_CK4(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_OTH_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_OTH_CK3(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_OTH_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_OTH_CK2(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_OTH_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_OTH_CK1(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_OTH_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_OTH_RK5(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_OTH_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_OTH_RK4(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_OTH_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_OTH_RK3(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_OTH_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_OTH_RK2(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_OTH_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMC_IMP_ANI_XLD_OTH_RK1(mat, subset, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ANI_XLD_OTH_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit, blocking, iteration

    interface setMatChol

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_UXD_ONO_CK5(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_UXD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_UXD_ONO_CK4(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_UXD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_UXD_ONO_CK3(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_UXD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_UXD_ONO_CK2(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_UXD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_UXD_ONO_CK1(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_UXD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_UXD_ONO_RK5(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_UXD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_UXD_ONO_RK4(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_UXD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_UXD_ONO_RK3(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_UXD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_UXD_ONO_RK2(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_UXD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_UXD_ONO_RK1(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_UXD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_XLD_ONO_CK5(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_XLD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_XLD_ONO_CK4(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_XLD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_XLD_ONO_CK3(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_XLD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_XLD_ONO_CK2(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_XLD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_XLD_ONO_CK1(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_XLD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_XLD_ONO_RK5(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_XLD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_XLD_ONO_RK4(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_XLD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_XLD_ONO_RK3(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_XLD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_XLD_ONO_RK2(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_XLD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setMC_IMP_ABI_XLD_ONO_RK1(mat, subset, info, control, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABI_XLD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit, blocking, recursion

    interface setMatChol

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_UXD_ONO_CK5(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_UXD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_UXD_ONO_CK4(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_UXD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_UXD_ONO_CK3(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_UXD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_UXD_ONO_CK2(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_UXD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_UXD_ONO_CK1(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_UXD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_UXD_ONO_RK5(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_UXD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_UXD_ONO_RK4(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_UXD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_UXD_ONO_RK3(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_UXD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_UXD_ONO_RK2(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_UXD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_UXD_ONO_RK1(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_UXD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_XLD_ONO_CK5(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_XLD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_XLD_ONO_CK4(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_XLD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_XLD_ONO_CK3(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_XLD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_XLD_ONO_CK2(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_XLD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_XLD_ONO_CK1(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_XLD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_XLD_ONO_RK5(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_XLD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_XLD_ONO_RK4(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_XLD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_XLD_ONO_RK3(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_XLD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_XLD_ONO_RK2(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_XLD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setMC_IMP_ABR_XLD_ONO_RK1(mat, subset, info, control)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_IMP_ABR_XLD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! explicit, default

    interface setMatChol

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMC_EXP_AXX_UXD_ONO_CK5(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_UXD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMC_EXP_AXX_UXD_ONO_CK4(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_UXD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMC_EXP_AXX_UXD_ONO_CK3(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_UXD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMC_EXP_AXX_UXD_ONO_CK2(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_UXD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMC_EXP_AXX_UXD_ONO_CK1(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_UXD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMC_EXP_AXX_UXD_ONO_RK5(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_UXD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMC_EXP_AXX_UXD_ONO_RK4(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_UXD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMC_EXP_AXX_UXD_ONO_RK3(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_UXD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMC_EXP_AXX_UXD_ONO_RK2(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_UXD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMC_EXP_AXX_UXD_ONO_RK1(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_UXD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMC_EXP_AXX_XLD_ONO_CK5(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_XLD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMC_EXP_AXX_XLD_ONO_CK4(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_XLD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMC_EXP_AXX_XLD_ONO_CK3(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_XLD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMC_EXP_AXX_XLD_ONO_CK2(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_XLD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMC_EXP_AXX_XLD_ONO_CK1(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_XLD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMC_EXP_AXX_XLD_ONO_RK5(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_XLD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMC_EXP_AXX_XLD_ONO_RK4(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_XLD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMC_EXP_AXX_XLD_ONO_RK3(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_XLD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMC_EXP_AXX_XLD_ONO_RK2(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_XLD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMC_EXP_AXX_XLD_ONO_RK1(mat, subset, info, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_AXX_XLD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! explicit, unblocked, nothing

    interface setMatChol

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_ONO_CK5(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_ONO_CK4(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_ONO_CK3(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_ONO_CK2(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_ONO_CK1(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_ONO_RK5(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_ONO_RK4(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_ONO_RK3(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_ONO_RK2(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_ONO_RK1(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_ONO_CK5(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_ONO_CK4(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_ONO_CK3(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_ONO_CK2(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_ONO_CK1(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_ONO_RK5(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_ONO_RK4(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_ONO_RK3(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_ONO_RK2(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_ONO_RK1(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(nothing_type)  , intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! explicit, unblocked, transHerm

    interface setMatChol

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_OTH_CK5(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_OTH_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_OTH_CK4(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_OTH_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_OTH_CK3(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_OTH_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_OTH_CK2(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_OTH_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_OTH_CK1(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_OTH_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_OTH_RK5(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_OTH_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_OTH_RK4(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_OTH_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_OTH_RK3(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_OTH_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_OTH_RK2(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_OTH_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMC_EXP_ANI_UXD_OTH_RK1(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_UXD_OTH_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_OTH_CK5(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_OTH_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_OTH_CK4(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_OTH_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_OTH_CK3(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_OTH_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_OTH_CK2(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_OTH_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_OTH_CK1(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_OTH_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        complex(CKG)        , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKG)        , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_OTH_RK5(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_OTH_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_OTH_RK4(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_OTH_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_OTH_RK3(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_OTH_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_OTH_RK2(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_OTH_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMC_EXP_ANI_XLD_OTH_RK1(mat, subset, info, chol, operation, ndim, roff, coff, rofc, cofc)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ANI_XLD_OTH_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff, rofc, cofc
        real(RKG)           , intent(in)    , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKG)           , intent(inout) , contiguous    :: chol(1 - rofc :, 1 - cofc :)
        type(transHerm_type), intent(in)                    :: operation
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! explicit, blocking, iteration

    interface setMatChol

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_UXD_ONO_CK5(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_UXD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_UXD_ONO_CK4(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_UXD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_UXD_ONO_CK3(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_UXD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_UXD_ONO_CK2(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_UXD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_UXD_ONO_CK1(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_UXD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_UXD_ONO_RK5(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_UXD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_UXD_ONO_RK4(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_UXD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_UXD_ONO_RK3(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_UXD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_UXD_ONO_RK2(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_UXD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_UXD_ONO_RK1(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_UXD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_XLD_ONO_CK5(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_XLD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_XLD_ONO_CK4(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_XLD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_XLD_ONO_CK3(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_XLD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_XLD_ONO_CK2(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_XLD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_XLD_ONO_CK1(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_XLD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_XLD_ONO_RK5(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_XLD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_XLD_ONO_RK4(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_XLD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_XLD_ONO_RK3(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_XLD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_XLD_ONO_RK2(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_XLD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setMC_EXP_ABI_XLD_ONO_RK1(mat, subset, info, control, ndim, roff, coff, bdim)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABI_XLD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(iteration_type), intent(in)                    :: control
        integer(IK)         , intent(in)    , optional      :: bdim
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! explicit, blocking, recursion

    interface setMatChol

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_UXD_ONO_CK5(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_UXD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_UXD_ONO_CK4(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_UXD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_UXD_ONO_CK3(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_UXD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_UXD_ONO_CK2(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_UXD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_UXD_ONO_CK1(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_UXD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_UXD_ONO_RK5(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_UXD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_UXD_ONO_RK4(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_UXD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_UXD_ONO_RK3(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_UXD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_UXD_ONO_RK2(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_UXD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_UXD_ONO_RK1(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_UXD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_XLD_ONO_CK5(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_XLD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_XLD_ONO_CK4(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_XLD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_XLD_ONO_CK3(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_XLD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_XLD_ONO_CK2(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_XLD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_XLD_ONO_CK1(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_XLD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        complex(CKG)        , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_XLD_ONO_RK5(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_XLD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_XLD_ONO_RK4(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_XLD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_XLD_ONO_RK3(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_XLD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_XLD_ONO_RK2(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_XLD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setMC_EXP_ABR_XLD_ONO_RK1(mat, subset, info, control, ndim, roff, coff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMC_EXP_ABR_XLD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        integer(IK)         , intent(out)                   :: info
        integer(IK)         , intent(in)                    :: ndim, roff, coff
        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)   , intent(in)                    :: subset
        type(recursion_type), intent(in)                    :: control
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!contains
!
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!>  \cond excluded
!   >  \brief
!   >  Return the Cholesky factorization of the input positive-definite matrix.
!   >
!   >  \param[in]      ndim            :   The size of the input square matrix - `ndim` by `ndim`.
!   >  \param[inout]   ChoUppCovLow    :   The input square positive-definite (covariance) matrix. Only the lower half and diagonal is needed.
!   >                                       The lower half will be rewritten with the lower triangle of the Cholesky Factorization.
!   >  \param[out]     dia          :   The diagonal elements of the Cholesky factorization.
!   >
!   >  \remark
!   >  If `ndim = 1`, `ChoUppCovLow` will not be touched, only `sqrt(ChoUppCovLow)` will be written to the output `dia`.
!   >
!   >  \warning
!   >  If Cholesky factorization fails, `dia(1)` will be set to `-1` to indicate error on return.
!   >
!   >  \remark
!   >  This routine is particularly useful for computing the Cholesky upper triangle of lower-triangle
!   >  inverse matrices returned by [getInvLowFromChoLow](@ref getInvLowFromChoLow).
!   >
!   >  \details
!   >  Returns in the upper triangle of `ChoUppCovLow`, the Cholesky factorization \f$L^T\f$ of \f$\ms{ChoUppCovLow} = L.L^T\f$.
!   >  On input, the lower triangle of `ChoUppCovLow` must be given, which remains intact on output.<br>
!   >
!   >  \final
!   >
!   >  \author
!   >  Amir Shahmoradi, Friday 1:53 AM, August 27, 2021, Dallas, TX
!    pure subroutine setChoUpp(ndim,ChoUppCovLow,dia)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setChoUpp
!#endif
!        use pm_kind, only: RK, IK
!        implicit none
!        integer(IK), intent(in)    :: ndim
!        real(RKG)  , intent(inout) :: ChoUppCovLow(ndim,ndim) ! Upper triangle + diagonal is input matrix, lower is output.
!        real(RKG)  , intent(out)   :: dia(ndim)
!        real(RKG)                  :: summ
!        integer(IK)                :: i
!        do i=1,ndim
!            summ = ChoUppCovLow(i,i) - dot_product(ChoUppCovLow(1:i-1,i),ChoUppCovLow(1:i-1,i))
!            if (summ <= 0._ONO_RK) then
!                dia(1) = -1._ONO_RK
!                return
!            end if
!            dia(i) = sqrt(summ)
!            ChoUppCovLow(i,i+1:ndim) = ( ChoUppCovLow(i+1:ndim,i) - matmul(ChoUppCovLow(1:i-1,i),ChoUppCovLow(1:i-1,i+1:ndim)) ) / dia(i)
!        end do
!    end subroutine setChoUpp
!>  \endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_matrixChol ! LCOV_EXCL_LINE