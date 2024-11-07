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
!>  This module contains procedures and generic interfaces relevant to the partially LU Pivoted decomposition of matrix operations and linear algebra.
!>
!>  \details
!>  In numerical analysis and linear algebra, lower–upper (LU) decomposition or factorization
!>  factors a matrix as the product of a lower triangular matrix and an upper triangular matrix.<br>
!>  The product sometimes also involves a **permutation matrix**.<br>
!>  LU decomposition can be viewed as the matrix form of **Gaussian elimination**.<br>
!>  Computers usually solve square systems of linear equations using LU decomposition,
!>  and it is also a key step when **inverting a matrix** or **computing the determinant** of a matrix.<br>
!>  The LU decomposition was introduced by the Polish astronomer Tadeusz Banachiewicz in 1938.<br>
!>  The LU factorization is sometimes also referred to as **LR decomposition**, factoring into left and right triangular matrices.<br>
!>
!>  Definition
!>  ==========
!>
!>  Let \f$A\f$ be a square matrix.<br>
!>  An LU factorization refers to the factorization of \f$A\f$, with proper row and/or column orderings or permutations,
!>  into two factors – a lower triangular matrix \f$L\f$ and an upper triangular matrix \f$U\f$:<br>
!>  \f{equation}{
!>      A = LU ~.
!>  \f}
!>  In the lower triangular matrix all elements above the diagonal are zero.<br>
!>  In the upper triangular matrix all the elements below the diagonal are zero.<br>
!>  For example, for a \f$3 \times 3\f$ matrix \f$A\f$, its LU decomposition looks like this:<br>
!>  \f{equation}{
!>      \begin{bmatrix}
!>          a_{11}&a_{12}&a_{13} \\
!>          a_{21}&a_{22}&a_{23} \\
!>          a_{31}&a_{32}&a_{33}
!>      \end{bmatrix} =
!>      \begin{bmatrix}
!>          \ell_{11}&0&0\\
!>          \ell_{21}&\ell_{22}&0\\
!>          \ell_{31}&\ell_{32}&\ell_{33}
!>      \end{bmatrix}
!>      \begin{bmatrix}
!>          u_{11}&u_{12}&u_{13}\\
!>          0&u_{22}&u_{23}\\
!>          0&0&u_{33}
!>      \end{bmatrix} ~.
!>  \f}
!>  Without a proper ordering or permutations in the matrix, the factorization may fail to materialize.<br>
!>  For example, it is easy to verify (by expanding the matrix multiplication) that \f$a_{11}=\ell_{11}u_{11}\f$.<br>
!>  If \f$a_{11}=0\f$, then at least one of \f$\ell_{11}\f$ and \f$u_{11}\f$ has to be zero, which implies that either \f$L\f$ or \f$U\f$ is singular.<br>
!>  This is impossible if \f$A\f$ is nonsingular (invertible).<br>
!>  This is a procedural problem.<br>
!>  It can be removed by simply reordering the rows of \f$A\f$ so that the first element of the permuted matrix is nonzero.<br>
!>  The same problem in subsequent factorization steps can be removed the same way; see the basic procedure below.<br>
!>
!>  LU factorization with partial pivoting
!>  --------------------------------------
!>
!>  It turns out that a proper permutation in rows (or columns) is sufficient for LU factorization.<br>
!>  LU factorization with partial pivoting (LUP) refers often to LU factorization with row permutations only:<br>
!>  \f{equation}{
!>      PA = LU ~,
!>  \f}
!>  where \f$L\f$ and \f$U\f$ are again lower and upper triangular matrices, and \f$P\f$ is a **permutation matrix**, which, when left-multiplied to \f$A\f$, reorders the rows of \f$A\f$.<br>
!>  It turns out that all square matrices can be factorized in this form, and the factorization is numerically stable in practice.<br>
!>  This makes LUP decomposition a useful technique in practice.<br>
!>
!>  LU factorization with full pivoting
!>  -----------------------------------
!>
!>  An LU factorization with full pivoting involves both row and column permutations:<br>
!>  \f{equation}{
!>      PAQ = LU ~,
!>  \f}
!>  where \f$L\f$, \f$U\f$ and \f$P\f$ are defined as before, and \f$Q\f$ is a permutation matrix that reorders the columns of \f$A\f$.<br>
!>
!>  \note
!>  To solve systems of equations using the output LU factorization from the generic interfaces of this module, see [pm_matrixMulTri](@ref pm_matrixMulTri).<br>
!>
!>  \see
!>  [pm_matrixInv](@ref pm_matrixInv)<br>
!>  [pm_matrixChol](@ref pm_matrixChol)<br>
!>
!>  \test
!>  [test_pm_matrixLUP](@ref test_pm_matrixLUP)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_matrixLUP

    use pm_kind, only: SK, IK, LK
    !use pm_control, only: iteration, iteration_type
    !use pm_control, only: recursion, recursion_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_matrixLUP"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the LU-Pivoted decomposition of the input square matrix `mat(ndim,ndim)`.<br>
    !>
    !>  \details
    !>  The LUP factorization takes the form,
    !>  \f{equation}{
    !>      \text{mat} = P * L * U
    !>  \f}
    !>  where,<br>
    !>  <ol>
    !>      <li>    \f$P\f$ is the permutation matrix,
    !>      <li>    \f$L\f$ is lower triangular with unit diagonal elements (that are missing in the output `mat`),
    !>      <li>    \f$U\f$ is upper triangular of the decomposition.
    !>  </ol>
    !>  The abbreviation **LUP** stands for the LU factorization with partial Pivoting.<br>
    !>
    !>  \param[inout]   mat     :   The input `contiguous` square matrix of shape `(ndim, ndim)` of,<br>
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              containing the matrix whose LUP factorization must be returned.<br>
    !>                              On output, `mat` is completely overwritten by its LU Pivoted decomposition.<br>
    !>  \param[out]     rperm   :   The output `contiguous` vector of size `ndim` of type `integer` of default kind \IK containing the pivot indices (row permutations).<br>
    !>                              For `1 <= i <= ndim`, the `i`th row of the matrix is interchanged with row `rperm(i)`.<br>
    !>  \param[out]     info    :   The output scalar `integer` of default kind \IK that is non-zero if
    !>                              a singular matrix is detected, indicating the LUP decomposition failure.<br>
    !>
    !>  \interface{setMatLUP}
    !>  \code{.F90}
    !>
    !>      use pm_matrixLUP, only: setMatLUP
    !>
    !>      call setMatLUP(mat(1 : ndim, 1 : ndim), rperm(1 : ndim), info)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(mat, 1) == size(mat, 2))` must hold for the corresponding input arguments (unless dispatch to LAPACK is enabled).<br>
    !>  The condition `size(rperm) == shape(mat, 1, IK)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  This routine can be used in combination with the generic interfaces of [pm_matrixMulTri](@ref pm_matrixMulTri)
    !>  to iteratively solve systems of linear equations or to [invert matrices](@ref pm_matrixInv).<br>
    !>
    !>  \see
    !>  [getMatInv](@ref pm_matrixInv::getMatInv)<br>
    !>  [setMatInv](@ref pm_matrixInv::setMatInv)<br>
    !>  [getMatChol](@ref pm_matrixChol::getMatChol)<br>
    !>  [setMatChol](@ref pm_matrixChol::setMatChol)<br>
    !>  [setMatLUP](@ref pm_matrixLUP::setMatLUP)<br>
    !>
    !>  \lapack{3.11}
    !>  `SGETRF`, `DGETRF`, `CGETRF`, and `ZGETRF`.<br>
    !>
    !>  \example{setMatLUP}
    !>  \include{lineno} example/pm_matrixLUP/setMatLUP/main.F90
    !>  \compilef{setMatLUP}
    !>  \output{setMatLUP}
    !>  \include{lineno} example/pm_matrixLUP/setMatLUP/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixLUP](@ref test_pm_matrixLUP)
    !>
    !>  \final{setMatLUP}
    !>
    !>  \author
    !>  \AmirShahmoradi, Apr 21, 2017, 1:43 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! implicit interface.

    interface setMatLUP

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatLUP_IMP_SQM_CK5(mat, rperm, info) ! , parity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_SQM_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
       !complex(CKG)        , intent(out)   , optional      :: parity
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatLUP_IMP_SQM_CK4(mat, rperm, info) ! , parity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_SQM_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
       !complex(CKG)        , intent(out)   , optional      :: parity
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatLUP_IMP_SQM_CK3(mat, rperm, info) ! , parity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_SQM_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
       !complex(CKG)        , intent(out)   , optional      :: parity
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatLUP_IMP_SQM_CK2(mat, rperm, info) ! , parity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_SQM_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
       !complex(CKG)        , intent(out)   , optional      :: parity
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatLUP_IMP_SQM_CK1(mat, rperm, info) ! , parity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_SQM_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
       !complex(CKG)        , intent(out)   , optional      :: parity
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatLUP_IMP_SQM_RK5(mat, rperm, info) ! , parity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_SQM_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
       !real(RKG)           , intent(out)   , optional      :: parity
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatLUP_IMP_SQM_RK4(mat, rperm, info) ! , parity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_SQM_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
       !real(RKG)           , intent(out)   , optional      :: parity
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatLUP_IMP_SQM_RK3(mat, rperm, info) ! , parity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_SQM_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
       !real(RKG)           , intent(out)   , optional      :: parity
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatLUP_IMP_SQM_RK2(mat, rperm, info) ! , parity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_SQM_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
       !real(RKG)           , intent(out)   , optional      :: parity
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatLUP_IMP_SQM_RK1(mat, rperm, info) ! , parity)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_SQM_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
       !real(RKG)           , intent(out)   , optional      :: parity
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!#if RK5_ENABLED
!    PURE module subroutine setMatLUP_IMP_ITE_RK5(mat, rperm, info, control)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_ITE_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(iteration_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setMatLUP_IMP_ITE_RK4(mat, rperm, info, control)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_ITE_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(iteration_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setMatLUP_IMP_ITE_RK3(mat, rperm, info, control)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_ITE_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(iteration_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setMatLUP_IMP_ITE_RK2(mat, rperm, info, control)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_ITE_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(iteration_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setMatLUP_IMP_ITE_RK1(mat, rperm, info, control)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_ITE_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(iteration_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE recursive module subroutine setMatLUP_IMP_REC_RK5(mat, rperm, info, control)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_REC_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(recursion_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE recursive module subroutine setMatLUP_IMP_REC_RK4(mat, rperm, info, control)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_REC_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(recursion_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE recursive module subroutine setMatLUP_IMP_REC_RK3(mat, rperm, info, control)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_REC_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(recursion_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE recursive module subroutine setMatLUP_IMP_REC_RK2(mat, rperm, info, control)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_REC_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(recursion_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE recursive module subroutine setMatLUP_IMP_REC_RK1(mat, rperm, info, control)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_IMP_REC_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(recursion_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!    ! explicit interface.
!
!    interface setMatLUP
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setMatLUP_EXP_SQM_RK5(mat, rperm, info, parity, nrow, ncol, roff, coff)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_EXP_SQM_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        integer(IK)         , intent(in)                    :: nrow, ncol, roff, coff
!        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!       !real(RKG)           , intent(out)   , optional      :: parity
!        integer(IK)         , intent(out)                   :: info
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setMatLUP_EXP_SQM_RK4(mat, rperm, info, parity, nrow, ncol, roff, coff)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_EXP_SQM_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        integer(IK)         , intent(in)                    :: nrow, ncol, roff, coff
!        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!       !real(RKG)           , intent(out)   , optional      :: parity
!        integer(IK)         , intent(out)                   :: info
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setMatLUP_EXP_SQM_RK3(mat, rperm, info, parity, nrow, ncol, roff, coff)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_EXP_SQM_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        integer(IK)         , intent(in)                    :: nrow, ncol, roff, coff
!        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!       !real(RKG)           , intent(out)   , optional      :: parity
!        integer(IK)         , intent(out)                   :: info
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setMatLUP_EXP_SQM_RK2(mat, rperm, info, parity, nrow, ncol, roff, coff)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_EXP_SQM_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        integer(IK)         , intent(in)                    :: nrow, ncol, roff, coff
!        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!       !real(RKG)           , intent(out)   , optional      :: parity
!        integer(IK)         , intent(out)                   :: info
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setMatLUP_EXP_SQM_RK1(mat, rperm, info, parity, nrow, ncol, roff, coff)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_EXP_SQM_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        integer(IK)         , intent(in)                    :: nrow, ncol, roff, coff
!        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!       !real(RKG)           , intent(out)   , optional      :: parity
!        integer(IK)         , intent(out)                   :: info
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE module subroutine setMatLUP_EXP_ITE_RK5(mat, rperm, info, control, nrow, ncol, roff, coff)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_EXP_ITE_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        integer(IK)         , intent(in)                    :: nrow, ncol, roff, coff
!        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(iteration_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE module subroutine setMatLUP_EXP_ITE_RK4(mat, rperm, info, control, nrow, ncol, roff, coff)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_EXP_ITE_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        integer(IK)         , intent(in)                    :: nrow, ncol, roff, coff
!        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(iteration_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE module subroutine setMatLUP_EXP_ITE_RK3(mat, rperm, info, control, nrow, ncol, roff, coff)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_EXP_ITE_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        integer(IK)         , intent(in)                    :: nrow, ncol, roff, coff
!        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(iteration_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE module subroutine setMatLUP_EXP_ITE_RK2(mat, rperm, info, control, nrow, ncol, roff, coff)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_EXP_ITE_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        integer(IK)         , intent(in)                    :: nrow, ncol, roff, coff
!        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(iteration_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE module subroutine setMatLUP_EXP_ITE_RK1(mat, rperm, info, control, nrow, ncol, roff, coff)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_EXP_ITE_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        integer(IK)         , intent(in)                    :: nrow, ncol, roff, coff
!        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(iteration_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!#if RK5_ENABLED
!    PURE recursive module subroutine setMatLUP_EXP_REC_RK5(mat, rperm, info, control, nrow, ncol, roff, coff)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_EXP_REC_RK5
!#endif
!        use pm_kind, only: RKG => RK5
!        integer(IK)         , intent(in)                    :: nrow, ncol, roff, coff
!        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(recursion_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!#if RK4_ENABLED
!    PURE recursive module subroutine setMatLUP_EXP_REC_RK4(mat, rperm, info, control, nrow, ncol, roff, coff)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_EXP_REC_RK4
!#endif
!        use pm_kind, only: RKG => RK4
!        integer(IK)         , intent(in)                    :: nrow, ncol, roff, coff
!        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(recursion_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!#if RK3_ENABLED
!    PURE recursive module subroutine setMatLUP_EXP_REC_RK3(mat, rperm, info, control, nrow, ncol, roff, coff)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_EXP_REC_RK3
!#endif
!        use pm_kind, only: RKG => RK3
!        integer(IK)         , intent(in)                    :: nrow, ncol, roff, coff
!        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(recursion_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!#if RK2_ENABLED
!    PURE recursive module subroutine setMatLUP_EXP_REC_RK2(mat, rperm, info, control, nrow, ncol, roff, coff)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_EXP_REC_RK2
!#endif
!        use pm_kind, only: RKG => RK2
!        integer(IK)         , intent(in)                    :: nrow, ncol, roff, coff
!        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(recursion_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!#if RK1_ENABLED
!    PURE recursive module subroutine setMatLUP_EXP_REC_RK1(mat, rperm, info, control, nrow, ncol, roff, coff)
!#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
!        !DEC$ ATTRIBUTES DLLEXPORT :: setMatLUP_EXP_REC_RK1
!#endif
!        use pm_kind, only: RKG => RK1
!        integer(IK)         , intent(in)                    :: nrow, ncol, roff, coff
!        real(RKG)           , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
!        integer(IK)         , intent(out)   , contiguous    :: rperm(:)
!        integer(IK)         , intent(out)                   :: info
!        type(recursion_type), intent(in)                    :: control
!    end subroutine
!#endif
!
!    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    end interface
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_matrixLUP ! LCOV_EXCL_LINE