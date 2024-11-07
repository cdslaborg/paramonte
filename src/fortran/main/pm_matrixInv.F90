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
!>  This module contains abstract and concrete derived types and procedures related to the inversion of square matrices.<br>
!>
!>  Inversion operation
!>  ===================
!>
!>  In linear algebra, an \f$n\f$-by-\f$n\f$ **square** matrix \f$A\f$ is called **invertible** (also **nonsingular**, **nondegenerate**),
!>  if there exists an \f$n\f$-by-\f$n\f$ square matrix \f$B\f$ such that,
!>  \f{equation}{
!>      \mathbf{AB} = \mathbf{BA} = \mathbf{I}_{n}~,
!>  \f}
!>  where \f$I_n\f$ denotes the \f$n\f$-by-\f$n\f$ identity matrix and the multiplication used is ordinary matrix multiplication.<br>
!>  If this is the case, then the matrix \f$B\f$ is uniquely determined by \f$A\f$, and is called the **(multiplicative) inverse** of \f$A\f$, denoted by \f$A^{−1}\f$.<br>
!>  source inversion is the process of finding the matrix \f$B\f$ that satisfies the prior equation for a given invertible matrix \f$A\f$.<br>
!>  Over a field, a square matrix that is **not invertible** is called **singular** or **degenerate**.<br>
!>  A square matrix with entries in a field is **singular if and only if its determinant is zero**.<br>
!>  **Singular matrices are rare** in the sense that if a square matrix entries are randomly selected from any bounded region on the number line or complex plane,
!>  the probability that the matrix is singular is \f$0\f$, that is, it will *almost never* be singular.<br>
!>  **Non-square matrices do not have an inverse**.<br>
!>  However, in some cases such a matrix may have a **left inverse** or **right inverse**.<br>
!>  If \f$A\f$ is \f$m\f$-by-\f$n\f$ and the rank of \f$A\f$ is equal to \f$n\f$ (\f$n \leq m\f$),
!>  then \f$A\f$ has a **left inverse**, an \f$n\f$-by-\f$m\f$ matrix \f$B\f$ such that \f$BA = I_n\f$.<br>
!>  If \f$A\f$ has rank \f$m\f$ (\f$m \leq n\f$), then it has a **right inverse**, an \f$n\f$-by-\f$m\f$ matrix \f$B\f$ such that \f$AB = I_m\f$.<br>
!>
!>  Inverse matrix properties
!>  -------------------------
!>
!>  The following properties hold for an invertible matrix \f$A\f$:<br>
!>  <ol>
!>      <li>    \f$(\mathbf{A}^{-1})^{-1} = \mathbf{A}\f$.<br>
!>      <li>    \f$(k\mathbf{A})^{-1} = k^{-1}\mathbf{A}^{-1}\f$ for nonzero scalar \f$k\f$.<br>
!>      <li>    \f$(\mathbf{Ax})^{+} = \mathbf{x}^{+}\mathbf{A}^{-1}\f$ if \f$A\f$ has orthonormal columns, where \f$+\f$ denotes the Moore–Penrose inverse and \f$x\f$ is a vector.<br>
!>      <li>    \f$(\mathbf{A}^{\mathrm{T}})^{-1} = (\mathbf{A}^{-1})^{\mathrm{T}}\f$.<br>
!>      <li>    For any invertible \f$n\f$-by-\f$n\f$ matrices \f$A\f$ and \f$B\f$, \f$(\mathbf{AB})^{-1} = \mathbf{B}^{-1}\mathbf{A}^{-1}\f$.<br>
!>              More generally, if \f$\mathbf{A}_1, \dots, \mathbf{A}_{k}\f$ are invertible \f$n\f$-by-\f$n\f$ matrices, then
!>              \f$(\mathbf{A}_{1}\mathbf{A}_{2} \cdots \mathbf{A}_{k-1}\mathbf{A}_{k})^{-1} = \mathbf{A}_{k}^{-1}\mathbf{A}_{k-1}^{-1} \cdots \mathbf{A}_{2}^{-1}\mathbf{A}_{1}^{-1}\f$.<br>
!>      <li>    \f$\det\mathbf{A}^{-1} = (\det \mathbf{A})^{-1}\f$.<br>
!>      <li>    The rows of the inverse matrix \f$V\f$ of a matrix \f$U\f$ are orthonormal to the columns of \f$U\f$.<br>
!>  </ol>
!>
!>  Inverse matrix computation
!>  --------------------------
!>
!>  The common approach to computing the inverse matrix stems from its definition.<br>
!>  The inverse matrix \f$A^{-1}\f$ of a square matrix \f$A\f$ is a square matrix such that \f$AA^{-1} = I\f$, where \f$I\f$ is the identity matrix.<br>
!>  Depending on the class of the square matrix \f$A\f$, there are several approaches that can be taken to compute its inverse.<br>
!>  However, all such methods attempt to factorize the matrix first (for example, [Cholesky decomposition](@ref pm_matrixChol)
!>  or [LU decomposition](@ref pm_matrixLUP)) and cast the problem into seeking the solution to a system of linear equations.<br>
!>  For example, for a general square matrix of shape \f$3\times 3\f$, the corresponding system of equations to solve would be:<br>
!>  \f{equation}{
!>      \begin{bmatrix}
!>          A_{11} & A_{12} & A_{13} \\
!>          A_{21} & A_{22} & A_{23} \\
!>          A_{31} & A_{32} & A_{33}
!>      \end{bmatrix}
!>      \begin{bmatrix}
!>          A_{11}^{-1} & A_{12}^{-1} & A_{13}^{-1} \\
!>          A_{21}^{-1} & A_{22}^{-1} & A_{23}^{-1} \\
!>          A_{31}^{-1} & A_{32}^{-1} & A_{33}^{-1}
!>      \end{bmatrix}
!>      =
!>      \begin{bmatrix}
!>          1 & 0 & 0 \\
!>          0 & 1 & 0 \\
!>          0 & 0 & 1
!>      \end{bmatrix}
!>      ~,
!>  \f}
!>  where the second matrix on the left hand side is the inverse of the first square matrix \f$A\f$.<br>
!>  The inverse matrix can be constructed as the collection of solutions to \f$n = 3\f$ systems of equations of the form \f$Ax=b\f$
!>  with different right hand sides \f$b\f$ matrices and \f$x\f$ representing \f$n^{\mathrm{th}}\f$ column of the inverse matrix.<br>
!>  The first system of equations to solve for the above problem would be:<br>
!>  \f{equation}{
!>      \begin{bmatrix}
!>          A_{11} & A_{12} & A_{13} \\
!>          A_{21} & A_{22} & A_{23} \\
!>          A_{31} & A_{32} & A_{33}
!>      \end{bmatrix}
!>      \begin{bmatrix}
!>          A_{11}^{-1} \\
!>          A_{21}^{-1} \\
!>          A_{31}^{-1}
!>      \end{bmatrix}
!>      =
!>      \begin{bmatrix}
!>          1 \\
!>          0 \\
!>          0 \\
!>      \end{bmatrix}
!>      ~.
!>  \f}
!>  The second column of the inverse can be computed by changing \f$b\f$ to \f$[0,1,0]^T\f$, the third column with \f$[0,0,1]^T\f$, and so on.<br>
!>  The task of computing the inverse of the matrix is now reduced to solving a series of systems of linear equations.<br>
!>  This can be done if the matrix \f$A\f$ is factorized into lower and upper triangular matrices,<br>
!>  \f{equation}{
!>      \begin{bmatrix}
!>          A_{11} & A_{12} & A_{13} \\
!>          A_{21} & A_{22} & A_{23} \\
!>          A_{31} & A_{32} & A_{33}
!>      \end{bmatrix}
!>      =
!>      \begin{bmatrix}
!>          l_{11} & 0 & 0 \\
!>          l_{21} & l_{22} & 0 \\
!>          l_{31} & l_{32} & l_{33}
!>      \end{bmatrix}
!>      \begin{bmatrix}
!>          u_{11} & u_{12} & u_{13} \\
!>          0 & u_{22} & u_{23} \\
!>          0 & 0 & u_{33} ~,
!>      \end{bmatrix}
!>  \f}
!>  such that the above system can be written as,<br>
!>  \f{equation}{
!>      \begin{bmatrix}
!>          l_{11} & 0 & 0 \\
!>          l_{21} & l_{22} & 0 \\
!>          l_{31} & l_{32} & l_{33}
!>      \end{bmatrix}
!>      \left(
!>          \begin{bmatrix}
!>              u_{11} & u_{12} & u_{13} \\
!>              0 & u_{22} & u_{23} \\
!>              0 & 0 & u_{33}
!>          \end{bmatrix}
!>          \begin{bmatrix}
!>              A_{11}^{-1} \\
!>              A_{21}^{-1} \\
!>              A_{31}^{-1}
!>          \end{bmatrix}
!>      \right)
!>      =
!>      \begin{bmatrix}
!>          1 \\
!>          0 \\
!>          0 ~.
!>      \end{bmatrix}
!>  \f}
!>  The above form of equations implied that the two systems of triangular equations have to be solved to obtain one column of the inverse matrix.<br>
!>  However, this method is fast because only back- and forward-substitution are required to solve for the column vectors after the initial factorization of the matrix \f$A\f$.<br>
!>  The most common factorization methods used for \f$A\f$ to obtain its inverse are:<br>
!>  <ol>
!>      <li>    [Pivoted LU factorization](@ref pm_matrixLUP) for general square matrices.<br>
!>      <li>    [Cholesky factorization](@ref pm_matrixChol) for Symmetric/Hermitian Positive-Definite square matrices.<br>
!>  </ol>
!>
!>  Inversion-Transposition operation
!>  =================================
!>
!>  There are also special matrix operations that mix **inversion** with Symmetric
!>  and Hermitian each having corresponding [matrix classes](@ref pm_matrixClass):**<br>
!>  <ol>
!>      <li>    **[Orthogonal transposition](@ref pm_matrixTrans::transOrth_type)**<br>
!>              A square matrix whose transpose is equal to its inverse is called an Orthogonal matrix.<br>
!>              In other words, \f$A\f$ is Orthogonal if \f$\mathbf{A}^{\up{T}} = \mathbf{A}^{-1}\f$.<br>
!>              The corresponding transposition is called [Orthogonal](@ref pm_matrixTrans::transOrth_type) denoted by the operator \f$\cdot^{\up{-T}}\f$.<br>
!>      <li>    **[Unitary transposition](@ref pm_matrixTrans::transUnit_type)**<br>
!>              A square complex matrix whose transpose is equal to its conjugate inverse is called a Unitary matrix.<br>
!>              In other words, \f$A\f$ is Unitary if \f$\mathbf{A}^{\up{T}} = {\overline{\mathbf{A}^{-1}}}\f$.<br>
!>              The corresponding transposition is called [Unitary](@ref pm_matrixTrans::transUnit_type) denoted by the operator \f$\cdot^{\up{-H}}\f$.<br>
!>  </ol>
!>
!>  \see
!>  [pm_matrixDet](@ref pm_matrixDet)<br>
!>  [pm_matrixLUP](@ref pm_matrixLUP)<br>
!>  [pm_matrixChol](@ref pm_matrixChol)<br>
!>  [pm_matrixTrans](@ref pm_matrixTrans)<br>
!>
!>  \benchmarks
!>
!>  \benchmark{setMatInv, The runtime performance of [setMatInv](@ref pm_matrixInv::setMatInv) for various methods and inverse matrix subsets.}
!>  \include{lineno} benchmark/pm_matrixInv/setMatInv/main.F90
!>  \compilefb{setMatInv}
!>  \postprocb{setMatInv}
!>  \include{lineno} benchmark/pm_matrixInv/setMatInv/main.py
!>  \visb{setMatInv}
!>  \image html benchmark/pm_matrixInv/setMatInv/benchmark.setMatInv.runtime.png width=1000
!>  \image html benchmark/pm_matrixInv/setMatInv/benchmark.setMatInv.runtime.ratio.png width=1000
!>  \moralb{setMatInv}
!>      -#  The procedures under the generic interface [setMatInv](@ref pm_matrixInv::setMatInv)
!>          use an unblocked approach to computing the matrix inverse.<br>
!>          However, specifying an upper-triangular Cholesky factor along with lower-triangle for the matrix inverse
!>          can potentially result in faster calculations as all matrix operations within the algorithm become column-major
!>          meaning that all memory access become local.<br>
!>
!>  \test
!>  [test_pm_matrixInv](@ref test_pm_matrixInv)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_matrixInv

    use pm_kind, only: SK, IK
    use pm_matrixClass, only: matrix_type
    use pm_matrixClass, only: lup, lup_type
    use pm_matrixClass, only: square, square_type
    use pm_matrixClass, only: choUpp, choUpp_type
    use pm_matrixClass, only: choLow, choLow_type
    use pm_matrixClass, only: posdefmat, posdefmat_type
    use pm_matrixClass, only: upperDiag, upperDiag_type
    use pm_matrixClass, only: lowerDiag, lowerDiag_type
    use pm_matrixClass, only: upperUnit, upperUnit_type
    use pm_matrixClass, only: lowerUnit, lowerUnit_type
    use pm_matrixSubset, only: lowDia, lowDia_type
    use pm_matrixSubset, only: uppDia, uppDia_type
    use pm_matrixSubset, only: subset_type

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_matrixInv"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request inversion operation on a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [inversion](@ref pm_matrixInv::inversion)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [inversion](@ref pm_matrixInv::inversion)<br>
    !>  [inversion_type](@ref pm_matrixInv::inversion_type)<br>
    !>  [trans_type](@ref pm_matrixTrans::trans_type)<br>
    !>
    !>  \final{inversion_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: inversion_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [inversion_type](@ref pm_matrixInv::inversion_type) that is exclusively used
    !>  to request no transpose of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [nothing](@ref pm_array::nothing)<br>
    !>  [trans](@ref pm_matrixTrans::trans)<br>
    !>
    !>  \final{inversion}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(inversion_type), parameter :: inversion = inversion_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: inversion
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the full inverse of an input matrix of general or triangular form directly or through its input the LU/Cholesky factorization.<br>
    !>
    !>  \param[in]      mat         :   The input array of rank `2` of square shape `(1 : ndim, 1 : ndim)` of,
    !>                                  <ol>
    !>                                      <li>     type `complex` of kind \CKALL,
    !>                                      <li>     type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing a matrix or a factorization of a matrix whose inverse must be computed.<br>
    !>                                  The nature of the input `mat` is determined by the optional input argument `auxil`.<br>
    !>  \param[inout]   auxil       :   The input scalar constant that can be:<br>
    !>                                  <ol>
    !>                                      <li>    The output scalar of type `integer` of default kind \IK that is `0` is **if and only if** the algorithm succeeds to compute the inverse of the input general matrix.<br>
    !>                                              No assumption is made about the input `mat` other than being given in full square form.<br>
    !>                                              See the `info` argument of [setMatLUP](@ref pm_matrixLUP::setMatLUP) for the possible meanings of a non-zero output value for `auxil`.<br>
    !>                                      <li>    The output scalar of the same type and kind as the input `mat` containing the determinant of the output inverse matrix `inv`.<br>
    !>                                              No assumption is made about the input `mat` other than being given in full square form.<br>
    !>                                      <li>    The input scalar constant [upperDiag](@ref pm_matrixClass::upperDiag) implying that the input `mat` is an upper-diagonal triangular matrix whose inverse is to be computed.<br>
    !>                                      <li>    The input scalar constant [lowerDiag](@ref pm_matrixClass::lowerDiag) implying that the input `mat` is an lower-diagonal triangular matrix whose inverse is to be computed.<br>
    !>                                      <li>    The input `contiguous` vector of shape `(1 : ndim)` containing the row permutations corresponding to the LU factorization of the matrix whose inverse is to be computed.<br>
    !>                                              This argument is automatically returned as `rperm` along with the LU factorization from the generic interfaces of [pm_matrixLUP](@ref pm_matrixLUP).<br>
    !>                                              This value must be specified **if and only if** the input argument `mat` contains the LU factorization of the matrix whose inverse must be computed.<br>
    !>                                      <li>    The input scalar constant [choUpp](@ref pm_matrixClass::choUpp) implying that only the upper-diagonal Cholesky factorization of the positive-definite matrix whose inverse is to be computed.<br>
    !>                                      <li>    The input scalar constant [choLow](@ref pm_matrixClass::choLow) implying that only the lower-diagonal Cholesky factorization of the positive-definite matrix whose inverse is to be computed.<br>
    !>                                  </ol>
    !>                                  (**optional**. If missing, the input matrix is assumed to be a general square matrix whose inverse is to be computed. If the algorithm fails, the program will halt by calling `error stop`.)
    !>
    !>  \return
    !>  `inv`                       :   The output matrix of the same type, kind, and shape as the input `mat` containing the **full inverse matrix**.<br>
    !>
    !>  \interface{getMatInv}
    !>  \code{.F90}
    !>
    !>      use pm_matrixInv, only: getMatInv, choUpp, choLow, upperDiag, lowerDiag
    !>
    !>      inv(1:ndim, 1:ndim) = getMatInv(mat(1:ndim, 1:ndim))
    !>      inv(1:ndim, 1:ndim) = getMatInv(mat(1:ndim, 1:ndim), auxil)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(mat, 1) == size(mat, 2))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(size(auxil) == shape(mat))` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>  The procedures of this generic interface are always `impure` when the argument `auxil` has `intent(out)`.<br>
    !>
    !>  \see
    !>  [getMatInv](@ref pm_matrixInv::getMatInv)<br>
    !>  [setMatInv](@ref pm_matrixInv::setMatInv)<br>
    !>  [getMatChol](@ref pm_matrixChol::getMatChol)<br>
    !>  [setMatChol](@ref pm_matrixChol::setMatChol)<br>
    !>  [setMatLUP](@ref pm_matrixLUP::setMatLUP)<br>
    !>  [Intel Fortran LAPACK documentation](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-1/matrix-inversion-lapack-computational-routines.html)<br>
    !>
    !>  \lapack{3.11}
    !>  `SGETRI`, `DGETRI`, `CGETRI`, and `ZGETRI`.<br>
    !>  `SPFTRI`, `DPFTRI`, `CPFTRI`, and `ZPFTRI`.<br>
    !>  `SPOTRI`, `DPOTRI`, `CPOTRI`, and `ZPOTRI`.<br>
    !>  `STRTRI`, `DTRTRI`, `CTRTRI`, and `ZTRTRI`.<br>
    !>
    !>  \example{getMatInv}
    !>  \include{lineno} example/pm_matrixInv/getMatInv/main.F90
    !>  \compilef{getMatInv}
    !>  \output{getMatInv}
    !>  \include{lineno} example/pm_matrixInv/getMatInv/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixInv](@ref test_pm_matrixInv)<br>
    !>
    !>  \final{getMatInv}
    !>
    !>  \author
    !>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! implicit Def.

    interface getMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getMatInvDef_IMP_CK5(mat) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDef_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if CK4_ENABLED
    impure module function getMatInvDef_IMP_CK4(mat) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDef_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if CK3_ENABLED
    impure module function getMatInvDef_IMP_CK3(mat) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDef_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if CK2_ENABLED
    impure module function getMatInvDef_IMP_CK2(mat) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDef_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if CK1_ENABLED
    impure module function getMatInvDef_IMP_CK1(mat) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDef_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMatInvDef_IMP_RK5(mat) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDef_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMatInvDef_IMP_RK4(mat) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDef_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMatInvDef_IMP_RK3(mat) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDef_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMatInvDef_IMP_RK2(mat) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDef_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMatInvDef_IMP_RK1(mat) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDef_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit Det.

    interface getMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getMatInvDet_IMP_CK5(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDet_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(out)                   :: auxil
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if CK4_ENABLED
    impure module function getMatInvDet_IMP_CK4(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDet_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(out)                   :: auxil
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if CK3_ENABLED
    impure module function getMatInvDet_IMP_CK3(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDet_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(out)                   :: auxil
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if CK2_ENABLED
    impure module function getMatInvDet_IMP_CK2(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDet_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(out)                   :: auxil
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if CK1_ENABLED
    impure module function getMatInvDet_IMP_CK1(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDet_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(out)                   :: auxil
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMatInvDet_IMP_RK5(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDet_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(out)                   :: auxil
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if RK4_ENABLED
    impure module function getMatInvDet_IMP_RK4(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDet_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(out)                   :: auxil
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if RK3_ENABLED
    impure module function getMatInvDet_IMP_RK3(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDet_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(out)                   :: auxil
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if RK2_ENABLED
    impure module function getMatInvDet_IMP_RK2(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDet_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(out)                   :: auxil
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

#if RK1_ENABLED
    impure module function getMatInvDet_IMP_RK1(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvDet_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(out)                   :: auxil
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit Inf.

    interface getMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    impure module function getMatInvInf_IMP_CK5(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvInf_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(out)                   :: auxil
    end function
#endif

#if CK4_ENABLED
    impure module function getMatInvInf_IMP_CK4(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvInf_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(out)                   :: auxil
    end function
#endif

#if CK3_ENABLED
    impure module function getMatInvInf_IMP_CK3(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvInf_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(out)                   :: auxil
    end function
#endif

#if CK2_ENABLED
    impure module function getMatInvInf_IMP_CK2(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvInf_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(out)                   :: auxil
    end function
#endif

#if CK1_ENABLED
    impure module function getMatInvInf_IMP_CK1(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvInf_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(out)                   :: auxil
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    impure module function getMatInvInf_IMP_RK5(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvInf_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(out)                   :: auxil
    end function
#endif

#if RK4_ENABLED
    impure module function getMatInvInf_IMP_RK4(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvInf_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(out)                   :: auxil
    end function
#endif

#if RK3_ENABLED
    impure module function getMatInvInf_IMP_RK3(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvInf_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(out)                   :: auxil
    end function
#endif

#if RK2_ENABLED
    impure module function getMatInvInf_IMP_RK2(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvInf_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(out)                   :: auxil
    end function
#endif

#if RK1_ENABLED
    impure module function getMatInvInf_IMP_RK1(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvInf_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(out)                   :: auxil
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit upperDiag.

    interface getMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatInvCUD_IMP_CK5(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUD_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperDiag_type), intent(in)                    :: auxil
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatInvCUD_IMP_CK4(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUD_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperDiag_type), intent(in)                    :: auxil
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatInvCUD_IMP_CK3(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUD_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperDiag_type), intent(in)                    :: auxil
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatInvCUD_IMP_CK2(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUD_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperDiag_type), intent(in)                    :: auxil
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatInvCUD_IMP_CK1(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUD_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperDiag_type), intent(in)                    :: auxil
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatInvCUD_IMP_RK5(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUD_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperDiag_type), intent(in)                    :: auxil
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatInvCUD_IMP_RK4(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUD_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperDiag_type), intent(in)                    :: auxil
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatInvCUD_IMP_RK3(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUD_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperDiag_type), intent(in)                    :: auxil
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatInvCUD_IMP_RK2(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUD_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperDiag_type), intent(in)                    :: auxil
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatInvCUD_IMP_RK1(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUD_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperDiag_type), intent(in)                    :: auxil
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit lowerDiag.

    interface getMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatInvCLD_IMP_CK5(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLD_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerDiag_type), intent(in)                    :: auxil
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatInvCLD_IMP_CK4(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLD_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerDiag_type), intent(in)                    :: auxil
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatInvCLD_IMP_CK3(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLD_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerDiag_type), intent(in)                    :: auxil
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatInvCLD_IMP_CK2(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLD_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerDiag_type), intent(in)                    :: auxil
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatInvCLD_IMP_CK1(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLD_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerDiag_type), intent(in)                    :: auxil
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatInvCLD_IMP_RK5(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLD_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerDiag_type), intent(in)                    :: auxil
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatInvCLD_IMP_RK4(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLD_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerDiag_type), intent(in)                    :: auxil
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatInvCLD_IMP_RK3(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLD_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerDiag_type), intent(in)                    :: auxil
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatInvCLD_IMP_RK2(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLD_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerDiag_type), intent(in)                    :: auxil
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatInvCLD_IMP_RK1(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLD_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerDiag_type), intent(in)                    :: auxil
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit upperUnit.

    interface getMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatInvCUU_IMP_CK5(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUU_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperUnit_type), intent(in)                    :: auxil
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatInvCUU_IMP_CK4(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUU_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperUnit_type), intent(in)                    :: auxil
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatInvCUU_IMP_CK3(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUU_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperUnit_type), intent(in)                    :: auxil
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatInvCUU_IMP_CK2(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUU_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperUnit_type), intent(in)                    :: auxil
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatInvCUU_IMP_CK1(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUU_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperUnit_type), intent(in)                    :: auxil
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatInvCUU_IMP_RK5(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUU_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperUnit_type), intent(in)                    :: auxil
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatInvCUU_IMP_RK4(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUU_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperUnit_type), intent(in)                    :: auxil
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatInvCUU_IMP_RK3(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUU_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperUnit_type), intent(in)                    :: auxil
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatInvCUU_IMP_RK2(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUU_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperUnit_type), intent(in)                    :: auxil
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatInvCUU_IMP_RK1(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCUU_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(upperUnit_type), intent(in)                    :: auxil
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit lowerUnit.

    interface getMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatInvCLU_IMP_CK5(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLU_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerUnit_type), intent(in)                    :: auxil
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatInvCLU_IMP_CK4(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLU_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerUnit_type), intent(in)                    :: auxil
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatInvCLU_IMP_CK3(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLU_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerUnit_type), intent(in)                    :: auxil
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatInvCLU_IMP_CK2(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLU_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerUnit_type), intent(in)                    :: auxil
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatInvCLU_IMP_CK1(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLU_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerUnit_type), intent(in)                    :: auxil
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatInvCLU_IMP_RK5(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLU_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerUnit_type), intent(in)                    :: auxil
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatInvCLU_IMP_RK4(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLU_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerUnit_type), intent(in)                    :: auxil
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatInvCLU_IMP_RK3(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLU_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerUnit_type), intent(in)                    :: auxil
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatInvCLU_IMP_RK2(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLU_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerUnit_type), intent(in)                    :: auxil
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatInvCLU_IMP_RK1(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCLU_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(lowerUnit_type), intent(in)                    :: auxil
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit lup.

    interface getMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatInvLUP_IMP_CK5(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvLUP_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatInvLUP_IMP_CK4(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvLUP_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatInvLUP_IMP_CK3(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvLUP_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatInvLUP_IMP_CK2(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvLUP_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatInvLUP_IMP_CK1(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvLUP_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatInvLUP_IMP_RK5(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvLUP_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatInvLUP_IMP_RK4(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvLUP_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatInvLUP_IMP_RK3(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvLUP_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatInvLUP_IMP_RK2(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvLUP_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatInvLUP_IMP_RK1(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvLUP_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit choUpp.

    interface getMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatInvCCU_IMP_CK5(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCU_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choUpp_type)   , intent(in)                    :: auxil
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatInvCCU_IMP_CK4(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCU_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choUpp_type)   , intent(in)                    :: auxil
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatInvCCU_IMP_CK3(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCU_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choUpp_type)   , intent(in)                    :: auxil
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatInvCCU_IMP_CK2(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCU_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choUpp_type)   , intent(in)                    :: auxil
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatInvCCU_IMP_CK1(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCU_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choUpp_type)   , intent(in)                    :: auxil
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatInvCCU_IMP_RK5(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCU_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choUpp_type)   , intent(in)                    :: auxil
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatInvCCU_IMP_RK4(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCU_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choUpp_type)   , intent(in)                    :: auxil
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatInvCCU_IMP_RK3(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCU_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choUpp_type)   , intent(in)                    :: auxil
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatInvCCU_IMP_RK2(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCU_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choUpp_type)   , intent(in)                    :: auxil
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatInvCCU_IMP_RK1(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCU_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choUpp_type)   , intent(in)                    :: auxil
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit choLow.

    interface getMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatInvCCL_IMP_CK5(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCL_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choLow_type)   , intent(in)                    :: auxil
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatInvCCL_IMP_CK4(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCL_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choLow_type)   , intent(in)                    :: auxil
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatInvCCL_IMP_CK3(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCL_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choLow_type)   , intent(in)                    :: auxil
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatInvCCL_IMP_CK2(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCL_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choLow_type)   , intent(in)                    :: auxil
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatInvCCL_IMP_CK1(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCL_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)                                        :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choLow_type)   , intent(in)                    :: auxil
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatInvCCL_IMP_RK5(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCL_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choLow_type)   , intent(in)                    :: auxil
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatInvCCL_IMP_RK4(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCL_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choLow_type)   , intent(in)                    :: auxil
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatInvCCL_IMP_RK3(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCL_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choLow_type)   , intent(in)                    :: auxil
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatInvCCL_IMP_RK2(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCL_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choLow_type)   , intent(in)                    :: auxil
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatInvCCL_IMP_RK1(mat, auxil) result(inv)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatInvCCL_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)                                           :: inv(size(mat, 1, IK), size(mat, 2, IK))
        type(choLow_type)   , intent(in)                    :: auxil
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the full inverse of a general or triangular matrix or a subset of the
    !>  inverse of a positive-definite matrix complementary to its specified Cholesky factorization subset.<br>
    !>
    !>  \param[inout]   inv         :   The output or input/output matrix of the same type, kind, and shape as the input `mat` containing the full inverse matrix or a subset of it:<br>
    !>                                  <ol>
    !>                                      <li>    If the specified input `mat` is a,<br>
    !>                                              <ol>
    !>                                                  <li>    [lower-diagonal triangular matrix](@ref pm_matrixClass::lowerDiag_type),
    !>                                                  <li>    [upper-diagonal triangular matrix](@ref pm_matrixClass::upperDiag_type),
    !>                                                  <li>    [lower-unit-diagonal triangular matrix](@ref pm_matrixClass::lowerUnit_type),
    !>                                                  <li>    [upper-unit-diagonal triangular matrix](@ref pm_matrixClass::upperUnit_type),
    !>                                                  <li>    [LUP factorization](@ref pm_matrixLUP) of a general square matrix,
    !>                                                  <li>    [Cholesky factorization](@ref pm_matrixChol) of a positive-definite square matrix while the optional input argument `subset` is missing,
    !>                                              </ol>
    !>                                              then `inv` has `intent(out)` and shall contain the **full inverse of the matrix** corresponding to the specified `mat`.<br>
    !>                                      <li>    If the specified input `mat` is a [Cholesky factorization](@ref pm_matrixChol), and the input argument `subset` is present, then `inv` has `intent(inout)`.<br>
    !>                                              On output, only the specified `subset` of `inv` shall be updated with the corresponding inverse matrix triangle while the opposite triangle of `inv` remains untouched.<br>
    !>                                              This storage scheme for the inverse matrix is very useful for efficient packing both the Cholesky factorization and
    !>                                              inverse of a matrix in a single `contiguous` almost-squqare-shaped matrix (e.g., `inv(1 : ndim, 0 : ndim)`).<br>
    !>                                  </ol>
    !>  \param[in]      mat         :   The input array of rank `2` of square shape `(1 : ndim, 1 : ndim)` of,
    !>                                  <ol>
    !>                                      <li>     type `complex` of kind \CKALL,
    !>                                      <li>     type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing a.<br>
    !>                                  <ol>
    !>                                      <li>    [upper-diagonal triangular matrix](@ref pm_matrixClass::upperDiag_type),
    !>                                      <li>    [lower-diagonal triangular matrix](@ref pm_matrixClass::lowerDiag_type),
    !>                                      <li>    [upper-unit-diagonal triangular matrix](@ref pm_matrixClass::upperUnit_type),
    !>                                      <li>    [lower-unit-diagonal triangular matrix](@ref pm_matrixClass::lowerUnit_type),
    !>                                      <li>    [LUP factorization](@ref pm_matrixLUP) of a general square matrix,
    !>                                      <li>    [lower Cholesky factorization](@ref pm_matrixClass::choLow_type) of a positive-definite matrix,
    !>                                      <li>    [upper Cholesky factorization](@ref pm_matrixClass::choUpp_type) of a positive-definite matrix,
    !>                                  </ol>
    !>                                  whose inverse must be computed.<br>
    !>                                  The class of the input `mat` is determined by the input argument `auxil`.<br>
    !>  \param[inout]   auxil       :   The input scalar constant that can be:<br>
    !>                                  <ol>
    !>                                      <li>    The input scalar constant [upperDiag](@ref pm_matrixClass::upperDiag) implying that the input `mat` is an upper-diagonal triangular matrix whose inverse is to be computed.<br>
    !>                                      <li>    The input scalar constant [lowerDiag](@ref pm_matrixClass::lowerDiag) implying that the input `mat` is an lower-diagonal triangular matrix whose inverse is to be computed.<br>
    !>                                      <li>    The input `contiguous` vector of shape `(1 : ndim)` containing the row permutations corresponding to the LU factorization of the matrix whose inverse is to be computed.<br>
    !>                                              This argument is automatically returned as `rperm` along with the LU factorization from the generic interfaces of [pm_matrixLUP](@ref pm_matrixLUP).<br>
    !>                                              This value must be specified **if and only if** the input argument `mat` contains the LU factorization of the matrix whose inverse must be computed.<br>
    !>                                      <li>    The input scalar constant [choUpp](@ref pm_matrixClass::choUpp) implying that only the upper-diagonal Cholesky factorization of the positive-definite matrix whose inverse is to be computed.<br>
    !>                                      <li>    The input scalar constant [choLow](@ref pm_matrixClass::choLow) implying that only the lower-diagonal Cholesky factorization of the positive-definite matrix whose inverse is to be computed.<br>
    !>                                  </ol>
    !>  \param[in]      subset      :   The input scalar that can be either,
    !>                                  <ol>
    !>                                      <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) implying that only
    !>                                              the upper-diagonal triangular block of `inv` should be returned on output
    !>                                              while the lower triangular part of `inv` is not referenced within the algorithm.
    !>                                      <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) implying that
    !>                                              the lower-diagonal triangular block of `inv` should be returned on output
    !>                                              while the lower triangular part of `inv` is not referenced within the algorithm.
    !>                                  </ol>
    !>                                  This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.<br>
    !>                                  (**optional**. It can be present **only if** `auxil` is set to either [choLow](@ref pm_matrixClass::choLow) or [choUpp](@ref pm_matrixClass::choUpp).)
    !>
    !>  \interface{setMatInv}
    !>  \code{.F90}
    !>
    !>      use pm_matrixInv, only: setMatInv
    !>      use pm_matrixInv, only: choLow, choUpp
    !>      use pm_matrixInv, only: upperDiag, lowerDiag
    !>
    !>      call setMatInv(inv(1:ndim, 1:ndim), mat(1:ndim, 1:ndim), auxil)
    !>      call setMatInv(inv(1:ndim, 1:ndim), mat(1:ndim, 1:ndim), auxil, subset) ! only if subset = choLow, choUpp.
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(mat, 1) == size(mat, 2))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(shape(inv) == shape(mat))` must hold for the corresponding input arguments.<br>
    !>  The condition `size(auxil) == shape(mat, 1, IK)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getMatInv](@ref pm_matrixInv::getMatInv)<br>
    !>  [setMatInv](@ref pm_matrixInv::setMatInv)<br>
    !>  [getMatChol](@ref pm_matrixChol::getMatChol)<br>
    !>  [setMatChol](@ref pm_matrixChol::setMatChol)<br>
    !>  [setMatLUP](@ref pm_matrixLUP::setMatLUP)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [Intel Fortran LAPACK documentation](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-fortran/2023-1/matrix-inversion-lapack-computational-routines.html)<br>
    !>
    !>  \lapack{3.11}
    !>  `SGETRI`, `DGETRI`, `CGETRI`, and `ZGETRI`.<br>
    !>  `SPFTRI`, `DPFTRI`, `CPFTRI`, and `ZPFTRI`.<br>
    !>  `SPOTRI`, `DPOTRI`, `CPOTRI`, and `ZPOTRI`.<br>
    !>  `STRTRI`, `DTRTRI`, `CTRTRI`, and `ZTRTRI`.<br>
    !>
    !>  \example{setMatInv}
    !>  \include{lineno} example/pm_matrixInv/setMatInv/main.F90
    !>  \compilef{setMatInv}
    !>  \output{setMatInv}
    !>  \include{lineno} example/pm_matrixInv/setMatInv/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixInv](@ref test_pm_matrixInv)<br>
    !>
    !>  \todo
    !>  \phigh
    !>  The functionality of this generic interface must be extended to fully dispatch to LAPACK inverse matrix routines where appropriate.<br>
    !>
    !>  \final{setMatInv}
    !>
    !>  \author
    !>  \AmirShahmoradi, Apr 21, 2017, 1:54 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! implicit upperDiag invFUL.

    interface setMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInvCUD_IMP_CK5(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUD_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(upperDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInvCUD_IMP_CK4(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUD_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(upperDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInvCUD_IMP_CK3(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUD_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(upperDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInvCUD_IMP_CK2(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUD_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(upperDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInvCUD_IMP_CK1(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUD_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(upperDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInvCUD_IMP_RK5(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUD_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(upperDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInvCUD_IMP_RK4(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUD_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(upperDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInvCUD_IMP_RK3(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUD_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(upperDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInvCUD_IMP_RK2(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUD_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(upperDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInvCUD_IMP_RK1(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUD_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(upperDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit lowerDiag invFUL.

    interface setMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInvCLD_IMP_CK5(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLD_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(lowerDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInvCLD_IMP_CK4(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLD_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(lowerDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInvCLD_IMP_CK3(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLD_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(lowerDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInvCLD_IMP_CK2(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLD_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(lowerDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInvCLD_IMP_CK1(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLD_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(lowerDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInvCLD_IMP_RK5(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLD_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(lowerDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInvCLD_IMP_RK4(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLD_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(lowerDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInvCLD_IMP_RK3(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLD_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(lowerDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInvCLD_IMP_RK2(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLD_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(lowerDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInvCLD_IMP_RK1(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLD_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(lowerDiag_type), intent(in)                    :: auxil
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit upperUnit invFUL.

    interface setMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInvCUU_IMP_CK5(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUU_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(upperUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInvCUU_IMP_CK4(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUU_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(upperUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInvCUU_IMP_CK3(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUU_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(upperUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInvCUU_IMP_CK2(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUU_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(upperUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInvCUU_IMP_CK1(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUU_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(upperUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInvCUU_IMP_RK5(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUU_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(upperUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInvCUU_IMP_RK4(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUU_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(upperUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInvCUU_IMP_RK3(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUU_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(upperUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInvCUU_IMP_RK2(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUU_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(upperUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInvCUU_IMP_RK1(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCUU_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(upperUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit lowerUnit invFUL.

    interface setMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInvCLU_IMP_CK5(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLU_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(lowerUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInvCLU_IMP_CK4(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLU_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(lowerUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInvCLU_IMP_CK3(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLU_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(lowerUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInvCLU_IMP_CK2(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLU_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(lowerUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInvCLU_IMP_CK1(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLU_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(lowerUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInvCLU_IMP_RK5(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLU_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(lowerUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInvCLU_IMP_RK4(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLU_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(lowerUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInvCLU_IMP_RK3(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLU_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(lowerUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInvCLU_IMP_RK2(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLU_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(lowerUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInvCLU_IMP_RK1(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCLU_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(lowerUnit_type), intent(in)                    :: auxil
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit lup invFUL.

    interface setMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInvLUP_IMP_CK5(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvLUP_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(out)   , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInvLUP_IMP_CK4(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvLUP_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(out)   , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInvLUP_IMP_CK3(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvLUP_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(out)   , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInvLUP_IMP_CK2(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvLUP_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(out)   , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInvLUP_IMP_CK1(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvLUP_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(out)   , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInvLUP_IMP_RK5(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvLUP_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)   , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInvLUP_IMP_RK4(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvLUP_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)   , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInvLUP_IMP_RK3(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvLUP_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)   , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInvLUP_IMP_RK2(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvLUP_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)   , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInvLUP_IMP_RK1(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvLUP_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)   , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        integer(IK)         , intent(in)    , contiguous    :: auxil(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit choUpp invFUL.

    interface setMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInvCCU_FUL_IMP_CK5(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_FUL_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(out)   , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInvCCU_FUL_IMP_CK4(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_FUL_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(out)   , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInvCCU_FUL_IMP_CK3(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_FUL_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(out)   , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInvCCU_FUL_IMP_CK2(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_FUL_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(out)   , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInvCCU_FUL_IMP_CK1(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_FUL_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(out)   , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInvCCU_FUL_IMP_RK5(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_FUL_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)   , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInvCCU_FUL_IMP_RK4(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_FUL_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)   , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInvCCU_FUL_IMP_RK3(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_FUL_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)   , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInvCCU_FUL_IMP_RK2(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_FUL_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)   , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInvCCU_FUL_IMP_RK1(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_FUL_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)   , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit choLow invFUL.

    interface setMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInvCCL_FUL_IMP_CK5(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_FUL_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(out)   , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInvCCL_FUL_IMP_CK4(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_FUL_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(out)   , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInvCCL_FUL_IMP_CK3(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_FUL_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(out)   , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInvCCL_FUL_IMP_CK2(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_FUL_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(out)   , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInvCCL_FUL_IMP_CK1(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_FUL_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(out)   , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInvCCL_FUL_IMP_RK5(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_FUL_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(out)   , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInvCCL_FUL_IMP_RK4(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_FUL_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(out)   , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInvCCL_FUL_IMP_RK3(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_FUL_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(out)   , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInvCCL_FUL_IMP_RK2(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_FUL_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(out)   , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInvCCL_FUL_IMP_RK1(inv, mat, auxil)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_FUL_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(out)   , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit choUpp invXLD.

    interface setMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInvCCU_XLD_IMP_CK5(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_XLD_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInvCCU_XLD_IMP_CK4(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_XLD_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInvCCU_XLD_IMP_CK3(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_XLD_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInvCCU_XLD_IMP_CK2(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_XLD_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInvCCU_XLD_IMP_CK1(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_XLD_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInvCCU_XLD_IMP_RK5(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_XLD_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInvCCU_XLD_IMP_RK4(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_XLD_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInvCCU_XLD_IMP_RK3(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_XLD_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInvCCU_XLD_IMP_RK2(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_XLD_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInvCCU_XLD_IMP_RK1(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_XLD_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit choLow invXLD.

    interface setMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInvCCL_XLD_IMP_CK5(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_XLD_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInvCCL_XLD_IMP_CK4(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_XLD_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInvCCL_XLD_IMP_CK3(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_XLD_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInvCCL_XLD_IMP_CK2(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_XLD_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInvCCL_XLD_IMP_CK1(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_XLD_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInvCCL_XLD_IMP_RK5(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_XLD_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInvCCL_XLD_IMP_RK4(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_XLD_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInvCCL_XLD_IMP_RK3(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_XLD_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInvCCL_XLD_IMP_RK2(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_XLD_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInvCCL_XLD_IMP_RK1(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_XLD_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(lowDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit choUpp invUXD.

    interface setMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInvCCU_UXD_IMP_CK5(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_UXD_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInvCCU_UXD_IMP_CK4(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_UXD_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInvCCU_UXD_IMP_CK3(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_UXD_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInvCCU_UXD_IMP_CK2(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_UXD_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInvCCU_UXD_IMP_CK1(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_UXD_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInvCCU_UXD_IMP_RK5(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_UXD_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInvCCU_UXD_IMP_RK4(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_UXD_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInvCCU_UXD_IMP_RK3(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_UXD_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInvCCU_UXD_IMP_RK2(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_UXD_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInvCCU_UXD_IMP_RK1(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCU_UXD_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choUpp_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! implicit choLow invUXD.

    interface setMatInv

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatInvCCL_UXD_IMP_CK5(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_UXD_IMP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatInvCCL_UXD_IMP_CK4(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_UXD_IMP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatInvCCL_UXD_IMP_CK3(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_UXD_IMP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatInvCCL_UXD_IMP_CK2(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_UXD_IMP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatInvCCL_UXD_IMP_CK1(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_UXD_IMP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(inout) , contiguous    :: inv(:,:)
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatInvCCL_UXD_IMP_RK5(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_UXD_IMP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatInvCCL_UXD_IMP_RK4(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_UXD_IMP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatInvCCL_UXD_IMP_RK3(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_UXD_IMP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatInvCCL_UXD_IMP_RK2(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_UXD_IMP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatInvCCL_UXD_IMP_RK1(inv, mat, auxil, subset)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatInvCCL_UXD_IMP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: inv(:,:)
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        type(choLow_type)   , intent(in)                    :: auxil
        type(uppDia_type)   , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_matrixInv ! LCOV_EXCL_LINE