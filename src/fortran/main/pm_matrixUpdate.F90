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
!>  This module contains procedures and generic interfaces relevant to arbitrary-rank updates to vectors, general matrices, or
!>  Symmetric/Hermitian triangular matrices of type `integer`, `complex`, and `real` of arbitrary type-kind parameters.<br>
!>
!>  \details
!>
!>  Rank-1 matrix update
!>  --------------------
!>
!>  Let \f$\ms{matA}\f$ be a \f$n\times n\f$ matrix and \f$\ms{vecA}\f$ and \f$\ms{vecB}\f$ two \f$n\times 1\f$ column vectors.<br>
!>  The following transformation is called a **rank one update** to \f$\ms{matA}\f$,<br>
!>  \f{equation}{
!>      \ms{matA} := \ms{matA} + \alpha \ms{vecA}~\ms{vecB}^\mathrm{T} ~,
!>  \f}
!>  where \f$\alpha\f$ is an arbitrary scalar.<br>
!>  For complex matrices, the transformation can be also done as the following,
!>  \f{equation}{
!>      \ms{matA} := \ms{matA} + \alpha \ms{vecA}~\ms{vecB}^\mathrm{H} ~,
!>  \f}
!>  where \f$\mathrm{H}\f$ denotes the [Hermitian](@ref pm_matrixTrans) (or conjugate or Adjoint) transpose operation.<br>
!>
!>  \note
!>  The *rank one* name of the transformation stems from the fact that the rank of the \f$n \times n\f$ matrix \f$\ms{vecA}~\ms{vecB}^{T}\f$ is equal to \f$1\f$.<br>
!>  This is because a single vector, \f$\ms{vecA}\f$, spans all the columns of \f$\ms{vecA}~\ms{vecB}^{T}\f$.
!>
!>  <b>Motivation: Sherman-Morrison formula</b><br>
!>  The Rank One Update transform offers a fast method of indirectly computing the inverse of an update matrix.<br>
!>  In mathematics, in particular linear algebra, the Sherman–Morrison formula, named after Jack Sherman and Winifred J. Morrison,
!>  computes the inverse of the sum of an invertible matrix \f$\ms{matA}\f$ and the outer product, \f$\ms{vecA}~\ms{vecB}^{\mathrm{T}}\f$,
!>  of vectors \f$\ms{vecA}\f$ and \f$\ms{vecB}\f$.<br>
!>  Suppose \f$\ms{matA} \in \mathbb{R}^{n\times n}\f$ is an invertible square matrix and \f$\ms{vecA}, \ms{vecB}\in \mathbb{R}^{n}\f$ are column vectors.<br>
!>  Then \f$\ms{matA} + \ms{vecA}~\ms{vecB}^{\mathrm{T}}\f$ is invertible <b>if and only if</b> \f$1 + \ms{vecB}^{\mathrm{T}} \ms{matA}^{-1} \ms{vecA} \neq 0\f$.<br>
!>  In this case,<br>
!>  \f{equation}{
!>      \left(\ms{matA} + \ms{vecA}~\ms{vecB}^{\textsf{T}}\right)^{-1} =
!>      \ms{matA}^{-1} - \frac{\ms{matA}^{-1} \ms{vecA}~\ms{vecB}^{\textsf{T}}\ms{matA}^{-1}} {1 + \ms{vecB}^{\textsf{T}}\ms{matA}^{-1}\ms{vecA}} ~.
!>  \f}
!>  Here, \f$\ms{vecA}~\ms{vecB}^{\textsf{T}}\f$ is thSe outer product of two vectors \f$\ms{vecA}\f$ and \f$\ms{vecB}\f$.<br>
!>  The Sherman-Morrison formula offers computationally efficient way of computing the inverse of an updated matrix when the original inverse \f$\ms{matA}^{-1}\f$ is already available.<br>
!>  The computational complexity of a full recomputation of the inverse is \f$n^{3}\f$ as opposed to \f$n^{2}\f$ via the Sherman-Morrison formula.<br>
!>  The difference in computational complexity is particularly relevant for high-dimensional matrices.<br>
!>
!>  Triangular matrix update
!>  ------------------------
!>
!>  The procedures under the generic interface [setMatUpdateTriang](@ref pm_matrixUpdate::setMatUpdateTriang)
!>  return the result of arbitrary-rank Symmetric/Hermitian updates to triangular matrices in one of the following forms,
!>  \f{align*}{
!>      & \ms{mat} \leftarrow \alpha \ms{matA}~~ \ms{matA}^T + \beta \ms{mat} \\
!>      & \ms{mat} \leftarrow \alpha \ms{matA}^T \ms{matA}~~ + \beta \ms{mat} \\
!>      & \ms{mat} \leftarrow \alpha \ms{matA}~~ \ms{matA}^H + \beta \ms{mat} \\
!>      & \ms{mat} \leftarrow \alpha \ms{matA}^H \ms{matA}~~ + \beta \ms{mat}
!>  \f}
!>  where \f$\cdot^T\f$ and \f$\cdot^H\f$ represent Symmetric and Hermitian transpositions respectively.<br>
!>  Because the output matrix \f$\ms{mat}\f$ is Symmetric/Hermitian, only the upper or the lower triangles need to be computed.<br>
!>
!>  \see
!>  [netlib::LAPACK](http://netlib.org/lapack/)<br>
!>  The IBM Engineering and Scientific Subroutine Library.<br>
!>  Developer Reference for Intel® oneAPI Math Kernel Library - Fortran.<br>
!>  Dongarra, J. J.; DuCroz, J.; Hammarling, S.; Hanson, R. J. March 1988. *An Extended Set of Fortran Basic Linear Algebra Subprograms.* ACM Transactions on Mathematical Software , 14(1):1–17.<br>
!>
!>  \remark
!>  The generic interfaces of this module dispatch to optimized BLAS libraries where available.<br>
!>  Otherwise, reference implementations of BLAS algorithms are used.<br>
!>
!>  \warning
!>  This module is still a work in progress.<br>
!>  However, all available interfaces and functions of this module are fully functional and tested.<br>
!>  But the interfaces may change over time as more functionalities are added.<br>
!>  In particular, the generic interfaces [setMatUpdateR1](@ref pm_matrixUpdate::setMatUpdateR1)
!>  and [setMatUpdateTriang](@ref pm_matrixUpdate::setMatUpdateTriang) will be merged with
!>  [setMatUpdateR1](@ref pm_matrixUpdate::setMatUpdateR1) in future library releases.<br>
!>
!>  \lapack{3.11}
!>  `SSHRK`, `DSHRK`, `CSHRK`, `ZSHRK`.<br>
!>  `SGER`, `DGER`, `CGERU`, `ZGERU`, `CGERC`, and `ZGERC`.<br>
!>  `SSYRK`, `DSYRK`, `CSYRK`, `ZSYRK`, `CHERK`, and `ZHERK`.<br>
!>  Notably, the interfaces are also extended to support matrices of type `integer` of arbitrary kinds.<br>
!>
!>  \test
!>  [test_pm_matrixUpdateR](@ref test_pm_matrixUpdateR)
!>
!>  \todo
!>  \pmed
!>  A benchmark comparison of the procedures of this module with the default BLAS/LAPACK implementation would be informative.<br>
!>
!>  \todo
!>  \pmed
!>  A benchmark comparison of the procedures of this module with the default BLAS/LAPACK implementation would be informative.<br>
!>
!>  \finmain
!>
!>  \author
!>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX
!>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_matrixUpdate

    use pm_kind, only: SK, IK, LK
    use pm_matrixSubset, only: uppDia, uppDia_type
    use pm_matrixSubset, only: lowDia, lowDia_type
    use pm_matrixClass, only: symmetric, symmetric_type
    use pm_matrixClass, only: hermitian, hermitian_type
    use pm_matrixTrans, only: transSymm, transSymm_type
    use pm_matrixTrans, only: transHerm, transHerm_type
    use pm_matrixTrans, only: trans, trans_type
    use pm_array, only: nothing, nothing_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_matrixUpdate"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the result of arbitrary-rank Symmetric/Hermitian updates
    !>  to triangular matrices of type `integer`, `complex`, and `real` of arbitrary type-kind parameters.
    !>
    !>  \details
    !>  see the documentation of [pm_matrixUpdate](@ref pm_matrixUpdate) for more details.<br>
    !>
    !>  \param[inout]   mat         :   The input/output `contiguous` matrix of arbitrary shape `(:,:)` of the same type and kind as `matA`,
    !>                                  containing the subset of matrix \f$\ms{mat}\f$ in the triangular matrix update.<br>
    !>  \param[in]      subset      :   The input scalar that can be either,
    !>                                  <ol>
    !>                                      <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) implying that
    !>                                              the upper-diagonal triangular block of `mat` should be used for storage and updating while the lower triangle remains intact.
    !>                                      <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) implying that
    !>                                              the lower-diagonal triangular block of `mat` should be used for storage and updating while the upper triangle remains intact.
    !>                                  </ol>
    !>                                  This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.
    !>  \param[in]      matA        :   The input `contiguous` matrix of arbitrary shape `(:,:)` of,
    !>                                  <ol>
    !>                                      <li>    type `integer` of kind \IKALL, or
    !>                                      <li>    type `complex` of kind \CKALL, or
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the subset of matrix \f$\ms{matA}\f$ in the triangular matrix update.<br>
    !>  \param[in]      operationA  :   The input scalar `parameter` that can be,
    !>                                  <ol>
    !>                                      <li>    the constant [transSymm](@ref pm_matrixTrans::transSymm) if `matA` is of type `integer`, `complex`, or `real`,
    !>                                              implying the use of the Symmetric form of the specified subset of `matA` in the triangular matrix update.
    !>                                      <li>    the constant [transHerm](@ref pm_matrixTrans::transHerm) if `matA` is of type `complex` and `alpha` and `beta` are of type `real`,
    !>                                              implying the use of the Hermitian transpose form of the specified subset of `matA` in the triangular matrix update.
    !>                                  </ol>
    !>                                  Specifying this argument changes the shape of the subset of `matA` used in the triangular matrix update. See the description of the input argument `ndum`.<br>
    !>                                  This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.<br>
    !>                                  (**optional**. If missing, the specified subset of `matA` will be used as is, without transposition.)
    !>  \param[in]      alpha       :   The input scalar containing the coefficient \f$\alpha\f$ in the triangular matrix update.<br>
    !>                                  <ol>
    !>                                      <li>    If `matA` is of type `integer` or `real`, then `alpha` must of the same type and kind as `matA`.
    !>                                      <li>    If `matA` is of type `complex`, then `alpha` must be of the type `complex` or `real` of the same kind as `matA`.
    !>                                              <ol>
    !>                                                  <li>    If `alpha` is of type `complex`, then the resulting matrix update will be Symmetric of the form,
    !>                                                          \f{align*}{
    !>                                                              & \ms{mat} \leftarrow \alpha \ms{matA}~~ \ms{matA}^T + \beta \ms{mat} &&\text{if operationA is missing.} \\
    !>                                                              & \ms{mat} \leftarrow \alpha \ms{matA}^T \ms{matA}~~ + \beta \ms{mat} &&\text{if operationA = transSymm}
    !>                                                          \f}
    !>                                                  <li>    If `alpha` is of type    `real`, then the resulting matrix update will be Hermitian of the form,
    !>                                                          \f{align*}{
    !>                                                              & \ms{mat} \leftarrow \alpha \ms{matA}~~ \ms{matA}^H + \beta \ms{mat} &&\text{if operationA is missing.} \\
    !>                                                              & \ms{mat} \leftarrow \alpha \ms{matA}^H \ms{matA}~~ + \beta \ms{mat} &&\text{if operationA = transHerm}
    !>                                                          \f}
    !>                                              </ol>
    !>                                  </ol>
    !>  \param[in]      beta        :   The input scalar containing the coefficient \f$\beta\f$ in the triangular matrix update.<br>
    !>                                  <ol>
    !>                                      <li>    If `matA` is of type `integer` or `real`, then `beta` must of the same type and kind as `matA`.
    !>                                      <li>    If `matA` is of type `complex`, then `beta` must be of the type `complex` or `real` of the same kind as `matA`.
    !>                                              <ol>
    !>                                                  <li>    If `beta` is of type `complex`, then the resulting matrix update will be Symmetric of the form,
    !>                                                          \f{align*}{
    !>                                                              & \ms{mat} \leftarrow \beta \ms{matA}~~ \ms{matA}^T + \beta \ms{mat} &&\text{if operationA is missing.} \\
    !>                                                              & \ms{mat} \leftarrow \beta \ms{matA}^T \ms{matA}~~ + \beta \ms{mat} &&\text{if operationA = transSymm}
    !>                                                          \f}
    !>                                                  <li>    If `beta` is of type    `real`, then the resulting matrix update will be Hermitian of the form,
    !>                                                          \f{align*}{
    !>                                                              & \ms{mat} \leftarrow \beta \ms{matA}~~ \ms{matA}^H + \beta \ms{mat} &&\text{if operationA is missing.} \\
    !>                                                              & \ms{mat} \leftarrow \beta \ms{matA}^H \ms{matA}~~ + \beta \ms{mat} &&\text{if operationA = transHerm}
    !>                                                          \f}
    !>                                              </ol>
    !>                                  </ol>
    !>  \param[in]      ndim        :   The input non-negative scalar of type `integer` of default kind \IK,
    !>                                  containing the number of rows/columns (i.e., the rank) of `mat` used in the update `mat(roff + 1 : roff + ndim, coff + 1 : coff + ndim)`.<br>
    !>  \param[in]      ndum        :   The input non-negative scalar of type `integer` of default kind \IK,
    !>                                  containing the number of dummy rows or columns of `matA` used in the triangular matrix update.<br>
    !>                                  <ol>
    !>                                      <li>    For matrix `matA`,
    !>                                              <ol>
    !>                                                  <li>    If the input argument `operationA` is missing, then the subset `matA(roffA + 1 : roffA + ndim, coffA + 1 : coffA + ndum)` is used in the triangular matrix update.
    !>                                                  <li>    If the input argument `operationA` is [transSymm](@ref pm_matrixTrans::transSymm), then the subset `transpose(matA(roffA + 1 : roffA + ndum, coffA + 1 : coffA + ndim)` is used in the triangular matrix update.
    !>                                                  <li>    If the input argument `operationA` is [transHerm](@ref pm_matrixTrans::transHerm), then the subset `conjg(transpose(matA(roffA + 1 : roffA + ndum, coffA + 1 : coffA + ndim))` is used in the triangular matrix update.
    !>                                              </ol>
    !>                                  </ol>
    !>  \param[in]      roff        :   The input non-negative scalar of type `integer` of default kind \IK containing the <b>off</b>set from the first <b>r</b>ow of the input matrix `mat`,
    !>                                  such that `mat(1 + roff, 1 + coff)` marks the top-left corner of the block of `mat` used in the triangular matrix update.<br>
    !>  \param[in]      coff        :   The input non-negative scalar of type `integer` of default kind \IK containing the <b>off</b>set from the first <b>c</b>olumn of the input matrix `mat`,
    !>                                  such that `mat(1 + roff, 1 + coff)` marks the top-left corner of the block of `mat` used in the triangular matrix update.<br>
    !>  \param[in]      roffA       :   The input non-negative scalar of type `integer` of default kind \IK containing the <b>off</b>set from the first <b>r</b>ow of the input matrix `matA`,
    !>                                  such that `matA(1 + roffA, 1 + coffA)` marks the top-left corner of the block of `matA` used in the triangular matrix update.<br>
    !>  \param[in]      coffA       :   The input non-negative scalar of type `integer` of default kind \IK containing the <b>off</b>set from the first <b>c</b>olumn of the input matrix `matA`,
    !>                                  such that `matA(1 + roffA, 1 + coffA)` marks the top-left corner of the block of `matA` used in the triangular matrix update.<br>
    !>
    !>  \interface{setMatUpdate}
    !>  \code{.F90}
    !>
    !>      use pm_matrixUpdate, only: transSymm, transHerm ! Possible values of `operationA`
    !>      use pm_matrixUpdate, only: uppDia, lowDia ! subset: upper-diagonal, lower-diagonal
    !>      use pm_matrixUpdate, only: setMatUpdate
    !>
    !>      call setMatUpdate(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
    !>      call setMatUpdate(matA = transSymm, mat, subset, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
    !>      call setMatUpdate(matA = transHerm, mat, subset, alpha, beta, ndim, ndum, roff, coff, roffA, coffA) ! complex(*) :: matA, real(kind(matA)) :: alpha, beta
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `alpha%im == 0` must hold when `class = hermitian`.<br>
    !>  The condition `all(0 <= [roffA, coffA])` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= [roff, coff])` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0_IK <= [ndim, ndim, ndum])` must hold for the corresponding input arguments.<br>
    !>  The condition `all([roffA + ndim, coffA + ndum] <= shape(matA))` must hold for the corresponding input arguments when the input argument `operationA` is missing.<br>
    !>  The condition `all([roffA + ndum, coffA + ndim] <= shape(matA))` must hold for the corresponding input arguments when the input argument `operationA` is present.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \devnote
    !>  The optional input argument `operationA` in the current interface was originally designed be either missing or set to [trans](@ref pm_matrixTrans::trans) standing for *transposed*.<br>
    !>  While such convention offers a cleaner interface, it was later replaced with the two possible options [transSymm](@ref pm_matrixTrans::transSymm) and [transHerm](@ref pm_matrixTrans::transHerm).<br>
    !>  Although the new interface is uglier and verbose, it allows the possibility of extending this interface in the future without breaking the existing interface.<br>
    !>  The redundancy also offers an automatic compile-time check against semantic bugs, given that `alpha, beta` must be of type `real` when [transHerm](@ref pm_matrixTrans::transHerm) transposition is requested.<br>
    !>
    !>  \lapack{3.11}
    !>  `SSYRK`, `DSYRK`, `CSYRK`, `ZSYRK`, `CHERK`, and `ZHERK`.<br>
    !>  Notably, the interfaces are also extended to support matrices of type `integer` of arbitrary kinds.<br>
    !>
    !>  \lapackint{setMatUpdate}
    !>
    !>  BLAS/LAPACK interface arguments |   [setMatUpdate](@ref pm_matrixUpdate::setMatUpdate) interface arguments
    !>  --------------------------------|-------------------------------------------------------------------------------------------
    !>  `uplo = "U"`                    |   `subset =` [uppDia](@ref pm_matrixSubset::uppDia)
    !>  `uplo = "L"`                    |   `subset =` [lowDia](@ref pm_matrixSubset::lowDia)
    !>  `n`                             |   `ndim`
    !>  `k`                             |   `ndum`
    !>  `trans = "N"`                   |   `operationA` argument is missing.
    !>  `trans = "T"`                   |   `operationA =` [transSymm](@ref pm_matrixTrans::transSymm)
    !>  `trans = "C"`                   |   `operationA =` [transHerm](@ref pm_matrixTrans::transHerm)
    !>  `alpha`                         |   `alpha`
    !>  `a = A(i,j)`                    |   `matA = A, roffA = i - 1, coffA = j - 1`
    !>  `lda`                           |   NONE (passed implicitly).
    !>  `beta`                          |   `beta`
    !>  `c = C(i,j)`                    |   `mat = C, roff = i - 1, coff = j - 1`
    !>  `ldc`                           |   NONE (passed implicitly).
    !>
    !>  \see
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [transSymm](@ref pm_matrixTrans::transSymm)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [setMatUpdateR1](@ref pm_matrixUpdate::setMatUpdateR1)<br>
    !>
    !>  \example{setMatUpdate}
    !>  \include{lineno} example/pm_matrixUpdate/setMatUpdate/main.F90
    !>  \compilef{setMatUpdate}
    !>  \output{setMatUpdate}
    !>  \include{lineno} example/pm_matrixUpdate/setMatUpdate/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixUpdate](@ref test_pm_matrixUpdate)
    !>
    !>  \todo
    !>  \pmed
    !>  The input shape and offset arguments can be made optional by adding new interfaces to the module.<br>
    !>  This should be done only after performing relevant benchmarks with the current interface
    !>  to gauge whether the extension of this module is worth the effort.<br>
    !>
    !>  \finmain{setMatUpdate}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface setMatUpdate

    ! naming convention:
    ! shrk_EXP_LXX
    !                   ||||||
    !                   ||||||
    !                   ||||A/S/H => asis/Symmetric/Hermitian
    !                   |||L/U => LowDia/UppDia
    !                   ||O => offset
    !                   |S => shape
    !                   C => coefficients

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_ONO_IK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_ONO_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_ONO_IK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_ONO_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_ONO_IK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_ONO_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_ONO_IK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_ONO_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_ONO_IK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_ONO_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_ONO_CK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_ONO_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_ONO_CK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_ONO_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_ONO_CK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_ONO_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_ONO_CK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_ONO_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_ONO_CK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_ONO_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_ONO_RK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_ONO_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_ONO_RK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_ONO_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_ONO_RK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_ONO_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_ONO_RK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_ONO_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_ONO_RK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_ONO_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_OTP_IK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_OTP_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_OTP_IK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_OTP_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_OTP_IK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_OTP_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_OTP_IK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_OTP_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_OTP_IK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_OTP_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_OTP_CK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_OTP_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_OTP_CK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_OTP_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_OTP_CK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_OTP_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_OTP_CK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_OTP_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_OTP_CK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_OTP_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_OTP_RK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_OTP_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_OTP_RK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_OTP_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_OTP_RK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_OTP_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_OTP_RK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_OTP_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine shrk_ASS_CSM_SLD_OTP_RK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SLD_OTP_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_ONO_IK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_ONO_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_ONO_IK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_ONO_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_ONO_IK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_ONO_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_ONO_IK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_ONO_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_ONO_IK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_ONO_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_ONO_CK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_ONO_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_ONO_CK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_ONO_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_ONO_CK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_ONO_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_ONO_CK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_ONO_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_ONO_CK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_ONO_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_ONO_RK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_ONO_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_ONO_RK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_ONO_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_ONO_RK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_ONO_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_ONO_RK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_ONO_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_ONO_RK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_ONO_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_OTP_IK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_OTP_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_OTP_IK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_OTP_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_OTP_IK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_OTP_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_OTP_IK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_OTP_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_OTP_IK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_OTP_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_OTP_CK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_OTP_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_OTP_CK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_OTP_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_OTP_CK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_OTP_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_OTP_CK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_OTP_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_OTP_CK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_OTP_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_OTP_RK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_OTP_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_OTP_RK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_OTP_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_OTP_RK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_OTP_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_OTP_RK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_OTP_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine shrk_ASS_CSM_SUD_OTP_RK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CSM_SUD_OTP_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_ONO_IK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_ONO_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_ONO_IK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_ONO_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_ONO_IK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_ONO_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_ONO_IK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_ONO_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_ONO_IK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_ONO_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_ONO_CK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_ONO_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_ONO_CK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_ONO_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_ONO_CK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_ONO_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_ONO_CK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_ONO_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_ONO_CK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_ONO_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_ONO_RK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_ONO_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_ONO_RK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_ONO_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_ONO_RK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_ONO_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_ONO_RK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_ONO_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_ONO_RK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_ONO_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_OTP_IK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_OTP_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_OTP_IK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_OTP_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_OTP_IK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_OTP_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_OTP_IK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_OTP_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_OTP_IK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_OTP_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_OTP_CK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_OTP_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_OTP_CK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_OTP_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_OTP_CK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_OTP_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_OTP_CK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_OTP_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_OTP_CK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_OTP_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_OTP_RK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_OTP_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_OTP_RK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_OTP_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_OTP_RK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_OTP_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_OTP_RK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_OTP_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine shrk_ASS_CHM_SLD_OTP_RK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SLD_OTP_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_ONO_IK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_ONO_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_ONO_IK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_ONO_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_ONO_IK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_ONO_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_ONO_IK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_ONO_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_ONO_IK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_ONO_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_ONO_CK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_ONO_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_ONO_CK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_ONO_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_ONO_CK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_ONO_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_ONO_CK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_ONO_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_ONO_CK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_ONO_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_ONO_RK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_ONO_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_ONO_RK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_ONO_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_ONO_RK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_ONO_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_ONO_RK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_ONO_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_ONO_RK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_ONO_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_OTP_IK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_OTP_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_OTP_IK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_OTP_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_OTP_IK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_OTP_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_OTP_IK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_OTP_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_OTP_IK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_OTP_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_OTP_CK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_OTP_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_OTP_CK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_OTP_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_OTP_CK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_OTP_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_OTP_CK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_OTP_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_OTP_CK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_OTP_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_OTP_RK5(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_OTP_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_OTP_RK4(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_OTP_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_OTP_RK3(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_OTP_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_OTP_RK2(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_OTP_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine shrk_ASS_CHM_SUD_OTP_RK1(mat, class, subset, matA, operationA, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_ASS_CHM_SUD_OTP_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_ONO_IK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_ONO_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_ONO_IK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_ONO_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_ONO_IK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_ONO_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_ONO_IK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_ONO_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_ONO_IK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_ONO_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_ONO_CK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_ONO_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_ONO_CK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_ONO_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_ONO_CK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_ONO_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_ONO_CK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_ONO_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_ONO_CK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_ONO_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_ONO_RK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_ONO_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_ONO_RK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_ONO_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_ONO_RK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_ONO_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_ONO_RK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_ONO_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_ONO_RK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_ONO_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_OTP_IK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_OTP_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_OTP_IK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_OTP_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_OTP_IK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_OTP_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_OTP_IK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_OTP_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_OTP_IK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_OTP_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_OTP_CK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_OTP_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_OTP_CK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_OTP_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_OTP_CK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_OTP_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_OTP_CK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_OTP_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_OTP_CK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_OTP_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_OTP_RK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_OTP_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_OTP_RK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_OTP_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_OTP_RK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_OTP_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_OTP_RK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_OTP_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine shrk_EXP_CSM_SLD_OTP_RK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SLD_OTP_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_ONO_IK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_ONO_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_ONO_IK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_ONO_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_ONO_IK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_ONO_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_ONO_IK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_ONO_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_ONO_IK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_ONO_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_ONO_CK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_ONO_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_ONO_CK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_ONO_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_ONO_CK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_ONO_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_ONO_CK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_ONO_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_ONO_CK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_ONO_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_ONO_RK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_ONO_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_ONO_RK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_ONO_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_ONO_RK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_ONO_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_ONO_RK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_ONO_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_ONO_RK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_ONO_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_OTP_IK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_OTP_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_OTP_IK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_OTP_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_OTP_IK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_OTP_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_OTP_IK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_OTP_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_OTP_IK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_OTP_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_OTP_CK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_OTP_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_OTP_CK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_OTP_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_OTP_CK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_OTP_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_OTP_CK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_OTP_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_OTP_CK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_OTP_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_OTP_RK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_OTP_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_OTP_RK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_OTP_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_OTP_RK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_OTP_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_OTP_RK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_OTP_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine shrk_EXP_CSM_SUD_OTP_RK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CSM_SUD_OTP_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(symmetric_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_ONO_IK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_ONO_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_ONO_IK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_ONO_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_ONO_IK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_ONO_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_ONO_IK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_ONO_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_ONO_IK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_ONO_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_ONO_CK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_ONO_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_ONO_CK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_ONO_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_ONO_CK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_ONO_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_ONO_CK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_ONO_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_ONO_CK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_ONO_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_ONO_RK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_ONO_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_ONO_RK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_ONO_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_ONO_RK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_ONO_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_ONO_RK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_ONO_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_ONO_RK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_ONO_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_OTP_IK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_OTP_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_OTP_IK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_OTP_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_OTP_IK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_OTP_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_OTP_IK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_OTP_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_OTP_IK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_OTP_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_OTP_CK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_OTP_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_OTP_CK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_OTP_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_OTP_CK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_OTP_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_OTP_CK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_OTP_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_OTP_CK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_OTP_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_OTP_RK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_OTP_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_OTP_RK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_OTP_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_OTP_RK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_OTP_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_OTP_RK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_OTP_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine shrk_EXP_CHM_SLD_OTP_RK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SLD_OTP_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(lowDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_ONO_IK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_ONO_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_ONO_IK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_ONO_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_ONO_IK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_ONO_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_ONO_IK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_ONO_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_ONO_IK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_ONO_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_ONO_CK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_ONO_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_ONO_CK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_ONO_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_ONO_CK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_ONO_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_ONO_CK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_ONO_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_ONO_CK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_ONO_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_ONO_RK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_ONO_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_ONO_RK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_ONO_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_ONO_RK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_ONO_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_ONO_RK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_ONO_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_ONO_RK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_ONO_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(nothing_type)      , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_OTP_IK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_OTP_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_OTP_IK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_OTP_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_OTP_IK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_OTP_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_OTP_IK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_OTP_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_OTP_IK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_OTP_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_OTP_CK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_OTP_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_OTP_CK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_OTP_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_OTP_CK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_OTP_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_OTP_CK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_OTP_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_OTP_CK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_OTP_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_OTP_RK5(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_OTP_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_OTP_RK4(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_OTP_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_OTP_RK3(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_OTP_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_OTP_RK2(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_OTP_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine shrk_EXP_CHM_SUD_OTP_RK1(mat, class, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: shrk_EXP_CHM_SUD_OTP_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)                    :: alpha, beta
        type(trans_type)        , intent(in)                    :: operationA
        type(uppDia_type)       , intent(in)                    :: subset
        type(hermitian_type)    , intent(in)                    :: class
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the rank-1 update of the input matrix `mat` using the input vectors `vecA` and `vecB`.
    !>
    !>  \details
    !>  see the documentation of [pm_matrixUpdate](@ref pm_matrixUpdate) for more details.<br>
    !>
    !>  \param[inout]   mat         :   The input/output `contiguous` matrix of arbitrary shape `(:,:)` of,<br>
    !>                                  <ol>
    !>                                      <li>    type `integer` of kind \IKALL,<br>
    !>                                      <li>    type `complex` of kind \CKALL,<br>
    !>                                      <li>    type `real` of kind \RKALL,<br>
    !>                                  </ol>
    !>                                  containing the matrix to be updated on output.<br>
    !>                                  The matrix does not have to be square.<br>
    !>  \param[in]      vecA        :   The input `contiguous` vector of the same type and kind as `mat`, of either
    !>                                  <ol>
    !>                                      <li>    the same size as `size(mat, 1)`, or<br>
    !>                                      <li>    the specific size such that the condition `(size(vecA, 1) - 1) / abs(incA) + 1 + roff == size(mat, 1)`.<br>
    !>                                  </ol>
    !>  \param[in]      vecB        :   The input `contiguous` vector of the same type and kind as `mat`, of either
    !>                                  <ol>
    !>                                      <li>    the same size as `size(mat, 2)`, or<br>
    !>                                      <li>    the specific size such that the condition `(size(vecB, 1) - 1) / abs(incB) + 1 + roff == size(mat, 2)`.<br>
    !>                                  </ol>
    !>  \param[in]      operationB  :   The input scalar `parameter` object [transHerm](@ref pm_matrixTrans::transHerm) of type [transHerm_type](@ref pm_matrixTrans::transHerm_type).<br>
    !>                                  This argument is merely used to require complex conjugate (Hermitian) transpose in the computation of the outer product: \f$\ms{vecA}~\ms{vecB}^{\mathrm{H}}\f$.<br>
    !>                                  This argument is currently available for and relevant only to input arguments `mat` of type `complex`.<br>
    !>                                  (**optional**. If missing, the real transpose (without conjugation) will be used instead.)
    !>  \param[in]      alpha       :   The input scalar of the same type and kind as `mat`,
    !>                                  representing the factor by which the outer product \f$\ms{vecA}~\ms{vecB}^{\mathrm{T}}\f$ must be multiplied.<br>
    !>                                  (**optional**, default = `1.`)
    !>  \param[in]      incA        :   The input scalar of type `integer` of default kind \IK, representing the increment along the `vecA` vector.<br>
    !>                                  It can be any value other than `0`. An input value `incA < 0` is equivalent to inverting `vecA` and resetting `incA = abs(incA)`.<br>
    !>  \param[in]      incA        :   The input scalar of type `integer` of default kind \IK, representing the increment along the `vecB` vector.<br>
    !>                                  It can be any value other than `0`. An input value `incB < 0` is equivalent to inverting `vecB` and resetting `incB = abs(incB)`.<br>
    !>  \param[in]      roff        :   The input scalar of type `integer` of default kind \IK, representing the offset from the starting row of `mat` to which the update must occur.<br>
    !>
    !>  \interface{setMatUpdateR1}
    !>  \code{.F90}
    !>
    !>      use pm_matrixUpdate, only: setMatUpdateR1
    !>
    !>      call setMatUpdateR1(mat, vecA, vecB, incA, incB, offset)
    !>      call setMatUpdateR1(mat, vecA, vecB, alpha, incA, incB, offset)
    !>      call setMatUpdateR1(mat, vecA, vecB, operationB, incA, incB, offset) ! only complex arguments (cgerc, zgerc).
    !>      call setMatUpdateR1(mat, vecA, vecB, operationB, alpha, incA, incB, offset) ! only complex arguments (cgerc, zgerc).
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `incA /= 0_IK` must hold for the corresponding input arguments.<br>
    !>  The condition `incB /= 0_IK` must hold for the corresponding input arguments.<br>
    !>  The condition `0_IK <= offset` must hold for the corresponding input arguments.<br>
    !>  The condition `size(mat, 2, IK) == (size(vecB, 1, IK) - 1_IK) / abs(incB) + 1_IK` must hold for the corresponding input arguments.<br>
    !>  The condition `size(mat, 1, IK) >= (size(vecA, 1, IK) - 1_IK) / abs(incA) + 1_IK + offset` must hold for the corresponding input arguments.<br>
    !>  The condition `size(vecA, 1, IK) == incA * ((size(vecA, 1, IK) - 1_IK) / abs(incA) + 1_IK)` must hold for the corresponding input arguments.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \lapack{3.11}
    !>  `SGER`, `DGER`, `CGERU`, `ZGERU`, `CGERC`, and `ZGERC`.<br>
    !>  <ol>
    !>      <li>    The interface corresponds to `CGERC` and `ZGERC` when the input optional argument `operationA = transHerm_type()` is present.<br>
    !>      <li>    Otherwise, the interface corresponds to `SGER`, `DGER`, `CGERU`, `ZGERU`.<br>
    !>  </ol>
    !>
    !>  \see
    !>  [netlib::LAPACK](http://netlib.org/lapack/)<br>
    !>  The IBM Engineering and Scientific Subroutine Library.<br>
    !>  Developer Reference for Intel® oneAPI Math Kernel Library - Fortran.<br>
    !>  Dongarra, J. J.; DuCroz, J.; Hammarling, S.; Hanson, R. J. March 1988. *An Extended Set of Fortran Basic Linear Algebra Subprograms.* ACM Transactions on Mathematical Software , 14(1):1–17.<br>
    !>
    !>  \example{setMatUpdateR1}
    !>  \include{lineno} example/pm_matrixUpdate/setMatUpdateR1/main.F90
    !>  \compilef{setMatUpdateR1}
    !>  \output{setMatUpdateR1}
    !>  \include{lineno} example/pm_matrixUpdate/setMatUpdateR1/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixUpdate](@ref test_pm_matrixUpdate)
    !>
    !>  \todo
    !>  \pmed
    !>  Benchmarks comparing this interface with LAPACK routines and conventional approach should be added to the documentation.
    !>
    !>  \finmain{setMatUpdateR1}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 11:34 PM, August 10, 2021, Dallas, TX
    interface setMatUpdateR1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatUpdateR1F_IK5(mat, vecA, vecB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1F_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: incA, incB, roff
        integer(IKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatUpdateR1F_IK4(mat, vecA, vecB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1F_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: incA, incB, roff
        integer(IKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatUpdateR1F_IK3(mat, vecA, vecB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1F_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: incA, incB, roff
        integer(IKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatUpdateR1F_IK2(mat, vecA, vecB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1F_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: incA, incB, roff
        integer(IKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatUpdateR1F_IK1(mat, vecA, vecB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1F_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: incA, incB, roff
        integer(IKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatUpdateR1H_CK5(mat, vecA, vecB, operationB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1H_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatUpdateR1H_CK4(mat, vecA, vecB, operationB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1H_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatUpdateR1H_CK3(mat, vecA, vecB, operationB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1H_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatUpdateR1H_CK2(mat, vecA, vecB, operationB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1H_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatUpdateR1H_CK1(mat, vecA, vecB, operationB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1H_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatUpdateR1F_CK5(mat, vecA, vecB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1F_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatUpdateR1F_CK4(mat, vecA, vecB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1F_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatUpdateR1F_CK3(mat, vecA, vecB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1F_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatUpdateR1F_CK2(mat, vecA, vecB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1F_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatUpdateR1F_CK1(mat, vecA, vecB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1F_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatUpdateR1F_RK5(mat, vecA, vecB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1F_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: incA, incB, roff
        real(RKC)               , intent(in)    , contiguous    :: vecA(:), vecB(:)
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatUpdateR1F_RK4(mat, vecA, vecB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1F_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: incA, incB, roff
        real(RKC)               , intent(in)    , contiguous    :: vecA(:), vecB(:)
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatUpdateR1F_RK3(mat, vecA, vecB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1F_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: incA, incB, roff
        real(RKC)               , intent(in)    , contiguous    :: vecA(:), vecB(:)
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatUpdateR1F_RK2(mat, vecA, vecB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1F_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: incA, incB, roff
        real(RKC)               , intent(in)    , contiguous    :: vecA(:), vecB(:)
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatUpdateR1F_RK1(mat, vecA, vecB, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1F_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: incA, incB, roff
        real(RKC)               , intent(in)    , contiguous    :: vecA(:), vecB(:)
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatUpdateR1A_IK5(mat, vecA, vecB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1A_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: incA, incB, roff
        integer(IKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)                    :: alpha
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatUpdateR1A_IK4(mat, vecA, vecB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1A_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: incA, incB, roff
        integer(IKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)                    :: alpha
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatUpdateR1A_IK3(mat, vecA, vecB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1A_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: incA, incB, roff
        integer(IKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)                    :: alpha
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatUpdateR1A_IK2(mat, vecA, vecB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1A_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: incA, incB, roff
        integer(IKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)                    :: alpha
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatUpdateR1A_IK1(mat, vecA, vecB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1A_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: incA, incB, roff
        integer(IKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        integer(IKC)            , intent(inout) , contiguous    :: mat(:,:)
        integer(IKC)            , intent(in)                    :: alpha
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatUpdateR1AH_CK5(mat, vecA, vecB, operationB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1AH_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)                    :: alpha
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatUpdateR1AH_CK4(mat, vecA, vecB, operationB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1AH_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)                    :: alpha
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatUpdateR1AH_CK3(mat, vecA, vecB, operationB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1AH_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)                    :: alpha
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatUpdateR1AH_CK2(mat, vecA, vecB, operationB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1AH_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)                    :: alpha
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatUpdateR1AH_CK1(mat, vecA, vecB, operationB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1AH_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)                    :: alpha
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatUpdateR1A_CK5(mat, vecA, vecB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1A_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)                    :: alpha
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatUpdateR1A_CK4(mat, vecA, vecB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1A_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)                    :: alpha
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatUpdateR1A_CK3(mat, vecA, vecB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1A_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)                    :: alpha
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatUpdateR1A_CK2(mat, vecA, vecB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1A_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)                    :: alpha
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatUpdateR1A_CK1(mat, vecA, vecB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1A_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: incA, incB, roff
        complex(CKC)            , intent(in)    , contiguous    :: vecA(:), vecB(:)
        complex(CKC)            , intent(inout) , contiguous    :: mat(:,:)
        complex(CKC)            , intent(in)                    :: alpha
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatUpdateR1A_RK5(mat, vecA, vecB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1A_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: incA, incB, roff
        real(RKC)               , intent(in)    , contiguous    :: vecA(:), vecB(:)
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)                    :: alpha
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatUpdateR1A_RK4(mat, vecA, vecB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1A_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: incA, incB, roff
        real(RKC)               , intent(in)    , contiguous    :: vecA(:), vecB(:)
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)                    :: alpha
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatUpdateR1A_RK3(mat, vecA, vecB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1A_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: incA, incB, roff
        real(RKC)               , intent(in)    , contiguous    :: vecA(:), vecB(:)
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)                    :: alpha
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatUpdateR1A_RK2(mat, vecA, vecB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1A_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: incA, incB, roff
        real(RKC)               , intent(in)    , contiguous    :: vecA(:), vecB(:)
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)                    :: alpha
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatUpdateR1A_RK1(mat, vecA, vecB, alpha, incA, incB, roff)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateR1A_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: incA, incB, roff
        real(RKC)               , intent(in)    , contiguous    :: vecA(:), vecB(:)
        real(RKC)               , intent(inout) , contiguous    :: mat(:,:)
        real(RKC)               , intent(in)                    :: alpha
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the result of arbitrary-rank Symmetric/Hermitian updates
    !>  to triangular matrices of type `integer`, `complex`, and `real` of arbitrary type-kind parameters.
    !>
    !>  \details
    !>  see the documentation of [pm_matrixUpdate](@ref pm_matrixUpdate) for more details.<br>
    !>
    !>  \param[inout]   mat         :   The input/output `contiguous` matrix of arbitrary shape `(:,:)` of the same type and kind as `matA`,
    !>                                  containing the subset of matrix \f$\texttt{mat}\f$ in the triangular matrix update.<br>
    !>  \param[in]      subset      :   The input scalar that can be either,
    !>                                  <ol>
    !>                                      <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) implying that
    !>                                              the upper-diagonal triangular block of `mat` should be used for storage and updating while the lower triangle remains intact.
    !>                                      <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) implying that
    !>                                              the lower-diagonal triangular block of `mat` should be used for storage and updating while the upper triangle remains intact.
    !>                                  </ol>
    !>                                  This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.
    !>  \param[in]      matA        :   The input `contiguous` matrix of arbitrary shape `(:,:)` of,
    !>                                  <ol>
    !>                                      <li>    type `integer` of kind \IKALL, or
    !>                                      <li>    type `complex` of kind \CKALL, or
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the subset of matrix \f$\texttt{matA}\f$ in the triangular matrix update.<br>
    !>  \param[in]      operationA  :   The input scalar `parameter` that can be,
    !>                                  <ol>
    !>                                      <li>    the constant [transSymm](@ref pm_matrixTrans::transSymm) if `matA` is of type `integer`, `complex`, or `real`,
    !>                                              implying the use of the Symmetric form of the specified subset of `matA` in the triangular matrix update.
    !>                                      <li>    the constant [transHerm](@ref pm_matrixTrans::transHerm) if `matA` is of type `complex` and `alpha` and `beta` are of type `real`,
    !>                                              implying the use of the Hermitian transpose form of the specified subset of `matA` in the triangular matrix update.
    !>                                  </ol>
    !>                                  Specifying this argument changes the shape of the subset of `matA` used in the triangular matrix update. See the description of the input argument `ndum`.<br>
    !>                                  This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.<br>
    !>                                  (**optional**. If missing, the specified subset of `matA` will be used as is, without transposition.)
    !>  \param[in]      alpha       :   The input scalar containing the coefficient \f$\alpha\f$ in the triangular matrix update.<br>
    !>                                  <ol>
    !>                                      <li>    If `matA` is of type `integer` or `real`, then `alpha` must of the same type and kind as `matA`.
    !>                                      <li>    If `matA` is of type `complex`, then `alpha` must be of the type `complex` or `real` of the same kind as `matA`.
    !>                                              <ol>
    !>                                                  <li>    If `alpha` is of type `complex`, then the resulting matrix update will be Symmetric of the form,
    !>                                                          \f{align*}{
    !>                                                              & \texttt{mat} \leftarrow \alpha \texttt{matA}~~ \texttt{matA}^T + \beta \texttt{mat} &&\text{if operationA is missing.} \\
    !>                                                              & \texttt{mat} \leftarrow \alpha \texttt{matA}^T \texttt{matA}~~ + \beta \texttt{mat} &&\text{if operationA = transSymm}
    !>                                                          \f}
    !>                                                  <li>    If `alpha` is of type    `real`, then the resulting matrix update will be Hermitian of the form,
    !>                                                          \f{align*}{
    !>                                                              & \texttt{mat} \leftarrow \alpha \texttt{matA}~~ \texttt{matA}^H + \beta \texttt{mat} &&\text{if operationA is missing.} \\
    !>                                                              & \texttt{mat} \leftarrow \alpha \texttt{matA}^H \texttt{matA}~~ + \beta \texttt{mat} &&\text{if operationA = transHerm}
    !>                                                          \f}
    !>                                              </ol>
    !>                                  </ol>
    !>  \param[in]      beta        :   The input scalar containing the coefficient \f$\beta\f$ in the triangular matrix update.<br>
    !>                                  <ol>
    !>                                      <li>    If `matA` is of type `integer` or `real`, then `beta` must of the same type and kind as `matA`.
    !>                                      <li>    If `matA` is of type `complex`, then `beta` must be of the type `complex` or `real` of the same kind as `matA`.
    !>                                              <ol>
    !>                                                  <li>    If `beta` is of type `complex`, then the resulting matrix update will be Symmetric of the form,
    !>                                                          \f{align*}{
    !>                                                              & \texttt{mat} \leftarrow \beta \texttt{matA}~~ \texttt{matA}^T + \beta \texttt{mat} &&\text{if operationA is missing.} \\
    !>                                                              & \texttt{mat} \leftarrow \beta \texttt{matA}^T \texttt{matA}~~ + \beta \texttt{mat} &&\text{if operationA = transSymm}
    !>                                                          \f}
    !>                                                  <li>    If `beta` is of type    `real`, then the resulting matrix update will be Hermitian of the form,
    !>                                                          \f{align*}{
    !>                                                              & \texttt{mat} \leftarrow \beta \texttt{matA}~~ \texttt{matA}^H + \beta \texttt{mat} &&\text{if operationA is missing.} \\
    !>                                                              & \texttt{mat} \leftarrow \beta \texttt{matA}^H \texttt{matA}~~ + \beta \texttt{mat} &&\text{if operationA = transHerm}
    !>                                                          \f}
    !>                                              </ol>
    !>                                  </ol>
    !>  \param[in]      ndim        :   The input non-negative scalar of type `integer` of default kind \IK,
    !>                                  containing the number of rows/columns (i.e., the rank) of `mat` used in the update `mat(roff + 1 : roff + ndim, coff + 1 : coff + ndim)`.<br>
    !>  \param[in]      ndum        :   The input non-negative scalar of type `integer` of default kind \IK,
    !>                                  containing the number of dummy rows or columns of `matA` used in the triangular matrix update.<br>
    !>                                  <ol>
    !>                                      <li>    For matrix `matA`,
    !>                                              <ol>
    !>                                                  <li>    If the input argument `operationA` is missing, then the subset `matA(roffA + 1 : roffA + ndim, coffA + 1 : coffA + ndum)` is used in the triangular matrix update.
    !>                                                  <li>    If the input argument `operationA` is [transSymm](@ref pm_matrixTrans::transSymm), then the subset `transpose(matA(roffA + 1 : roffA + ndum, coffA + 1 : coffA + ndim)` is used in the triangular matrix update.
    !>                                                  <li>    If the input argument `operationA` is [transHerm](@ref pm_matrixTrans::transHerm), then the subset `conjg(transpose(matA(roffA + 1 : roffA + ndum, coffA + 1 : coffA + ndim))` is used in the triangular matrix update.
    !>                                              </ol>
    !>                                  </ol>
    !>  \param[in]      roffA       :   The input non-negative scalar of type `integer` of default kind \IK containing the <b>off</b>set from the first <b>r</b>ow of the input matrix `matA`,
    !>                                  such that `matA(1 + roffA, 1 + coffA)` marks the top-left corner of the block of `matA` used in the triangular matrix update.<br>
    !>  \param[in]      coffA       :   The input non-negative scalar of type `integer` of default kind \IK containing the <b>off</b>set from the first <b>c</b>olumn of the input matrix `matA`,
    !>                                  such that `matA(1 + roffA, 1 + coffA)` marks the top-left corner of the block of `matA` used in the triangular matrix update.<br>
    !>  \param[in]      roff       :   The input non-negative scalar of type `integer` of default kind \IK containing the <b>off</b>set from the first <b>r</b>ow of the input matrix `mat`,
    !>                                  such that `mat(1 + roff, 1 + coff)` marks the top-left corner of the block of `mat` used in the triangular matrix update.<br>
    !>  \param[in]      coff       :   The input non-negative scalar of type `integer` of default kind \IK containing the <b>off</b>set from the first <b>c</b>olumn of the input matrix `mat`,
    !>                                  such that `mat(1 + roff, 1 + coff)` marks the top-left corner of the block of `mat` used in the triangular matrix update.<br>
    !>
    !>  \interface{setMatUpdateTriang}
    !>  \code{.F90}
    !>
    !>      use pm_matrixUpdate, only: transSymm, transHerm ! Possible values of `operationA`
    !>      use pm_matrixUpdate, only: uppDia, lowDia ! subset: upper-diagonal, lower-diagonal
    !>      use pm_matrixUpdate, only: setMatUpdateTriang
    !>
    !>      call setMatUpdateTriang(matA, mat, subset, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
    !>      call setMatUpdateTriang(matA, operationA, mat, subset, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
    !>      call setMatUpdateTriang(matA, operationA, mat, subset, alpha, beta, ndim, ndum, roff, coff, roffA, coffA) ! complex(*) :: matA, real(kind(matA)) :: alpha, beta
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `all(0 <= [roffA, coffA])` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= [roff, coff])` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0_IK <= [ndim, ndim, ndum])` must hold for the corresponding input arguments.<br>
    !>  The condition `all([roffA + ndim, coffA + ndum] <= shape(matA))` must hold for the corresponding input arguments when the input argument `operationA` is missing.<br>
    !>  The condition `all([roffA + ndum, coffA + ndim] <= shape(matA))` must hold for the corresponding input arguments when the input argument `operationA` is present.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \devnote
    !>  The optional input argument `operationA` in the current interface was originally designed be either missing or set to [trans](@ref pm_matrixTrans::trans) standing for *transposed*.<br>
    !>  While such convention offers a cleaner interface, it was later replaced with the two possible options [transSymm](@ref pm_matrixTrans::transSymm) and [transHerm](@ref pm_matrixTrans::transHerm).<br>
    !>  Although the new interface is uglier and verbose, it allows the possibility of extending this interface in the future without breaking the existing interface.<br>
    !>  The redundancy also offers an automatic compile-time check against semantic bugs, given that `alpha, beta` must be of type `real` when [transHerm](@ref pm_matrixTrans::transHerm) transposition is requested.<br>
    !>
    !>  \lapack{3.11}
    !>  `SSYRK`, `DSYRK`, `CSYRK`, `ZSYRK`, `CHERK`, and `ZHERK`.<br>
    !>  Notably, the interfaces are also extended to support matrices of type `integer` of arbitrary kinds.<br>
    !>
    !>  \lapackint{setMatUpdateTriang}
    !>
    !>  BLAS/LAPACK interface arguments |   [setMatUpdateTriang](@ref pm_matrixUpdate::setMatUpdateTriang) interface arguments
    !>  --------------------------------|-------------------------------------------------------------------------------------------
    !>  `uplo = "U"`                    |   `subset =` [uppDia](@ref pm_matrixSubset::uppDia)
    !>  `uplo = "L"`                    |   `subset =` [lowDia](@ref pm_matrixSubset::lowDia)
    !>  `n`                             |   `ndim`
    !>  `k`                             |   `ndum`
    !>  `trans = "N"`                   |   `operationA` argument is missing.
    !>  `trans = "T"`                   |   `operationA =` [transSymm](@ref pm_matrixTrans::transSymm)
    !>  `trans = "C"`                   |   `operationA =` [transHerm](@ref pm_matrixTrans::transHerm)
    !>  `alpha`                         |   `alpha`
    !>  `a = A(i,j)`                    |   `matA = A, roffA = i - 1, coffA = j - 1`
    !>  `lda`                           |   NONE (passed implicitly).
    !>  `beta`                          |   `beta`
    !>  `c = C(i,j)`                    |   `mat = C, roff = i - 1, coff = j - 1`
    !>  `ldc`                           |   NONE (passed implicitly).
    !>
    !>  \see
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [transSymm](@ref pm_matrixTrans::transSymm)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [setMatUpdateR1](@ref pm_matrixUpdate::setMatUpdateR1)<br>
    !>
    !>  \example{setMatUpdateTriang}
    !>  \include{lineno} example/pm_matrixUpdate/setMatUpdateTriang/main.F90
    !>  \compilef{setMatUpdateTriang}
    !>  \output{setMatUpdateTriang}
    !>  \include{lineno} example/pm_matrixUpdate/setMatUpdateTriang/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixUpdate](@ref test_pm_matrixUpdate)
    !>
    !>  \todo
    !>  \pmed
    !>  The input shape and offset arguments can be made optional by adding new interfaces to the module.<br>
    !>  This should be done only after performing relevant benchmarks with the current interface
    !>  to gauge whether the extension of this module is worth the effort.<br>
    !>
    !>  \finmain{setMatUpdateTriang}
    !>
    !>  \author
    !>  \FatemehBagheri, Tuesday 08:49 PM, August 10, 2021, Dallas, TX
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface setMatUpdateTriang

    ! naming convention:
    ! setMatUpdateTriangCSOLXX
    !                   ||||||
    !                   ||||||
    !                   ||||A/S/H => asis/Symmetric/Hermitian
    !                   |||L/U => LowDia/UppDia
    !                   ||O => offset
    !                   |S => shape
    !                   C => coefficients

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAS_IK5(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAS_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAS_IK4(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAS_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAS_IK3(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAS_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAS_IK2(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAS_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAS_IK1(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAS_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAS_CK5(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAS_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAS_CK4(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAS_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAS_CK3(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAS_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAS_CK2(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAS_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAS_CK1(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAS_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAS_RK5(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAS_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAS_RK4(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAS_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAS_RK3(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAS_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAS_RK2(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAS_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAS_RK1(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAS_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLSA_IK5(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLSA_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLSA_IK4(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLSA_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLSA_IK3(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLSA_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLSA_IK2(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLSA_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLSA_IK1(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLSA_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLSA_CK5(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLSA_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLSA_CK4(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLSA_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLSA_CK3(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLSA_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLSA_CK2(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLSA_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLSA_CK1(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLSA_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLSA_RK5(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLSA_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLSA_RK4(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLSA_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLSA_RK3(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLSA_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLSA_RK2(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLSA_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLSA_RK1(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLSA_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAS_IK5(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAS_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAS_IK4(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAS_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAS_IK3(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAS_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAS_IK2(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAS_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAS_IK1(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAS_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAS_CK5(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAS_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAS_CK4(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAS_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAS_CK3(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAS_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAS_CK2(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAS_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAS_CK1(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAS_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAS_RK5(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAS_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAS_RK4(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAS_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAS_RK3(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAS_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAS_RK2(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAS_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAS_RK1(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAS_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUSA_IK5(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUSA_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUSA_IK4(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUSA_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUSA_IK3(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUSA_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUSA_IK2(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUSA_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUSA_IK1(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUSA_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUSA_CK5(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUSA_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUSA_CK4(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUSA_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUSA_CK3(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUSA_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUSA_CK2(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUSA_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUSA_CK1(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUSA_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUSA_RK5(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUSA_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUSA_RK4(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUSA_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUSA_RK3(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUSA_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUSA_RK2(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUSA_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUSA_RK1(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUSA_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAH_CK5(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAH_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAH_CK4(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAH_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAH_CK3(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAH_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAH_CK2(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAH_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLAH_CK1(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLAH_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLHA_CK5(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLHA_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLHA_CK4(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLHA_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLHA_CK3(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLHA_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLHA_CK2(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLHA_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatUpdateTriangCSOLHA_CK1(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOLHA_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(lowDia_type)       , intent(in)                    :: subset
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAH_CK5(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAH_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAH_CK4(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAH_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAH_CK3(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAH_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAH_CK2(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAH_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUAH_CK1(mat, subset, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUAH_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUHA_CK5(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUHA_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUHA_CK4(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUHA_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUHA_CK3(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUHA_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUHA_CK2(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUHA_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatUpdateTriangCSOUHA_CK1(mat, subset, matA, operationA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatUpdateTriangCSOUHA_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, ndum, roff, coff, roffA, coffA
        real(CKC)               , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(inout) , contiguous    :: mat(1 - roff :, 1 - coff :)
        type(uppDia_type)       , intent(in)                    :: subset
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_matrixUpdate ! LCOV_EXCL_LINE