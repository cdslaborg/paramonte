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
!>  This module contains procedures and generic interfaces relevant to combined matrix-matrix or matrix-vector multiplication and addition.
!>
!>  \details
!>  The procedures under the generic interface [setMatMulAdd](@ref pm_matrixMulAdd::setMatMulAdd) of this module
!>  return the result of the multiplication of the input matrices `matA` and `matB` in one of the following forms,
!>  \f{align*}{
!>      & \ms{matC} \leftarrow \alpha \ms{matA}~~ \ms{matB} + \beta \ms{matC} && \ms{matC} \leftarrow \alpha \ms{matA}~~ \ms{matB}^T + \beta \ms{matC} && \ms{matC} \leftarrow \alpha \ms{matA}~~ \ms{matB}^H + \beta \ms{matC} \\
!>      & \ms{matC} \leftarrow \alpha \ms{matA}^T \ms{matB} + \beta \ms{matC} && \ms{matC} \leftarrow \alpha \ms{matA}^T \ms{matB}^T + \beta \ms{matC} && \ms{matC} \leftarrow \alpha \ms{matA}^T \ms{matB}^H + \beta \ms{matC} \\
!>      & \ms{matC} \leftarrow \alpha \ms{matA}^H \ms{matB} + \beta \ms{matC} && \ms{matC} \leftarrow \alpha \ms{matA}^H \ms{matB}^T + \beta \ms{matC} && \ms{matC} \leftarrow \alpha \ms{matA}^H \ms{matB}^H + \beta \ms{matC}
!>  \f}
!>  where \f$\cdot^T\f$ represents a Symmetric transpose, \f$\cdot^H\f$ represents a Hermitian transpose,
!>  and `matA` or `matB` can be also specified as Symmetric/Hermitian upper/lower triangular matrices.<br>
!>
!>  The following figure illustrates the form of the **general matrix-matrix or matrix-vector multiplication** depending on the input values.<br>
!>  <ol>
!>      <li>    Default Multiplication (no transposition involved).
!>              \image html pm_matrixMulAdd@TNA_TNB.png width=1000
!>              <br>
!>      <li>    Same as the default but with `operationA` set to [transSymm](@ref pm_matrixTrans::transSymm) or [transHerm](@ref pm_matrixTrans::transHerm).<br>
!>              \image html pm_matrixMulAdd@TSA_TNB.png width=1000
!>              <br>
!>      <li>    Same as the default but with `operationB` set to [transSymm](@ref pm_matrixTrans::transSymm) or [transHerm](@ref pm_matrixTrans::transHerm).<br>
!>              \image html pm_matrixMulAdd@TNA_TSB.png width=1000
!>  </ol>
!>  For symmetric or Hermitian matrices, only the upper or lower triangles of the corresponding matrices are required and referenced.<br>
!>  Although the implementation is custom-defined for Symmetric/Hermitian matrices, the multiplication is defined as in the figure of case 1 above.<br>
!>
!>  \note
!>  For triangular matrix-matrix or matrix-vector multiplications, see [pm_matrixMulTri](@ref pm_matrixMulTri).<br>
!>
!>  \lapack{3.11}
!>  `SAXPY`, `DAXPY`, `CAXPY`, `ZAXPY`
!>  `SGEMV`, `DGEMV`, `CGEMV`, `ZGEMV`
!>  `SSPMV`, `DSPMV`, `CHPMV`, `ZHPMV`, `SSYMV`, `DSYMV`, `CHEMV`, `ZHEMV`,
!>  `SGEMM`, `DGEMM`, `CGEMM`, `ZGEMM`, `SSYMM`, `DSYMM`, `CSYMM`, `ZSYMM`, `CHEMM`, `ZHEMM`.<br>
!>  In particular multiplications of matrices of type `integer` of arbitrary kinds are also possible.<br>
!>
!>  \see
!>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
!>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
!>  [symmetric](@ref pm_matrixClass::symmetric)<br>
!>  [hermitian](@ref pm_matrixClass::hermitian)<br>
!>  [transSymm](@ref pm_matrixTrans::transSymm)<br>
!>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
!>  [setMatMulTri](@ref pm_matrixMulTri::setMatMulTri)<br>
!>
!>  \benchmarks
!>
!>  \benchmark{setMatMulAdd, The runtime performance of [setMatMulAdd](@ref pm_matrixMulAdd::setMatMulAdd) vs. other other approaches.}
!>  \include{lineno} benchmark/pm_matrixMulAdd/setMatMulAdd/main.F90
!>  \compilefb{pm_matrixMulAdd/setMatMulAdd}
!>  \postprocb{setMatMulAdd}
!>  \include{lineno} benchmark/pm_matrixMulAdd/setMatMulAdd/main.py
!>  \visb{setMatMulAdd}
!>  \image html benchmark/pm_matrixMulAdd/setMatMulAdd/benchmark.setMatMulAdd.runtime.png width=1000
!>  \image html benchmark/pm_matrixMulAdd/setMatMulAdd/benchmark.setMatMulAdd.runtime.ratio.png width=1000
!>  \moralb{setMatMulAdd}
!>      -#  The procedures under the generic interface [setMatMulAdd](@ref pm_matrixMulAdd::setMatMulAdd) are enhancements to the original BLAS routines.<br>
!>          However, they are currently not optimized for cache-efficient data access on various hardware.<br>
!>          Similarly, [setMatMulAdd](@ref pm_matrixMulAdd::setMatMulAdd) does not utilize blocking methods to improve cache access.<br>
!>      -#  For practical routine daily usages, [setMatMulAdd](@ref pm_matrixMulAdd::setMatMulAdd) offers a nice generic matrix multiplication interface.<br>
!>          However, the current implementations are suboptimal to hardware-tuned BLAS libraries such as OpenBLAS and MKL.<br>
!>          Such performance differences will be most striking for large matrices where cache-efficiency dominates the performance bottleneck.<br>
!>
!>  \test
!>  [test_pm_matrixMulAdd](@ref test_pm_matrixMulAdd)<br>
!>
!>  \todo
!>  \pvhigh
!>  The following BLAS Band-matrix routines must be added to this module:
!>  <ol>
!>      <li>    `SGBMV`, `DGBMV`, `CGBMV`, and `ZGBMV` (Matrix-Vector Product for a General Band Matrix, Its Transpose, or Its Conjugate Transpose).<br>
!>      <li>    `SSBMV`, `DSBMV`, `CHBMV`, and `ZHBMV` (Matrix-Vector Product for a Real Symmetric or Complex Hermitian Band Matrix).<br>
!>      <li>    `STBMV`, `DTBMV`, `CTBMV`, and `ZTBMV` (Matrix-Vector Product for a Triangular Band Matrix, Its Transpose, or Its Conjugate Transpose).<br>
!>  </ol>
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_matrixMulAdd

    use pm_kind, only: SK, IK, LK, RKS, RKD
    use pm_matrixSubset, only: lowDia, lowDia_type
    use pm_matrixSubset, only: uppDia, uppDia_type
    use pm_matrixPack, only: lfpack, lfpack_type
    use pm_matrixClass, only: symmetric, symmetric_type
    use pm_matrixClass, only: hermitian, hermitian_type
    use pm_matrixTrans, only: transSymm, transSymm_type
    use pm_matrixTrans, only: transHerm, transHerm_type
    use pm_blas, only: blasSYMM, blasHEMM, blasGEMM
    use pm_blas, only: blasSPMV, blasHPMV
    use pm_blas, only: blasSYMV, blasHEMV
    use pm_blas, only: blasGEMV

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_matrixMulAdd"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the result of the multiplication of the input matrices/vector `matA` and `matB` in the user-specified form.<br>
    !>
    !>  \details
    !>  Return the result of the multiplication of the input matrices `matA` and `matB` in one of the following forms,
    !>  \f{align*}{
    !>      & \ms{matC} \leftarrow \alpha \ms{matA}~~ \ms{matB} + \beta \ms{matC} && \ms{matC} \leftarrow \alpha \ms{matA}~~ \ms{matB}^T + \beta \ms{matC} && \ms{matC} \leftarrow \alpha \ms{matA}~~ \ms{matB}^H + \beta \ms{matC} \\
    !>      & \ms{matC} \leftarrow \alpha \ms{matA}^T \ms{matB} + \beta \ms{matC} && \ms{matC} \leftarrow \alpha \ms{matA}^T \ms{matB}^T + \beta \ms{matC} && \ms{matC} \leftarrow \alpha \ms{matA}^T \ms{matB}^H + \beta \ms{matC} \\
    !>      & \ms{matC} \leftarrow \alpha \ms{matA}^H \ms{matB} + \beta \ms{matC} && \ms{matC} \leftarrow \alpha \ms{matA}^H \ms{matB}^T + \beta \ms{matC} && \ms{matC} \leftarrow \alpha \ms{matA}^H \ms{matB}^H + \beta \ms{matC}
    !>  \f}
    !>  where \f$\cdot^T\f$ represents a Symmetric transpose, and \f$\cdot^H\f$ represents a Hermitian transpose.<br>
    !>
    !>  The following figure illustrates the form of the **general matrix-matrix-matrix or matrix-vector multiplication** depending on the input values.<br>
    !>  <ul>
    !>      <li>    Default Multiplication (no transposition involved).<br>
    !>              \image html pm_matrixMulAdd@TNA_TNB.png width=1000
    !>              <br>
    !>      <li>    Same as the default but with `operationA` set to [transSymm](@ref pm_matrixTrans::transSymm) or [transHerm](@ref pm_matrixTrans::transHerm).<br>
    !>              \image html pm_matrixMulAdd@TSA_TNB.png width=1000
    !>              <br>
    !>      <li>    Same as the default but with `operationB` set to [transSymm](@ref pm_matrixTrans::transSymm) or [transHerm](@ref pm_matrixTrans::transHerm).<br>
    !>              \image html pm_matrixMulAdd@TNA_TSB.png width=1000
    !>  </ul>
    !>
    !>  \note
    !>  For triangular matrix-matrix or matrix-vector multiplications, see [pm_matrixMulTri](@ref pm_matrixMulTri).<br>
    !>
    !>  \param[in]      matA        :   The input `contiguous` matrix of either shape `(:,:)` or `(:)`
    !>                                  (in [Linear Full Packing (LFP)](@ref pm_matrixPack) known as the **LAPACK packed storage** format) of either,
    !>                                  <ol>
    !>                                      <li>    type `integer` of kind \IKALL,
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing the subset of matrix \f$\ms{matA}\f$ to be used in the matrix-matrix or matrix-vector multiplication.<br>
    !>                                  <ol>
    !>                                      <li>    If the input argument `packA` is present, then `matA` must be a vector of shape `(1 : ndim * (ndim + 1) / 2)`
    !>                                              containing the upper/lower triangle of a Symmetric/Hermitian matrix in [Linear Full Packing (LFP)](@ref pm_matrixPack)
    !>                                              format (also known as the **LAPACK packed storage format**) to be used in the multiplication.<br>
    !>                                      <li>    If `matA` is a matrix (of shape `(:,:)`), then `matA` is assumed to contain a general full matrix or a subset of a Symmetric/Hermitian
    !>                                              matrix in regular ([Rectangular Default Packing (RDP)](@ref pm_matrixPack) format with the optionally-specified `subsetA` format.<br>
    !>                                  </ol>
    !>  \param[in]      classA      :   The input scalar that can be either,
    !>                                  <ol>
    !>                                      <li>    the constant [symmetric](@ref pm_matrixClass::symmetric) for `matA` of any type of any kind,
    !>                                              implying that `matA` is a Symmetric matrix with the specified `subsetA`.
    !>                                      <li>    the constant [hermitian](@ref pm_matrixClass::hermitian) for `matA` of any type of any kind,
    !>                                              implying that `matA` is a Hermitian (conjugate) matrix with the specified `subsetA`.
    !>                                  </ol>
    !>                                  This argument is merely a convenience to resolve the different procedure functionalities within this generic interface.<br>
    !>                                  (**optional**, default = [general matrix](@ref pm_matrixClass::genrecmat_type). It must be present **if and only if** the input argument `subsetA` is also present.)
    !>  \param[in]      subsetA     :   The input scalar that can be either,
    !>                                  <ol>
    !>                                      <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) implying that
    !>                                              the upper-diagonal triangular subset of the input (Symmetric/Hermitian) `matA` should be used in the multiplication without referencing the lower triangle.<br>
    !>                                      <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) implying that
    !>                                              the lower-diagonal triangular subset of the input (Symmetric/Hermitian) `matA` should be used in the multiplication without referencing the upper triangle.<br>
    !>                                  </ol>
    !>                                  This optional argument is merely a convenience to resolve the different procedure functionalities within this generic interface.<br>
    !>                                  If missing, the input `matA` is assumed to a generic matrix whose full specified subset will be used in multiplication.<br>
    !>                                  (**optional**. It can must be present **if and only if** the input argument `classA` is present.)
    !>  \param[in]      packA       :   The input scalar that can be:
    !>                                  <ol>
    !>                                      <li>    The constant [lfpack](@ref pm_matrixPack::lfpack) signifying the [Linear Full Packing (LFP)](@ref pm_matrixPack)
    !>                                              format (also known as the **LAPACK packed storage format** of the input subset of the matrix represented by `matA`.<br>
    !>                                  </ol>
    !>                                  (**optional**, default = [rdpack](@ref pm_matrixPack::rdpack) (Rectangular Default Packing). It must be present **if and only if** the input arguments `classA` and `subsetA` are both present.)
    !>  \param[in]      operationA  :   The input scalar that can be either,
    !>                                  <ol>
    !>                                      <li>    the constant [transSymm](@ref pm_matrixTrans::transSymm) for `matA` of any type of any kind,
    !>                                              implying that a Symmetric transpose of the specified subset of `matA` is to be used in the matrix-matrix or matrix-vector multiplication.
    !>                                      <li>    the constant [transHerm](@ref pm_matrixTrans::transHerm) for `matA` of any type of any kind,
    !>                                              implying that a Hermitian transpose of the specified subset of `matA` is to be used in the matrix-matrix or matrix-vector multiplication.
    !>                                  </ol>
    !>                                  Specifying this argument changes the shape of the subset of `matA` used in the matrix-matrix or matrix-vector multiplication.<br>
    !>                                  See the description of the input argument `ndum`.<br>
    !>                                  This argument is merely a convenience to resolve the different procedure functionalities within this generic interface.<br>
    !>                                  (**optional**. If missing, the specified subset of `matA` will be used as is, without any operations performed on it.)
    !>  \param[in]      matB        :   The input `contiguous` vector of shape `(:)` or matrix of shape `(:,:)` of the same type and kind as `matA`,
    !>                                  containing the subset of vector/matrix \f$\ms{matB}\f$ in the matrix-matrix or matrix-vector multiplication.<br>
    !>  \param[in]      classB      :   The input scalar that can be either,
    !>                                  <ol>
    !>                                      <li>    the constant [symmetric](@ref pm_matrixClass::symmetric) for `matB` of any type of any kind,
    !>                                              implying that `matB` is a Symmetric matrix with the specified `subsetB`.
    !>                                      <li>    the constant [hermitian](@ref pm_matrixClass::hermitian) for `matB` of any type of any kind,
    !>                                              implying that `matB` is a Hermitian (conjugate) matrix with the specified `subsetB`.
    !>                                  </ol>
    !>                                  This argument is merely a convenience to resolve the different procedure functionalities within this generic interface.<br>
    !>                                  (**optional**. It must be present **if and only if** the input argument `subsetB` is also present.)
    !>  \param[in]      subsetB     :   The input scalar that can be either,
    !>                                  <ol>
    !>                                      <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) implying that
    !>                                              the upper-diagonal triangular block of the input Symmetric/Hermitian `matB` should be used in the multiplication without referencing the lower triangle.
    !>                                      <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) implying that
    !>                                              the lower-diagonal triangular block of the input Symmetric/Hermitian `matB` should be used in the multiplication without referencing the upper triangle.
    !>                                  </ol>
    !>                                  This optional argument is merely a convenience to resolve the different procedure functionalities within this generic interface.<br>
    !>                                  If missing, the input `matB` is assumed to a generic matrix whose full specified subset will be used in multiplication.<br>
    !>                                  (**optional**. It can be present **only if** either `classB` is present.)
    !>  \param[in]      operationB  :   The input scalar that can be either,
    !>                                  <ol>
    !>                                      <li>    the constant [transSymm](@ref pm_matrixTrans::transSymm) for `matB` of any type of any kind,
    !>                                              implying that a Symmetric transpose of the specified subset of `matB` is to be used in the matrix-matrix or matrix-vector multiplication.
    !>                                      <li>    the constant [transHerm](@ref pm_matrixTrans::transHerm) for `matB` of any type of any kind,
    !>                                              implying that a Hermitian transpose of the specified subset of `matB` is to be used in the matrix-matrix or matrix-vector multiplication.
    !>                                  </ol>
    !>                                  Specifying this argument changes the shape of the subset of `matB` used in the matrix-matrix or matrix-vector multiplication.<br>
    !>                                  See the description of the input argument `ndum`.<br>
    !>                                  This argument is merely a convenience to resolve the different procedure functionalities within this generic interface.<br>
    !>                                  (**optional**. If missing, the specified subset of `matB` will be used as is, without any operations performed on it.)
    !>  \param[inout]   matC        :   The input/output `contiguous` vector of shape `(:)` or matrix of shape `(:,:)` of the same type and kind as `matA`,
    !>                                  containing the subset of vector/matrix \f$\ms{matC}\f$ in the matrix-matrix or matrix-vector multiplication.<br>
    !>  \param[in]      alpha       :   The input scalar of the same type and kind as `matA`, containing the coefficient \f$\alpha\f$ in the matrix-matrix or matrix-vector multiplication.<br>
    !>                                  (**optional**, default = `1`. It must be present if any of the input arguments `nrow`, `ncol`, `ndum`, `roffA`, `coffA`, `roffB`, `coffB`, `roffC`, `coffC`, `incB`, `incC` are also present.)
    !                                   It cannot be present if `matA` is a scalar value, in which case, `alpha` can be multiplied and merged with `matA` prior to calling this interface.)
    !>  \param[in]      beta        :   The input scalar of the same type and kind as `matA`, containing the coefficient \f$\beta\f$ in the matrix-matrix or matrix-vector multiplication.<br>
    !>                                  (**optional**, default = `1`. It must be present if any of the input arguments `nrow`, `ncol`, `ndum`, `roffA`, `coffA`, `roffB`, `coffB`, `roffC`, `coffC`, `incB`, `incC` are also present.)
    !>  \param[in]      nrow        :   The input non-negative scalar of type `integer` of default kind \IK,
    !>                                  containing the number of rows of `matC` used in the matrix-matrix or matrix-vector multiplication `matC(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol)`.<br>
    !>                                  (**optional**, default = `size(matC, 1)`. It must be present if any of the input arguments `nrow`, `ncol`, `ndum`, `roffA`, `coffA`, `roffB`, `coffB`, `roffC`, `coffC`, `incB`, `incC` are also present.)
    !>  \param[in]      ncol        :   The input non-negative scalar of type `integer` of default kind \IK,
    !>                                  containing the number of columns of `matC` used in the matrix-matrix or matrix-vector multiplication `matC(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol)`.<br>
    !>                                  (**optional**, default = `size(matC, 2)`. It must be present if any of the input arguments `nrow`, `ncol`, `ndum`, `roffA`, `coffA`, `roffB`, `coffB`, `roffC`, `coffC`, `incB`, `incC` are also present.)
    !>  \param[in]      ndum        :   The input non-negative scalar of type `integer` of default kind \IK,
    !>                                  containing the number of dummy rows or columns of `matA` and `matB` used in the matrix-matrix or matrix-vector multiplication.<br>
    !>                                  <ol>
    !>                                      <li>    For matrix `matA`,
    !>                                              <ol>
    !>                                                  <li>    If the input argument `operationA` is missing, then the subset `matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum)` is used in the matrix-matrix or matrix-vector multiplication.
    !>                                                  <li>    If the input argument `operationA` is [transSymm](@ref pm_matrixTrans::transSymm), then the subset `transpose(matA(roffA + 1 : roffA + ndum, coffA + 1 : coffA + nrow)` is used in the matrix-matrix or matrix-vector multiplication.
    !>                                                  <li>    If the input argument `operationA` is [transHerm](@ref pm_matrixTrans::transHerm), then the subset `conjg(transpose(matA(roffA + 1 : roffA + ndum, coffA + 1 : coffA + nrow))` is used in the matrix-matrix or matrix-vector multiplication.
    !>                                              </ol>
    !>                                      <li>    For matrix `matB`,
    !>                                              <ol>
    !>                                                  <li>    If the input argument `operationB` is missing, then the subset `matA(roffB + 1 : roffB + nrow, coffB + 1 : coffB + ndum)` is used in the matrix-matrix or matrix-vector multiplication.
    !>                                                  <li>    If the input argument `operationB` is [transSymm](@ref pm_matrixTrans::transSymm), then the subset `transpose(matA(roffB + 1 : roffB + ndum, coffB + 1 : coffB + ncol)` is used in the matrix-matrix or matrix-vector multiplication.
    !>                                                  <li>    If the input argument `operationB` is [transHerm](@ref pm_matrixTrans::transHerm), then the subset `conjg(transpose(matA(roffB + 1 : roffB + ncol, coffB + 1 : coffB + ndum))` is used in the matrix-matrix or matrix-vector multiplication.
    !>                                              </ol>
    !>                                  </ol>
    !>                                  (**optional**. It must be missing if the input arguments `subsetA` or `subsetB`  are present, that is, when either `matA` or `matB` is a square upper/lower-triangular Symmetric matrix.)
    !>  \param[in]      ndim        :   The input non-negative scalar of type `integer` of default kind \IK,
    !>                                  containing the number of rows and columns of the square Symmetric/Hermitian subset of `matA` used in the matrix-matrix or matrix-vector multiplication `matA(roffA + 1 : roffA + ndim, coffA + 1 : coffA + ndim)`.<br>
    !>                                  (**optional**, default = `size(matA, 1)`. It can be present **only if** the arguments `matB` and `matC` are of rank `1` and the input arguments `incB` and `incC` are present.)
    !>  \param[in]      roffA       :   The input non-negative scalar of type `integer` of default kind \IK containing the <b>off</b>set from the first <b>r</b>ow of the input matrix `matA`,
    !>                                  such that `matA(1 + roffA, 1 + coffA)` marks the top-left corner of the block of `matA` used in the matrix-matrix or matrix-vector multiplication.<br>
    !>                                  (**optional**, default = `0`)
    !>  \param[in]      coffA       :   The input non-negative scalar of type `integer` of default kind \IK containing the <b>off</b>set from the first <b>c</b>olumn of the input matrix `matA`,
    !>                                  such that `matA(1 + roffA, 1 + coffA)` marks the top-left corner of the block of `matA` used in the matrix-matrix or matrix-vector multiplication.<br>
    !>                                  (**optional**, default = `0`)
    !>  \param[in]      roffB       :   The input non-negative scalar of type `integer` of default kind \IK containing the <b>off</b>set from the first <b>r</b>ow of the input matrix `matB`,
    !>                                  such that `matB(1 + roffB, 1 + coffB)` marks the top-left corner of the block of `matB` used in the matrix-matrix or matrix-vector multiplication.<br>
    !>                                  (**optional**, default = `0`)
    !>  \param[in]      coffB       :   The input non-negative scalar of type `integer` of default kind \IK containing the <b>off</b>set from the first <b>c</b>olumn of the input matrix `matB`,
    !>                                  such that `matB(1 + roffB, 1 + coffB)` marks the top-left corner of the block of `matB` used in the matrix-matrix or matrix-vector multiplication.<br>
    !>                                  (**optional**, default = `0`)
    !>  \param[in]      roffC       :   The input non-negative scalar of type `integer` of default kind \IK containing the <b>off</b>set from the first <b>r</b>ow of the input matrix `matC`,
    !>                                  such that `matC(1 + roffC, 1 + coffC)` marks the top-left corner of the block of `matC` used in the matrix-matrix or matrix-vector multiplication.<br>
    !>                                  (**optional**, default = `0`)
    !>  \param[in]      coffC       :   The input non-negative scalar of type `integer` of default kind \IK containing the <b>off</b>set from the first <b>c</b>olumn of the input matrix `matC`,
    !>                                  such that `matC(1 + roffC, 1 + coffC)` marks the top-left corner of the block of `matC` used in the matrix-matrix or matrix-vector multiplication.<br>
    !>                                  (**optional**, default = `0`)
    !>  \param[in]      incB        :   The input non-zero scalar `integer` of default kind \IK, containing the stride of the input argument `matB(:)`.<br>
    !>                                  <ol>
    !>                                      <li>    A positive value implies the multiplication to be performed on the subset `matB(1 : 1 + (ndim - 1) * incB : incB)`.
    !>                                      <li>    A negative value implies the multiplication to be performed on the subset `matB(1 + (1 - ndim) * incB : 1 : incB)`.
    !>                                  </ol>
    !>                                  (**optional**, default = `1`. It can be present only if `matB` and `matC` are of rank `1`.)
    !>  \param[in]      incC        :   The input non-zero scalar `integer` of default kind \IK, containing the stride of the input/output argument `matC(:)`.<br>
    !>                                  <ol>
    !>                                      <li>    A positive value implies the multiplication to be performed on the subset `matC(1 : 1 + (ndim - 1) * incC : incC)`.<br>
    !>                                      <li>    A negative value implies the multiplication to be performed on the subset `matC(1 + (1 - ndim) * incC : 1 : incC)`.<br>
    !>                                  </ol>
    !>                                  (**optional**, default = `1`. It can be present only if `matB` and `matC` are of rank `1`.)
    !>
    !>  \interface{setMatMulAdd}
    !>  \code{.F90}
    !>
    !>      use pm_matrixMulAdd, only: setMatMulAdd
    !>      use pm_matrixMulAdd, only: uppDiaA, lowDiaA
    !>      use pm_matrixMulAdd, only: uppDiaB, lowDiaB
    !>      use pm_matrixMulAdd, only: transSymm, transHerm
    !>      use pm_matrixMulAdd, only: symmetric, hermitian
    !>
    !>      ! BLAS - LEVEL 2: ?GEMV - SIMPLIFIED INTERFACE.
    !>
    !>      call setMatMulAdd(matA(:,:), matB(:), matC(:), alpha = alpha, beta = beta)
    !>      call setMatMulAdd(matA(:,:), operationA, matB(:), matC(:), alpha = alpha, beta = beta) ! operationA = transSymm/transHerm
    !>
    !>      ! BLAS - LEVEL 2: ?GEMV - contiguous interface.
    !>
    !>      call setMatMulAdd(matA(:,:), matB(:), matC(:), alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
    !>      call setMatMulAdd(matA(:,:), operationA, matB(:), matC(:), alpha, beta, nrow, ncol, roffA, coffA, incB, incC) ! operationA = transSymm/transHerm
    !>
    !>      ! BLAS - LEVEL 3: ?GEMM - SIMPLIFIED INTERFACE.
    !>
    !>      call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha = alpha, beta = beta)
    !>      call setMatMulAdd(matA(:,:), operationA, matB(:,:), matC(:,:), alpha = alpha, beta = beta) ! operationA = transSymm/transHerm
    !>      call setMatMulAdd(matA(:,:), matB(:,:), operationB, matC(:,:), alpha = alpha, beta = beta) ! operationB = transSymm/transHerm
    !>      call setMatMulAdd(matA(:,:), operationA, matB(:,:), operationB, matC(:,:), alpha = alpha, beta = beta) ! operationA/operationB = transSymm/transHerm
    !>
    !>      ! BLAS - LEVEL 3: ?GEMM - contiguous interface.
    !>
    !>      call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
    !>      call setMatMulAdd(matA(:,:), operationA, matB(:,:), matC(:,:), alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC) ! operationA = transSymm/transHerm
    !>      call setMatMulAdd(matA(:,:), matB(:,:), operationB, matC(:,:), alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC) ! operationB = transSymm/transHerm
    !>      call setMatMulAdd(matA(:,:), operationA, matB(:,:), operationB, matC(:,:), alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC) ! operationA/operationB = transSymm/transHerm
    !>
    !>      ! BLAS - LEVEL 2: ?SYMV / ?HEMV - SIMPLIFIED INTERFACE.
    !>
    !>      call setMatMulAdd(matA(:,:), classA, subsetA, matB(:), matC(:), alpha = alpha, beta = beta) ! classA = symmetric/hermitian, subsetA = uppDia/lowDia
    !>
    !>      ! BLAS - LEVEL 2: ?SYMV / ?HEMV - contiguous interface.
    !>
    !>      call setMatMulAdd(matA(:,:), classA, subsetA, matB(:), matC(:), alpha, beta, ndim, roffA, coffA, incB, incC) ! classA = symmetric/hermitian, subsetA = uppDia/lowDia
    !>
    !>      ! BLAS - LEVEL 3: ?SYMM / ?HEMM - SIMPLIFIED INTERFACE.
    !>
    !>      call setMatMulAdd(matA(:,:), classA, subsetA, matB(:,:), matC(:,:), alpha = alpha, beta = beta) ! classA = symmetric/hermitian, subsetA = uppDia/lowDia
    !>      call setMatMulAdd(matA(:,:), matB(:,:), classB, subsetB, matC(:,:), alpha = alpha, beta = beta) ! classB = symmetric/hermitian, subsetB = uppDia/lowDia
    !>
    !>      ! BLAS - LEVEL 3: ?SYMM / ?HEMM - contiguous interface.
    !>
    !>      call setMatMulAdd(matA(:,:), classA, subsetA, matB(:,:), matC(:,:), alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC) ! classA = symmetric/hermitian, subsetA = uppDia/lowDia
    !>      call setMatMulAdd(matA(:,:), matB(:,:), classB, subsetB, matC(:,:), alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC) ! classB = symmetric/hermitian, subsetB = uppDia/lowDia
    !>
    !>      ! BLAS - LEVEL 2: ?SPMV / ?HPMV - SIMPLIFIED INTERFACE.
    !>
    !>      call setMatMulAdd(matA(:), classA, subsetA, packA, matB(:), matC(:), alpha = alpha, beta = beta) ! classA = symmetric/hermitian, subsetA = uppDia/lowDia, packA = lfpack
    !>
    !>      ! BLAS - LEVEL 2: ?SPMV / ?HPMV - contiguous interface.
    !>
    !>      call setMatMulAdd(matA(:), classA, subsetA, packA, matB(:), matC(:), alpha, beta, ndim, incB, incC) ! classA = symmetric/hermitian, subsetA = uppDia/lowDia, packA = lfpack
    !>      !
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 <= ndum` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 /= [incB, incC])` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= [nrow, ncol])` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= [roffA, coffA])` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= [roffB, coffB])` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= [roffC, coffC])` must hold for the corresponding input arguments.<br>
    !>  The condition `all(0 <= [ndim, roffA, coffA])` must hold for the corresponding input arguments.<br>
    !>  The condition `abs(incB) * (ndim - 1) < size(matB(:))` must hold for the corresponding input arguments.<br>
    !>  The condition `abs(incC) * (ndim - 1) < size(matC(:))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(ndim + [roffA, coffA] <= shape(matA))` must hold for the corresponding input arguments.<br>
    !>  The condition `all(ndim == [shape(matA(:,:)), shape(matB(:)), shape(matC(:))])` must hold for the corresponding input arguments.<br>
    !>  The condition `all([roffA + nrow, coffA + ndum] <= shape(matA))` must hold for the corresponding input arguments when the input argument `operationA` is missing.<br>
    !>  The condition `all([roffA + ndum, coffA + nrow] <= shape(matA))` must hold for the corresponding input arguments when the input argument `operationA` is present.<br>
    !>  The condition `all([roffB + ndum, coffB + ncol] <= shape(matB))` must hold for the corresponding input arguments when the input argument `operationB` is missing.<br>
    !>  The condition `all([roffB + ncol, coffB + ndum] <= shape(matB))` must hold for the corresponding input arguments when the input argument `operationB` is present.<br>
    !>  All inequalities in the above conditions become equalities when the optional shape arguments `ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC` are missing.<br>
    !>  There are many more runtime bound and sanity checks that are performed within the procedures, but not listed above.<br>
    !>  \vericons
    !>
    !>  \warnpure
    !>
    !>  \lapack{3.11}
    !>  `SGEMV`, `DGEMV`, `CGEMV`, `ZGEMV`,
    !>  `SSPMV`, `DSPMV`, `CSPMV`, `ZSPMV`,
    !>  `SHPMV`, `DHPMV`, `CHPMV`, `ZHPMV`,
    !>  `SSYMV`, `DSYMV`, `CSYMV`, `ZSYMV`,
    !>  `SHEMV`, `DHEMV`, `CHEMV`, `ZHEMV`,
    !>  `SSYMM`, `DSYMM`, `CSYMM`, `ZSYMM`,
    !>  `SHEMM`, `DHEMM`, `CHEMM`, `ZHEMM`,
    !>  `SGEMM`, `DGEMM`, `CGEMM`, `ZGEMM`.<br>
    !>  In particular multiplications of matrices of type `integer` of arbitrary kinds are also possible.<br>
    !>
    !   \lapackint{setMatMulAdd}
    ! 
    !   \warning This table must be updated with the correct modern interfaces.
    !   BLAS/LAPACK interface                                                                                                                                                       |   [setMatMulAdd](@ref pm_matrixMulAdd::setMatMulAdd) interface
    !   ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------
    !   `call SSPMV / DSPMV / CHPMV / ZHPMV (uplo = "U", n = ndim, alpha, ap = matA(i:), x = matB(k), incx = incB, beta, y = matC(p), incy = incC)`                                 |   `call setMatMulAdd(matA(i :), classA = symmetric, subsetA = uppDiaA, matB(k :), matC(p :), alpha, beta, ndim, incB, incC)`
    !   `call SSPMV / DSPMV / CHPMV / ZHPMV (uplo = "L", n = ndim, alpha, ap = matA(i:), x = matB(k), incx = incB, beta, y = matC(p), incy = incC)`                                 |   `call setMatMulAdd(matA(i :), classA = symmetric, subsetA = lowDiaA, matB(k :), matC(p :), alpha, beta, ndim, incB, incC)`
    !   `call                 CHPMV / ZHPMV (uplo = "U", n = ndim, alpha, ap = matA(i:), x = matB(k), incx = incB, beta, y = matC(p), incy = incC)`                                 |   `call setMatMulAdd(matA(i :), classA = hermitian, subsetA = uppDiaA, matB(k :), matC(p :), alpha, beta, ndim, incB, incC)`
    !   `call                 CHPMV / ZHPMV (uplo = "L", n = ndim, alpha, ap = matA(i:), x = matB(k), incx = incB, beta, y = matC(p), incy = incC)`                                 |   `call setMatMulAdd(matA(i :), classA = hermitian, subsetA = lowDiaA, matB(k :), matC(p :), alpha, beta, ndim, incB, incC)`
    !   `call SSYMV / DSYMV / CHEMV / ZHEMV (uplo = "U", n = ndim, alpha, a = matA(i,j), lda, x = matB(k:), incx = incB, beta, y = matC(p:), incy = incC)`                          |   `call setMatMulAdd(matA(:,:), classA = symmetric, subsetA = uppDiaA, matB(k :), matC(p :), alpha, beta, ndim, roffA = i - 1, coffA = j - 1, incB, incC)`
    !   `call SSYMV / DSYMV / CHEMV / ZHEMV (uplo = "L", n = ndim, alpha, a = matA(i,j), lda, x = matB(k:), incx = incB, beta, y = matC(p:), incy = incC)`                          |   `call setMatMulAdd(matA(:,:), classA = symmetric, subsetA = lowDiaA, matB(k :), matC(p :), alpha, beta, ndim, roffA = i - 1, coffA = j - 1, incB, incC)`
    !   `call                 CHEMV / ZHEMV (uplo = "U", n = ndim, alpha, a = matA(i,j), lda, x = matB(k:), incx = incB, beta, y = matC(p:), incy = incC)`                          |   `call setMatMulAdd(matA(:,:), classA = hermitian, subsetA = uppDiaA, matB(k :), matC(p :), alpha, beta, ndim, roffA = i - 1, coffA = j - 1, incB, incC)`
    !   `call                 CHEMV / ZHEMV (uplo = "L", n = ndim, alpha, a = matA(i,j), lda, x = matB(k:), incx = incB, beta, y = matC(p:), incy = incC)`                          |   `call setMatMulAdd(matA(:,:), classA = hermitian, subsetA = lowDiaA, matB(k :), matC(p :), alpha, beta, ndim, roffA = i - 1, coffA = j - 1, incB, incC)`
    !   `call SGEMM / DGEMM / CGEMM / ZGEMM (transa = "N", transb = "N", m = nrow, n = ncol, k = ndum, alpha, a = matA(i,j), lda, b = matB(k,l), ldb, beta, c = matC(p,q), ldc)`    |   `call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha, beta, nrow, ncol, ndum, roffA = i - 1, coffA = j - 1, roffB = k - 1, coffB = l - 1, roffC = p - 1, coffC = q - 1)`
    !   `call SGEMM / DGEMM / CGEMM / ZGEMM (transa = "N", transb = "T", m = nrow, n = ncol, k = ndum, alpha, a = matA(i,j), lda, b = matB(k,l), ldb, beta, c = matC(p,q), ldc)`    |   `call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha, beta, nrow, ncol, ndum, roffA = i - 1, coffA = j - 1, roffB = k - 1, coffB = l - 1, roffC = p - 1, coffC = q - 1, operationB = transSymm)`
    !   `call SGEMM / DGEMM / CGEMM / ZGEMM (transa = "T", transb = "N", m = nrow, n = ncol, k = ndum, alpha, a = matA(i,j), lda, b = matB(k,l), ldb, beta, c = matC(p,q), ldc)`    |   `call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha, beta, nrow, ncol, ndum, roffA = i - 1, coffA = j - 1, roffB = k - 1, coffB = l - 1, roffC = p - 1, coffC = q - 1, operationA = transSymm)`
    !   `call SGEMM / DGEMM / CGEMM / ZGEMM (transa = "T", transb = "T", m = nrow, n = ncol, k = ndum, alpha, a = matA(i,j), lda, b = matB(k,l), ldb, beta, c = matC(p,q), ldc)`    |   `call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha, beta, nrow, ncol, ndum, roffA = i - 1, coffA = j - 1, roffB = k - 1, coffB = l - 1, roffC = p - 1, coffC = q - 1, operationA = transSymm, operationB = transSymm)`
    !   `call                 CGEMM / ZGEMM (transa = "N", transb = "C", m = nrow, n = ncol, k = ndum, alpha, a = matA(i,j), lda, b = matB(k,l), ldb, beta, c = matC(p,q), ldc)`    |   `call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha, beta, nrow, ncol, ndum, roffA = i - 1, coffA = j - 1, roffB = k - 1, coffB = l - 1, roffC = p - 1, coffC = q - 1, operationB = transHerm)`
    !   `call                 CGEMM / ZGEMM (transa = "C", transb = "N", m = nrow, n = ncol, k = ndum, alpha, a = matA(i,j), lda, b = matB(k,l), ldb, beta, c = matC(p,q), ldc)`    |   `call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha, beta, nrow, ncol, ndum, roffA = i - 1, coffA = j - 1, roffB = k - 1, coffB = l - 1, roffC = p - 1, coffC = q - 1, operationA = transHerm)`
    !   `call                 CGEMM / ZGEMM (transa = "T", transb = "C", m = nrow, n = ncol, k = ndum, alpha, a = matA(i,j), lda, b = matB(k,l), ldb, beta, c = matC(p,q), ldc)`    |   `call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha, beta, nrow, ncol, ndum, roffA = i - 1, coffA = j - 1, roffB = k - 1, coffB = l - 1, roffC = p - 1, coffC = q - 1, operationA = transSymm, operationB = transHerm)`
    !   `call                 CGEMM / ZGEMM (transa = "C", transb = "T", m = nrow, n = ncol, k = ndum, alpha, a = matA(i,j), lda, b = matB(k,l), ldb, beta, c = matC(p,q), ldc)`    |   `call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha, beta, nrow, ncol, ndum, roffA = i - 1, coffA = j - 1, roffB = k - 1, coffB = l - 1, roffC = p - 1, coffC = q - 1, operationA = transHerm, operationB = transSymm)`
    !   `call                 CGEMM / ZGEMM (transa = "C", transb = "T", m = nrow, n = ncol, k = ndum, alpha, a = matA(i,j), lda, b = matB(k,l), ldb, beta, c = matC(p,q), ldc)`    |   `call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha, beta, nrow, ncol, ndum, roffA = i - 1, coffA = j - 1, roffB = k - 1, coffB = l - 1, roffC = p - 1, coffC = q - 1, operationA = transHerm, operationB = transHerm)`
    !   `call SSYMM / DSYMM / CSYMM / ZSYMM (side = "L", uplo = "U", m = nrow, n = ncol, alpha, a = matA(i,j), lda, b = matB(k,l), ldb, beta, c = matC(p,q), ldc)`                  |   `call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha, beta, nrow, ncol, roffA = i - 1, coffA = j - 1, roffB = k - 1, coffB = l - 1, roffC = p - 1, coffC = q - 1, subsetA = uppDiaA, classA = symmetric)`
    !   `call SSYMM / DSYMM / CSYMM / ZSYMM (side = "R", uplo = "U", m = nrow, n = ncol, alpha, a = matB(k,l), lda, b = matA(i,j), ldb, beta, c = matC(p,q), ldc)`                  |   `call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha, beta, nrow, ncol, roffA = i - 1, coffA = j - 1, roffB = k - 1, coffB = l - 1, roffC = p - 1, coffC = q - 1, subsetB = uppDiaB, classB = symmetric)`
    !   `call SSYMM / DSYMM / CSYMM / ZSYMM (side = "L", uplo = "L", m = nrow, n = ncol, alpha, a = matA(i,j), lda, b = matB(k,l), ldb, beta, c = matC(p,q), ldc)`                  |   `call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha, beta, nrow, ncol, roffA = i - 1, coffA = j - 1, roffB = k - 1, coffB = l - 1, roffC = p - 1, coffC = q - 1, subsetA = lowDiaA, classA = symmetric)`
    !   `call SSYMM / DSYMM / CSYMM / ZSYMM (side = "R", uplo = "L", m = nrow, n = ncol, alpha, a = matB(k,l), lda, b = matA(i,j), ldb, beta, c = matC(p,q), ldc)`                  |   `call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha, beta, nrow, ncol, roffA = i - 1, coffA = j - 1, roffB = k - 1, coffB = l - 1, roffC = p - 1, coffC = q - 1, subsetB = lowDiaB, classB = symmetric)`
    !   `call                 CHEMM / ZHEMM (side = "L", uplo = "U", m = nrow, n = ncol, alpha, a = matA(i,j), lda, b = matB(k,l), ldb, beta, c = matC(p,q), ldc)`                  |   `call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha, beta, nrow, ncol, roffA = i - 1, coffA = j - 1, roffB = k - 1, coffB = l - 1, roffC = p - 1, coffC = q - 1, subsetA = uppDiaA, classA = hermitian)`
    !   `call                 CHEMM / ZHEMM (side = "R", uplo = "U", m = nrow, n = ncol, alpha, a = matB(k,l), lda, b = matA(i,j), ldb, beta, c = matC(p,q), ldc)`                  |   `call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha, beta, nrow, ncol, roffA = i - 1, coffA = j - 1, roffB = k - 1, coffB = l - 1, roffC = p - 1, coffC = q - 1, subsetB = uppDiaB, classB = hermitian)`
    !   `call                 CHEMM / ZHEMM (side = "L", uplo = "L", m = nrow, n = ncol, alpha, a = matA(i,j), lda, b = matB(k,l), ldb, beta, c = matC(p,q), ldc)`                  |   `call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha, beta, nrow, ncol, roffA = i - 1, coffA = j - 1, roffB = k - 1, coffB = l - 1, roffC = p - 1, coffC = q - 1, subsetA = lowDiaA, classA = hermitian)`
    !   `call                 CHEMM / ZHEMM (side = "R", uplo = "L", m = nrow, n = ncol, alpha, a = matB(k,l), lda, b = matA(i,j), ldb, beta, c = matC(p,q), ldc)`                  |   `call setMatMulAdd(matA(:,:), matB(:,:), matC(:,:), alpha, beta, nrow, ncol, roffA = i - 1, coffA = j - 1, roffB = k - 1, coffB = l - 1, roffC = p - 1, coffC = q - 1, subsetB = lowDiaB, classB = hermitian)`
    !>
    !>  BLAS/LAPACK interface arguments |   [setMatMulAdd](@ref pm_matrixMulAdd::setMatMulAdd) interface arguments
    !>  --------------------------------|-------------------------------------------------------------------------------------------
    !>  `side = "L"`                    |   `classA` argument is present and `classB` argument is missing.
    !>  `side = "B"`                    |   `classA` argument is missing and `classB` argument is present.
    !>  `uplo = "U"`                    |   `subsetA = uppDiaA` or `subsetB = uppDiaB`.
    !>  `uplo = "L"`                    |   `subsetA = lowDiaA` or `subsetB = lowDiaB`.
    !>  `transa = "N"`                  |   `operationA` argument is missing.
    !>  `transa = "T"`                  |   `operationA =` [transSymm](@ref pm_matrixTrans::transSymm)
    !>  `transb = "C"`                  |   `operationA =` [transHerm](@ref pm_matrixTrans::transHerm)
    !>  `transb = "N"`                  |   `operationB` argument is missing.
    !>  `transb = "T"`                  |   `operationB =` [transSymm](@ref pm_matrixTrans::transSymm)
    !>  `transb = "C"`                  |   `operationB =` [transHerm](@ref pm_matrixTrans::transHerm)
    !>  `m`                             |   `nrow`
    !>  `n`                             |   `ncol`
    !>  `k`                             |   `ndum`
    !>  `alpha`                         |   `alpha`
    !>  `a = A(i,j)`                    |   `matA = A, roffA = i - 1, coffA = j - 1`
    !>  `lda`                           |   NONE (passed implicitly).
    !>  `b = B(i,j)`                    |   `matB = B, roffB = i - 1, coffB = j - 1`
    !>  `ldb`                           |   NONE (passed implicitly).
    !>  `beta`                          |   `beta`
    !>  `c = C(i,j)`                    |   `matC = C, roffC = i - 1, coffC = j - 1`
    !>  `ldc`                           |   NONE (passed implicitly).
    !>
    !>  \see
    !>  [lowDia](@ref pm_matrixSubset::lowDia)<br>
    !>  [uppDia](@ref pm_matrixSubset::uppDia)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [transSymm](@ref pm_matrixTrans::transSymm)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [setMatMulTri](@ref pm_matrixMulTri::setMatMulTri)<br>
    !>
    !>  \example{setMatMulAdd}
    !>  \include{lineno} example/pm_matrixMulAdd/setMatMulAdd/main.F90
    !>  \compilef{setMatMulAdd}
    !>  \output{setMatMulAdd}
    !>  \include{lineno} example/pm_matrixMulAdd/setMatMulAdd/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixMulAdd](@ref test_pm_matrixMulAdd)
    !>
    !>  \naming
    !>  \code{.F90}
    !>      XXXX_ASS_CNA_CNB_SFA_SFB_TNA_TNB_IK5()
    !>      |||| ||| ||| ||| ||| ||| ||| ||| |||
    !>      |||| ||| ||| ||| ||| ||| ||| ||| The type and kind of the input matrix arguments.
    !>      |||| ||| ||| ||| ||| ||| ||| The transposition of `matB`: TNB/TSB/THB => Transposition None / Transposition Symmetric / Transposition Hermitian
    !>      |||| ||| ||| ||| ||| ||| The transposition of `matA`: TNA/TSA/THA => Transposition None / Transposition Symmetric / Transposition Hermitian
    !>      |||| ||| ||| ||| ||| The subset of `matB` used: SFB/SUB/SLB => Subset Full / Subset Upper / Subset Lower
    !>      |||| ||| ||| ||| The subset of `matA` used: SFA/SUA/SLA => Subset Full / Subset Upper / Subset Lower
    !>      |||| ||| ||| The class of `matB`: CNB/CSB/CHB => General (None) / Symmetric / Hermitian
    !>      |||| ||| The class of `matA`: CNA/CSA/CHA => General (None) / Symmetric / Hermitian
    !>      |||| The explicitness of the array bounds: ASS/EXP => ASSumed bounds passed / EXPlicit bounds passed
    !>      The BLAS routine implemented: gemv, symv, hemv, gemm, symm, hemm, spmv, hpmv
    !>  \endcode
    !>
    !>  \finmain{setMatMulAdd}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX

    ! BLAS 2: `SGEMV`, `DGEMV`, `CGEMV`, `ZGEMV`

    interface setMatMulAdd

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TNA_TNB_IK5(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TNA_TNB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TNA_TNB_IK4(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TNA_TNB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TNA_TNB_IK3(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TNA_TNB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TNA_TNB_IK2(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TNA_TNB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TNA_TNB_IK1(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TNA_TNB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TNA_TNB_CK5(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TNA_TNB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TNA_TNB_CK4(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TNA_TNB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TNA_TNB_CK3(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TNA_TNB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TNA_TNB_CK2(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TNA_TNB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TNA_TNB_CK1(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TNA_TNB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TNA_TNB_RK5(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TNA_TNB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TNA_TNB_RK4(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TNA_TNB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TNA_TNB_RK3(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TNA_TNB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TNA_TNB_RK2(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TNA_TNB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TNA_TNB_RK1(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TNA_TNB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TSA_TNB_IK5(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TSA_TNB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TSA_TNB_IK4(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TSA_TNB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TSA_TNB_IK3(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TSA_TNB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TSA_TNB_IK2(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TSA_TNB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TSA_TNB_IK1(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TSA_TNB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TSA_TNB_CK5(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TSA_TNB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TSA_TNB_CK4(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TSA_TNB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TSA_TNB_CK3(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TSA_TNB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TSA_TNB_CK2(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TSA_TNB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TSA_TNB_CK1(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TSA_TNB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TSA_TNB_RK5(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TSA_TNB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TSA_TNB_RK4(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TSA_TNB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TSA_TNB_RK3(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TSA_TNB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TSA_TNB_RK2(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TSA_TNB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_TSA_TNB_RK1(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_TSA_TNB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_THA_TNB_IK5(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_THA_TNB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_THA_TNB_IK4(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_THA_TNB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_THA_TNB_IK3(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_THA_TNB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_THA_TNB_IK2(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_THA_TNB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_THA_TNB_IK1(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_THA_TNB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_THA_TNB_CK5(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_THA_TNB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_THA_TNB_CK4(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_THA_TNB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_THA_TNB_CK3(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_THA_TNB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_THA_TNB_CK2(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_THA_TNB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_THA_TNB_CK1(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_THA_TNB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_THA_TNB_RK5(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_THA_TNB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_THA_TNB_RK4(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_THA_TNB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_THA_TNB_RK3(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_THA_TNB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_THA_TNB_RK2(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_THA_TNB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemv_ASS_SFA_SFB_THA_TNB_RK1(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_ASS_SFA_SFB_THA_TNB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
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
    PURE module subroutine gemv_EXP_SFA_SFB_TNA_TNB_IK5(matA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TNA_TNB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TNA_TNB_IK4(matA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TNA_TNB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TNA_TNB_IK3(matA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TNA_TNB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TNA_TNB_IK2(matA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TNA_TNB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TNA_TNB_IK1(matA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TNA_TNB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TNA_TNB_CK5(matA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TNA_TNB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TNA_TNB_CK4(matA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TNA_TNB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TNA_TNB_CK3(matA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TNA_TNB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TNA_TNB_CK2(matA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TNA_TNB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TNA_TNB_CK1(matA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TNA_TNB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TNA_TNB_RK5(matA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TNA_TNB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TNA_TNB_RK4(matA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TNA_TNB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TNA_TNB_RK3(matA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TNA_TNB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TNA_TNB_RK2(matA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TNA_TNB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TNA_TNB_RK1(matA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TNA_TNB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TSA_TNB_IK5(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TSA_TNB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TSA_TNB_IK4(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TSA_TNB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TSA_TNB_IK3(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TSA_TNB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TSA_TNB_IK2(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TSA_TNB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TSA_TNB_IK1(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TSA_TNB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TSA_TNB_CK5(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TSA_TNB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TSA_TNB_CK4(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TSA_TNB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TSA_TNB_CK3(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TSA_TNB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TSA_TNB_CK2(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TSA_TNB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TSA_TNB_CK1(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TSA_TNB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TSA_TNB_RK5(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TSA_TNB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TSA_TNB_RK4(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TSA_TNB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TSA_TNB_RK3(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TSA_TNB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TSA_TNB_RK2(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TSA_TNB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_TSA_TNB_RK1(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_TSA_TNB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_THA_TNB_IK5(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_THA_TNB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_THA_TNB_IK4(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_THA_TNB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_THA_TNB_IK3(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_THA_TNB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_THA_TNB_IK2(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_THA_TNB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_THA_TNB_IK1(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_THA_TNB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_THA_TNB_CK5(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_THA_TNB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_THA_TNB_CK4(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_THA_TNB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_THA_TNB_CK3(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_THA_TNB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_THA_TNB_CK2(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_THA_TNB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_THA_TNB_CK1(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_THA_TNB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_THA_TNB_RK5(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_THA_TNB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_THA_TNB_RK4(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_THA_TNB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_THA_TNB_RK3(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_THA_TNB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_THA_TNB_RK2(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_THA_TNB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemv_EXP_SFA_SFB_THA_TNB_RK1(matA, operationA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemv_EXP_SFA_SFB_THA_TNB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! BLAS 2: `SSPMV`, `DSPMV`, `CSPMV`, `ZSPMV`

    interface setMatMulAdd

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine spmv_ASS_CSA_SUA_CNB_SFB_IK5(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SUA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine spmv_ASS_CSA_SUA_CNB_SFB_IK4(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SUA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine spmv_ASS_CSA_SUA_CNB_SFB_IK3(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SUA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine spmv_ASS_CSA_SUA_CNB_SFB_IK2(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SUA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine spmv_ASS_CSA_SUA_CNB_SFB_IK1(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SUA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine spmv_ASS_CSA_SUA_CNB_SFB_CK5(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SUA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine spmv_ASS_CSA_SUA_CNB_SFB_CK4(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SUA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine spmv_ASS_CSA_SUA_CNB_SFB_CK3(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SUA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine spmv_ASS_CSA_SUA_CNB_SFB_CK2(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SUA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine spmv_ASS_CSA_SUA_CNB_SFB_CK1(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SUA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine spmv_ASS_CSA_SUA_CNB_SFB_RK5(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SUA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine spmv_ASS_CSA_SUA_CNB_SFB_RK4(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SUA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine spmv_ASS_CSA_SUA_CNB_SFB_RK3(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SUA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine spmv_ASS_CSA_SUA_CNB_SFB_RK2(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SUA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine spmv_ASS_CSA_SUA_CNB_SFB_RK1(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SUA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine spmv_ASS_CSA_SLA_CNB_SFB_IK5(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SLA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine spmv_ASS_CSA_SLA_CNB_SFB_IK4(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SLA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine spmv_ASS_CSA_SLA_CNB_SFB_IK3(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SLA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine spmv_ASS_CSA_SLA_CNB_SFB_IK2(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SLA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine spmv_ASS_CSA_SLA_CNB_SFB_IK1(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SLA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine spmv_ASS_CSA_SLA_CNB_SFB_CK5(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SLA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine spmv_ASS_CSA_SLA_CNB_SFB_CK4(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SLA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine spmv_ASS_CSA_SLA_CNB_SFB_CK3(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SLA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine spmv_ASS_CSA_SLA_CNB_SFB_CK2(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SLA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine spmv_ASS_CSA_SLA_CNB_SFB_CK1(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SLA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine spmv_ASS_CSA_SLA_CNB_SFB_RK5(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SLA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine spmv_ASS_CSA_SLA_CNB_SFB_RK4(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SLA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine spmv_ASS_CSA_SLA_CNB_SFB_RK3(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SLA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine spmv_ASS_CSA_SLA_CNB_SFB_RK2(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SLA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine spmv_ASS_CSA_SLA_CNB_SFB_RK1(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_ASS_CSA_SLA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
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
    PURE module subroutine spmv_EXP_CSA_SUA_CNB_SFB_IK5(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SUA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine spmv_EXP_CSA_SUA_CNB_SFB_IK4(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SUA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine spmv_EXP_CSA_SUA_CNB_SFB_IK3(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SUA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine spmv_EXP_CSA_SUA_CNB_SFB_IK2(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SUA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine spmv_EXP_CSA_SUA_CNB_SFB_IK1(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SUA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine spmv_EXP_CSA_SUA_CNB_SFB_CK5(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SUA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine spmv_EXP_CSA_SUA_CNB_SFB_CK4(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SUA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine spmv_EXP_CSA_SUA_CNB_SFB_CK3(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SUA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine spmv_EXP_CSA_SUA_CNB_SFB_CK2(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SUA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine spmv_EXP_CSA_SUA_CNB_SFB_CK1(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SUA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine spmv_EXP_CSA_SUA_CNB_SFB_RK5(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SUA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine spmv_EXP_CSA_SUA_CNB_SFB_RK4(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SUA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine spmv_EXP_CSA_SUA_CNB_SFB_RK3(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SUA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine spmv_EXP_CSA_SUA_CNB_SFB_RK2(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SUA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine spmv_EXP_CSA_SUA_CNB_SFB_RK1(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SUA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine spmv_EXP_CSA_SLA_CNB_SFB_IK5(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SLA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine spmv_EXP_CSA_SLA_CNB_SFB_IK4(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SLA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine spmv_EXP_CSA_SLA_CNB_SFB_IK3(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SLA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine spmv_EXP_CSA_SLA_CNB_SFB_IK2(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SLA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine spmv_EXP_CSA_SLA_CNB_SFB_IK1(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SLA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine spmv_EXP_CSA_SLA_CNB_SFB_CK5(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SLA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine spmv_EXP_CSA_SLA_CNB_SFB_CK4(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SLA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine spmv_EXP_CSA_SLA_CNB_SFB_CK3(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SLA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine spmv_EXP_CSA_SLA_CNB_SFB_CK2(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SLA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine spmv_EXP_CSA_SLA_CNB_SFB_CK1(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SLA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine spmv_EXP_CSA_SLA_CNB_SFB_RK5(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SLA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine spmv_EXP_CSA_SLA_CNB_SFB_RK4(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SLA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine spmv_EXP_CSA_SLA_CNB_SFB_RK3(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SLA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine spmv_EXP_CSA_SLA_CNB_SFB_RK2(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SLA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine spmv_EXP_CSA_SLA_CNB_SFB_RK1(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: spmv_EXP_CSA_SLA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! BLAS 2: `SHPMV`, `DHPMV`, `CHPMV`, `ZHPMV`

    interface setMatMulAdd

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SUA_CNB_SFB_IK5(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SUA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SUA_CNB_SFB_IK4(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SUA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SUA_CNB_SFB_IK3(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SUA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SUA_CNB_SFB_IK2(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SUA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SUA_CNB_SFB_IK1(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SUA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SUA_CNB_SFB_CK5(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SUA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SUA_CNB_SFB_CK4(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SUA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SUA_CNB_SFB_CK3(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SUA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SUA_CNB_SFB_CK2(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SUA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SUA_CNB_SFB_CK1(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SUA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SUA_CNB_SFB_RK5(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SUA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SUA_CNB_SFB_RK4(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SUA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SUA_CNB_SFB_RK3(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SUA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SUA_CNB_SFB_RK2(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SUA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SUA_CNB_SFB_RK1(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SUA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SLA_CNB_SFB_IK5(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SLA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SLA_CNB_SFB_IK4(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SLA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SLA_CNB_SFB_IK3(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SLA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SLA_CNB_SFB_IK2(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SLA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SLA_CNB_SFB_IK1(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SLA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SLA_CNB_SFB_CK5(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SLA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SLA_CNB_SFB_CK4(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SLA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SLA_CNB_SFB_CK3(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SLA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SLA_CNB_SFB_CK2(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SLA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SLA_CNB_SFB_CK1(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SLA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SLA_CNB_SFB_RK5(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SLA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SLA_CNB_SFB_RK4(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SLA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SLA_CNB_SFB_RK3(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SLA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SLA_CNB_SFB_RK2(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SLA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine hpmv_ASS_CHA_SLA_CNB_SFB_RK1(matA, classA, subsetA, packA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_ASS_CHA_SLA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
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
    PURE module subroutine hpmv_EXP_CHA_SUA_CNB_SFB_IK5(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SUA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SUA_CNB_SFB_IK4(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SUA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SUA_CNB_SFB_IK3(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SUA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SUA_CNB_SFB_IK2(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SUA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SUA_CNB_SFB_IK1(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SUA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SUA_CNB_SFB_CK5(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SUA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SUA_CNB_SFB_CK4(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SUA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SUA_CNB_SFB_CK3(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SUA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SUA_CNB_SFB_CK2(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SUA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SUA_CNB_SFB_CK1(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SUA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SUA_CNB_SFB_RK5(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SUA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SUA_CNB_SFB_RK4(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SUA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SUA_CNB_SFB_RK3(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SUA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SUA_CNB_SFB_RK2(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SUA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SUA_CNB_SFB_RK1(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SUA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SLA_CNB_SFB_IK5(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SLA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SLA_CNB_SFB_IK4(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SLA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SLA_CNB_SFB_IK3(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SLA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SLA_CNB_SFB_IK2(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SLA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SLA_CNB_SFB_IK1(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SLA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SLA_CNB_SFB_CK5(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SLA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SLA_CNB_SFB_CK4(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SLA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SLA_CNB_SFB_CK3(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SLA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SLA_CNB_SFB_CK2(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SLA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SLA_CNB_SFB_CK1(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SLA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SLA_CNB_SFB_RK5(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SLA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SLA_CNB_SFB_RK4(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SLA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SLA_CNB_SFB_RK3(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SLA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SLA_CNB_SFB_RK2(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SLA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine hpmv_EXP_CHA_SLA_CNB_SFB_RK1(matA, classA, subsetA, packA, matB, matC, alpha, beta, ndim, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hpmv_EXP_CHA_SLA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
        type(lfpack_type)       , intent(in)                    :: packA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! BLAS 2: `SSYMV`, `DSYMV`, `CSYMV`, `ZSYMV`

    interface setMatMulAdd

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine symv_ASS_CSA_SUA_CNB_SFB_IK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SUA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine symv_ASS_CSA_SUA_CNB_SFB_IK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SUA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine symv_ASS_CSA_SUA_CNB_SFB_IK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SUA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine symv_ASS_CSA_SUA_CNB_SFB_IK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SUA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine symv_ASS_CSA_SUA_CNB_SFB_IK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SUA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine symv_ASS_CSA_SUA_CNB_SFB_CK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SUA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine symv_ASS_CSA_SUA_CNB_SFB_CK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SUA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine symv_ASS_CSA_SUA_CNB_SFB_CK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SUA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine symv_ASS_CSA_SUA_CNB_SFB_CK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SUA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine symv_ASS_CSA_SUA_CNB_SFB_CK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SUA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine symv_ASS_CSA_SUA_CNB_SFB_RK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SUA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine symv_ASS_CSA_SUA_CNB_SFB_RK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SUA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine symv_ASS_CSA_SUA_CNB_SFB_RK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SUA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine symv_ASS_CSA_SUA_CNB_SFB_RK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SUA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine symv_ASS_CSA_SUA_CNB_SFB_RK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SUA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine symv_ASS_CSA_SLA_CNB_SFB_IK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SLA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine symv_ASS_CSA_SLA_CNB_SFB_IK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SLA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine symv_ASS_CSA_SLA_CNB_SFB_IK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SLA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine symv_ASS_CSA_SLA_CNB_SFB_IK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SLA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine symv_ASS_CSA_SLA_CNB_SFB_IK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SLA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine symv_ASS_CSA_SLA_CNB_SFB_CK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SLA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine symv_ASS_CSA_SLA_CNB_SFB_CK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SLA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine symv_ASS_CSA_SLA_CNB_SFB_CK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SLA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine symv_ASS_CSA_SLA_CNB_SFB_CK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SLA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine symv_ASS_CSA_SLA_CNB_SFB_CK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SLA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine symv_ASS_CSA_SLA_CNB_SFB_RK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SLA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine symv_ASS_CSA_SLA_CNB_SFB_RK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SLA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine symv_ASS_CSA_SLA_CNB_SFB_RK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SLA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine symv_ASS_CSA_SLA_CNB_SFB_RK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SLA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine symv_ASS_CSA_SLA_CNB_SFB_RK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_ASS_CSA_SLA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
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
    PURE module subroutine symv_EXP_CSA_SUA_CNB_SFB_IK5(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SUA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine symv_EXP_CSA_SUA_CNB_SFB_IK4(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SUA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine symv_EXP_CSA_SUA_CNB_SFB_IK3(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SUA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine symv_EXP_CSA_SUA_CNB_SFB_IK2(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SUA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine symv_EXP_CSA_SUA_CNB_SFB_IK1(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SUA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine symv_EXP_CSA_SUA_CNB_SFB_CK5(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SUA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine symv_EXP_CSA_SUA_CNB_SFB_CK4(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SUA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine symv_EXP_CSA_SUA_CNB_SFB_CK3(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SUA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine symv_EXP_CSA_SUA_CNB_SFB_CK2(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SUA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine symv_EXP_CSA_SUA_CNB_SFB_CK1(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SUA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine symv_EXP_CSA_SUA_CNB_SFB_RK5(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SUA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine symv_EXP_CSA_SUA_CNB_SFB_RK4(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SUA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine symv_EXP_CSA_SUA_CNB_SFB_RK3(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SUA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine symv_EXP_CSA_SUA_CNB_SFB_RK2(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SUA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine symv_EXP_CSA_SUA_CNB_SFB_RK1(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SUA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine symv_EXP_CSA_SLA_CNB_SFB_IK5(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SLA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine symv_EXP_CSA_SLA_CNB_SFB_IK4(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SLA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine symv_EXP_CSA_SLA_CNB_SFB_IK3(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SLA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine symv_EXP_CSA_SLA_CNB_SFB_IK2(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SLA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine symv_EXP_CSA_SLA_CNB_SFB_IK1(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SLA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine symv_EXP_CSA_SLA_CNB_SFB_CK5(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SLA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine symv_EXP_CSA_SLA_CNB_SFB_CK4(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SLA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine symv_EXP_CSA_SLA_CNB_SFB_CK3(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SLA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine symv_EXP_CSA_SLA_CNB_SFB_CK2(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SLA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine symv_EXP_CSA_SLA_CNB_SFB_CK1(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SLA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine symv_EXP_CSA_SLA_CNB_SFB_RK5(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SLA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine symv_EXP_CSA_SLA_CNB_SFB_RK4(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SLA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine symv_EXP_CSA_SLA_CNB_SFB_RK3(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SLA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine symv_EXP_CSA_SLA_CNB_SFB_RK2(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SLA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine symv_EXP_CSA_SLA_CNB_SFB_RK1(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symv_EXP_CSA_SLA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! BLAS 2: `SHEMV`, `DHEMV`, `CHEMV`, `ZHEMV`

    interface setMatMulAdd

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine hemv_ASS_CHA_SUA_CNB_SFB_IK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SUA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine hemv_ASS_CHA_SUA_CNB_SFB_IK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SUA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine hemv_ASS_CHA_SUA_CNB_SFB_IK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SUA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine hemv_ASS_CHA_SUA_CNB_SFB_IK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SUA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine hemv_ASS_CHA_SUA_CNB_SFB_IK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SUA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine hemv_ASS_CHA_SUA_CNB_SFB_CK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SUA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine hemv_ASS_CHA_SUA_CNB_SFB_CK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SUA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine hemv_ASS_CHA_SUA_CNB_SFB_CK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SUA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine hemv_ASS_CHA_SUA_CNB_SFB_CK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SUA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine hemv_ASS_CHA_SUA_CNB_SFB_CK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SUA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine hemv_ASS_CHA_SUA_CNB_SFB_RK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SUA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine hemv_ASS_CHA_SUA_CNB_SFB_RK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SUA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine hemv_ASS_CHA_SUA_CNB_SFB_RK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SUA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine hemv_ASS_CHA_SUA_CNB_SFB_RK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SUA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine hemv_ASS_CHA_SUA_CNB_SFB_RK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SUA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine hemv_ASS_CHA_SLA_CNB_SFB_IK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SLA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine hemv_ASS_CHA_SLA_CNB_SFB_IK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SLA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine hemv_ASS_CHA_SLA_CNB_SFB_IK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SLA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine hemv_ASS_CHA_SLA_CNB_SFB_IK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SLA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine hemv_ASS_CHA_SLA_CNB_SFB_IK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SLA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine hemv_ASS_CHA_SLA_CNB_SFB_CK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SLA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine hemv_ASS_CHA_SLA_CNB_SFB_CK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SLA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine hemv_ASS_CHA_SLA_CNB_SFB_CK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SLA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine hemv_ASS_CHA_SLA_CNB_SFB_CK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SLA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine hemv_ASS_CHA_SLA_CNB_SFB_CK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SLA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine hemv_ASS_CHA_SLA_CNB_SFB_RK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SLA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine hemv_ASS_CHA_SLA_CNB_SFB_RK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SLA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine hemv_ASS_CHA_SLA_CNB_SFB_RK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SLA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine hemv_ASS_CHA_SLA_CNB_SFB_RK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SLA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine hemv_ASS_CHA_SLA_CNB_SFB_RK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_ASS_CHA_SLA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
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
    PURE module subroutine hemv_EXP_CHA_SUA_CNB_SFB_IK5(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SUA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine hemv_EXP_CHA_SUA_CNB_SFB_IK4(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SUA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine hemv_EXP_CHA_SUA_CNB_SFB_IK3(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SUA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine hemv_EXP_CHA_SUA_CNB_SFB_IK2(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SUA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine hemv_EXP_CHA_SUA_CNB_SFB_IK1(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SUA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine hemv_EXP_CHA_SUA_CNB_SFB_CK5(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SUA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine hemv_EXP_CHA_SUA_CNB_SFB_CK4(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SUA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine hemv_EXP_CHA_SUA_CNB_SFB_CK3(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SUA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine hemv_EXP_CHA_SUA_CNB_SFB_CK2(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SUA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine hemv_EXP_CHA_SUA_CNB_SFB_CK1(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SUA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine hemv_EXP_CHA_SUA_CNB_SFB_RK5(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SUA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine hemv_EXP_CHA_SUA_CNB_SFB_RK4(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SUA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine hemv_EXP_CHA_SUA_CNB_SFB_RK3(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SUA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine hemv_EXP_CHA_SUA_CNB_SFB_RK2(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SUA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine hemv_EXP_CHA_SUA_CNB_SFB_RK1(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SUA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine hemv_EXP_CHA_SLA_CNB_SFB_IK5(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SLA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine hemv_EXP_CHA_SLA_CNB_SFB_IK4(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SLA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine hemv_EXP_CHA_SLA_CNB_SFB_IK3(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SLA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine hemv_EXP_CHA_SLA_CNB_SFB_IK2(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SLA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine hemv_EXP_CHA_SLA_CNB_SFB_IK1(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SLA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine hemv_EXP_CHA_SLA_CNB_SFB_CK5(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SLA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine hemv_EXP_CHA_SLA_CNB_SFB_CK4(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SLA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine hemv_EXP_CHA_SLA_CNB_SFB_CK3(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SLA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine hemv_EXP_CHA_SLA_CNB_SFB_CK2(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SLA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine hemv_EXP_CHA_SLA_CNB_SFB_CK1(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SLA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine hemv_EXP_CHA_SLA_CNB_SFB_RK5(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SLA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine hemv_EXP_CHA_SLA_CNB_SFB_RK4(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SLA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine hemv_EXP_CHA_SLA_CNB_SFB_RK3(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SLA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine hemv_EXP_CHA_SLA_CNB_SFB_RK2(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SLA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine hemv_EXP_CHA_SLA_CNB_SFB_RK1(matA, classA, subsetA, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemv_EXP_CHA_SLA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: ndim, roffA, coffA, incB, incC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! BLAS 3: `SSYMM`, `DSYMM`, `CSYMM`, `ZSYMM`

    interface setMatMulAdd

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine symm_ASS_CSA_SUA_CNB_SFB_IK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SUA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine symm_ASS_CSA_SUA_CNB_SFB_IK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SUA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine symm_ASS_CSA_SUA_CNB_SFB_IK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SUA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine symm_ASS_CSA_SUA_CNB_SFB_IK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SUA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine symm_ASS_CSA_SUA_CNB_SFB_IK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SUA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine symm_ASS_CSA_SUA_CNB_SFB_CK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SUA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine symm_ASS_CSA_SUA_CNB_SFB_CK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SUA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine symm_ASS_CSA_SUA_CNB_SFB_CK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SUA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine symm_ASS_CSA_SUA_CNB_SFB_CK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SUA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine symm_ASS_CSA_SUA_CNB_SFB_CK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SUA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine symm_ASS_CSA_SUA_CNB_SFB_RK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SUA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine symm_ASS_CSA_SUA_CNB_SFB_RK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SUA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine symm_ASS_CSA_SUA_CNB_SFB_RK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SUA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine symm_ASS_CSA_SUA_CNB_SFB_RK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SUA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine symm_ASS_CSA_SUA_CNB_SFB_RK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SUA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine symm_ASS_CSA_SLA_CNB_SFB_IK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SLA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine symm_ASS_CSA_SLA_CNB_SFB_IK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SLA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine symm_ASS_CSA_SLA_CNB_SFB_IK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SLA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine symm_ASS_CSA_SLA_CNB_SFB_IK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SLA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine symm_ASS_CSA_SLA_CNB_SFB_IK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SLA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine symm_ASS_CSA_SLA_CNB_SFB_CK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SLA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine symm_ASS_CSA_SLA_CNB_SFB_CK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SLA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine symm_ASS_CSA_SLA_CNB_SFB_CK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SLA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine symm_ASS_CSA_SLA_CNB_SFB_CK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SLA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine symm_ASS_CSA_SLA_CNB_SFB_CK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SLA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine symm_ASS_CSA_SLA_CNB_SFB_RK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SLA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine symm_ASS_CSA_SLA_CNB_SFB_RK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SLA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine symm_ASS_CSA_SLA_CNB_SFB_RK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SLA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine symm_ASS_CSA_SLA_CNB_SFB_RK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SLA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine symm_ASS_CSA_SLA_CNB_SFB_RK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CSA_SLA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SUB_IK5(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SUB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SUB_IK4(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SUB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SUB_IK3(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SUB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SUB_IK2(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SUB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SUB_IK1(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SUB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SUB_CK5(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SUB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SUB_CK4(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SUB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SUB_CK3(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SUB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SUB_CK2(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SUB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SUB_CK1(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SUB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SUB_RK5(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SUB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SUB_RK4(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SUB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SUB_RK3(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SUB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SUB_RK2(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SUB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SUB_RK1(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SUB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SLB_IK5(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SLB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SLB_IK4(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SLB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SLB_IK3(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SLB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SLB_IK2(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SLB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SLB_IK1(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SLB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SLB_CK5(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SLB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SLB_CK4(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SLB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SLB_CK3(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SLB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SLB_CK2(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SLB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SLB_CK1(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SLB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SLB_RK5(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SLB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SLB_RK4(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SLB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SLB_RK3(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SLB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SLB_RK2(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SLB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine symm_ASS_CNA_SFA_CSB_SLB_RK1(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_ASS_CNA_SFA_CSB_SLB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
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
    PURE module subroutine symm_EXP_CSA_SUA_CNB_SFB_IK5(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SUA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine symm_EXP_CSA_SUA_CNB_SFB_IK4(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SUA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine symm_EXP_CSA_SUA_CNB_SFB_IK3(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SUA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine symm_EXP_CSA_SUA_CNB_SFB_IK2(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SUA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine symm_EXP_CSA_SUA_CNB_SFB_IK1(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SUA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine symm_EXP_CSA_SUA_CNB_SFB_CK5(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SUA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine symm_EXP_CSA_SUA_CNB_SFB_CK4(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SUA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine symm_EXP_CSA_SUA_CNB_SFB_CK3(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SUA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine symm_EXP_CSA_SUA_CNB_SFB_CK2(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SUA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine symm_EXP_CSA_SUA_CNB_SFB_CK1(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SUA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine symm_EXP_CSA_SUA_CNB_SFB_RK5(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SUA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine symm_EXP_CSA_SUA_CNB_SFB_RK4(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SUA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine symm_EXP_CSA_SUA_CNB_SFB_RK3(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SUA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine symm_EXP_CSA_SUA_CNB_SFB_RK2(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SUA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine symm_EXP_CSA_SUA_CNB_SFB_RK1(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SUA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine symm_EXP_CSA_SLA_CNB_SFB_IK5(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SLA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine symm_EXP_CSA_SLA_CNB_SFB_IK4(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SLA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine symm_EXP_CSA_SLA_CNB_SFB_IK3(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SLA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine symm_EXP_CSA_SLA_CNB_SFB_IK2(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SLA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine symm_EXP_CSA_SLA_CNB_SFB_IK1(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SLA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine symm_EXP_CSA_SLA_CNB_SFB_CK5(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SLA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine symm_EXP_CSA_SLA_CNB_SFB_CK4(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SLA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine symm_EXP_CSA_SLA_CNB_SFB_CK3(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SLA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine symm_EXP_CSA_SLA_CNB_SFB_CK2(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SLA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine symm_EXP_CSA_SLA_CNB_SFB_CK1(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SLA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine symm_EXP_CSA_SLA_CNB_SFB_RK5(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SLA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine symm_EXP_CSA_SLA_CNB_SFB_RK4(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SLA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine symm_EXP_CSA_SLA_CNB_SFB_RK3(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SLA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine symm_EXP_CSA_SLA_CNB_SFB_RK2(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SLA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine symm_EXP_CSA_SLA_CNB_SFB_RK1(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CSA_SLA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(symmetric_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SUB_IK5(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SUB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SUB_IK4(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SUB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SUB_IK3(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SUB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SUB_IK2(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SUB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SUB_IK1(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SUB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SUB_CK5(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SUB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SUB_CK4(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SUB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SUB_CK3(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SUB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SUB_CK2(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SUB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SUB_CK1(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SUB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SUB_RK5(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SUB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SUB_RK4(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SUB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SUB_RK3(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SUB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SUB_RK2(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SUB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SUB_RK1(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SUB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SLB_IK5(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SLB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SLB_IK4(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SLB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SLB_IK3(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SLB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SLB_IK2(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SLB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SLB_IK1(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SLB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SLB_CK5(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SLB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SLB_CK4(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SLB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SLB_CK3(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SLB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SLB_CK2(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SLB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SLB_CK1(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SLB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SLB_RK5(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SLB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SLB_RK4(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SLB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SLB_RK3(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SLB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SLB_RK2(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SLB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine symm_EXP_CNA_SFA_CSB_SLB_RK1(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: symm_EXP_CNA_SFA_CSB_SLB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(symmetric_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! BLAS 3: `SHEMM`, `DHEMM`, `CHEMM`, `ZHEMM`

    interface setMatMulAdd

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine hemm_ASS_CHA_SUA_CNB_SFB_IK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SUA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine hemm_ASS_CHA_SUA_CNB_SFB_IK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SUA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine hemm_ASS_CHA_SUA_CNB_SFB_IK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SUA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine hemm_ASS_CHA_SUA_CNB_SFB_IK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SUA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine hemm_ASS_CHA_SUA_CNB_SFB_IK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SUA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine hemm_ASS_CHA_SUA_CNB_SFB_CK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SUA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine hemm_ASS_CHA_SUA_CNB_SFB_CK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SUA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine hemm_ASS_CHA_SUA_CNB_SFB_CK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SUA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine hemm_ASS_CHA_SUA_CNB_SFB_CK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SUA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine hemm_ASS_CHA_SUA_CNB_SFB_CK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SUA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine hemm_ASS_CHA_SUA_CNB_SFB_RK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SUA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine hemm_ASS_CHA_SUA_CNB_SFB_RK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SUA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine hemm_ASS_CHA_SUA_CNB_SFB_RK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SUA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine hemm_ASS_CHA_SUA_CNB_SFB_RK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SUA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine hemm_ASS_CHA_SUA_CNB_SFB_RK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SUA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine hemm_ASS_CHA_SLA_CNB_SFB_IK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SLA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine hemm_ASS_CHA_SLA_CNB_SFB_IK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SLA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine hemm_ASS_CHA_SLA_CNB_SFB_IK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SLA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine hemm_ASS_CHA_SLA_CNB_SFB_IK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SLA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine hemm_ASS_CHA_SLA_CNB_SFB_IK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SLA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine hemm_ASS_CHA_SLA_CNB_SFB_CK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SLA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine hemm_ASS_CHA_SLA_CNB_SFB_CK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SLA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine hemm_ASS_CHA_SLA_CNB_SFB_CK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SLA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine hemm_ASS_CHA_SLA_CNB_SFB_CK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SLA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine hemm_ASS_CHA_SLA_CNB_SFB_CK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SLA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine hemm_ASS_CHA_SLA_CNB_SFB_RK5(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SLA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine hemm_ASS_CHA_SLA_CNB_SFB_RK4(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SLA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine hemm_ASS_CHA_SLA_CNB_SFB_RK3(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SLA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine hemm_ASS_CHA_SLA_CNB_SFB_RK2(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SLA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine hemm_ASS_CHA_SLA_CNB_SFB_RK1(matA, classA, subsetA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CHA_SLA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SUB_IK5(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SUB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SUB_IK4(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SUB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SUB_IK3(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SUB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SUB_IK2(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SUB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SUB_IK1(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SUB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SUB_CK5(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SUB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SUB_CK4(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SUB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SUB_CK3(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SUB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SUB_CK2(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SUB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SUB_CK1(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SUB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SUB_RK5(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SUB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SUB_RK4(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SUB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SUB_RK3(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SUB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SUB_RK2(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SUB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SUB_RK1(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SUB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SLB_IK5(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SLB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SLB_IK4(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SLB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SLB_IK3(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SLB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SLB_IK2(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SLB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SLB_IK1(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SLB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SLB_CK5(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SLB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SLB_CK4(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SLB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SLB_CK3(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SLB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SLB_CK2(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SLB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SLB_CK1(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SLB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SLB_RK5(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SLB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SLB_RK4(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SLB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SLB_RK3(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SLB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SLB_RK2(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SLB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine hemm_ASS_CNA_SFA_CHB_SLB_RK1(matA, matB, classB, subsetB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_ASS_CNA_SFA_CHB_SLB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
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
    PURE module subroutine hemm_EXP_CHA_SUA_CNB_SFB_IK5(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SUA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine hemm_EXP_CHA_SUA_CNB_SFB_IK4(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SUA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine hemm_EXP_CHA_SUA_CNB_SFB_IK3(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SUA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine hemm_EXP_CHA_SUA_CNB_SFB_IK2(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SUA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine hemm_EXP_CHA_SUA_CNB_SFB_IK1(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SUA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine hemm_EXP_CHA_SUA_CNB_SFB_CK5(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SUA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine hemm_EXP_CHA_SUA_CNB_SFB_CK4(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SUA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine hemm_EXP_CHA_SUA_CNB_SFB_CK3(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SUA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine hemm_EXP_CHA_SUA_CNB_SFB_CK2(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SUA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine hemm_EXP_CHA_SUA_CNB_SFB_CK1(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SUA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine hemm_EXP_CHA_SUA_CNB_SFB_RK5(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SUA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine hemm_EXP_CHA_SUA_CNB_SFB_RK4(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SUA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine hemm_EXP_CHA_SUA_CNB_SFB_RK3(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SUA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine hemm_EXP_CHA_SUA_CNB_SFB_RK2(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SUA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine hemm_EXP_CHA_SUA_CNB_SFB_RK1(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SUA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine hemm_EXP_CHA_SLA_CNB_SFB_IK5(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SLA_CNB_SFB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine hemm_EXP_CHA_SLA_CNB_SFB_IK4(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SLA_CNB_SFB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine hemm_EXP_CHA_SLA_CNB_SFB_IK3(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SLA_CNB_SFB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine hemm_EXP_CHA_SLA_CNB_SFB_IK2(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SLA_CNB_SFB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine hemm_EXP_CHA_SLA_CNB_SFB_IK1(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SLA_CNB_SFB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine hemm_EXP_CHA_SLA_CNB_SFB_CK5(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SLA_CNB_SFB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine hemm_EXP_CHA_SLA_CNB_SFB_CK4(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SLA_CNB_SFB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine hemm_EXP_CHA_SLA_CNB_SFB_CK3(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SLA_CNB_SFB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine hemm_EXP_CHA_SLA_CNB_SFB_CK2(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SLA_CNB_SFB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine hemm_EXP_CHA_SLA_CNB_SFB_CK1(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SLA_CNB_SFB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine hemm_EXP_CHA_SLA_CNB_SFB_RK5(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SLA_CNB_SFB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine hemm_EXP_CHA_SLA_CNB_SFB_RK4(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SLA_CNB_SFB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine hemm_EXP_CHA_SLA_CNB_SFB_RK3(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SLA_CNB_SFB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine hemm_EXP_CHA_SLA_CNB_SFB_RK2(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SLA_CNB_SFB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine hemm_EXP_CHA_SLA_CNB_SFB_RK1(matA, classA, subsetA, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CHA_SLA_CNB_SFB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetA
        type(hermitian_type)    , intent(in)                    :: classA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SUB_IK5(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SUB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SUB_IK4(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SUB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SUB_IK3(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SUB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SUB_IK2(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SUB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SUB_IK1(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SUB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SUB_CK5(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SUB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SUB_CK4(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SUB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SUB_CK3(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SUB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SUB_CK2(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SUB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SUB_CK1(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SUB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SUB_RK5(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SUB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SUB_RK4(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SUB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SUB_RK3(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SUB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SUB_RK2(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SUB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SUB_RK1(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SUB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(uppDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SLB_IK5(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SLB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SLB_IK4(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SLB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SLB_IK3(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SLB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SLB_IK2(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SLB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SLB_IK1(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SLB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SLB_CK5(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SLB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SLB_CK4(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SLB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SLB_CK3(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SLB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SLB_CK2(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SLB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SLB_CK1(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SLB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SLB_RK5(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SLB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SLB_RK4(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SLB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SLB_RK3(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SLB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SLB_RK2(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SLB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine hemm_EXP_CNA_SFA_CHB_SLB_RK1(matA, matB, classB, subsetB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: hemm_EXP_CNA_SFA_CHB_SLB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(lowDia_type)       , intent(in)                    :: subsetB
        type(hermitian_type)    , intent(in)                    :: classB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! BLAS 3: `SGEMM`, `DGEMM`, `CGEMM`, `ZGEMM`

    interface setMatMulAdd

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TNB_IK5(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TNB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TNB_IK4(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TNB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TNB_IK3(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TNB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TNB_IK2(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TNB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TNB_IK1(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TNB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TNB_CK5(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TNB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TNB_CK4(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TNB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TNB_CK3(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TNB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TNB_CK2(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TNB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TNB_CK1(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TNB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TNB_RK5(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TNB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TNB_RK4(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TNB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TNB_RK3(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TNB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TNB_RK2(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TNB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TNB_RK1(matA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TNB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TNB_IK5(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TNB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TNB_IK4(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TNB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TNB_IK3(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TNB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TNB_IK2(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TNB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TNB_IK1(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TNB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TNB_CK5(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TNB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TNB_CK4(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TNB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TNB_CK3(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TNB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TNB_CK2(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TNB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TNB_CK1(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TNB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TNB_RK5(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TNB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TNB_RK4(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TNB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TNB_RK3(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TNB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TNB_RK2(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TNB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TNB_RK1(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TNB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TNB_IK5(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TNB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TNB_IK4(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TNB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TNB_IK3(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TNB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TNB_IK2(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TNB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TNB_IK1(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TNB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TNB_CK5(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TNB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TNB_CK4(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TNB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TNB_CK3(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TNB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TNB_CK2(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TNB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TNB_CK1(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TNB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TNB_RK5(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TNB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TNB_RK4(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TNB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TNB_RK3(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TNB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TNB_RK2(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TNB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TNB_RK1(matA, operationA, matB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TNB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TSB_IK5(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TSB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TSB_IK4(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TSB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TSB_IK3(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TSB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TSB_IK2(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TSB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TSB_IK1(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TSB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TSB_CK5(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TSB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TSB_CK4(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TSB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TSB_CK3(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TSB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TSB_CK2(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TSB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TSB_CK1(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TSB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TSB_RK5(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TSB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TSB_RK4(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TSB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TSB_RK3(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TSB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TSB_RK2(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TSB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_TSB_RK1(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_TSB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TSB_IK5(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TSB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TSB_IK4(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TSB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TSB_IK3(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TSB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TSB_IK2(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TSB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TSB_IK1(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TSB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TSB_CK5(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TSB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TSB_CK4(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TSB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TSB_CK3(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TSB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TSB_CK2(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TSB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TSB_CK1(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TSB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TSB_RK5(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TSB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TSB_RK4(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TSB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TSB_RK3(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TSB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TSB_RK2(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TSB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_TSB_RK1(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_TSB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TSB_IK5(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TSB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TSB_IK4(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TSB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TSB_IK3(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TSB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TSB_IK2(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TSB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TSB_IK1(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TSB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TSB_CK5(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TSB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TSB_CK4(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TSB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TSB_CK3(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TSB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TSB_CK2(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TSB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TSB_CK1(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TSB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TSB_RK5(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TSB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TSB_RK4(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TSB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TSB_RK3(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TSB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TSB_RK2(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TSB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_TSB_RK1(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_TSB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_THB_IK5(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_THB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_THB_IK4(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_THB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_THB_IK3(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_THB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_THB_IK2(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_THB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_THB_IK1(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_THB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_THB_CK5(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_THB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_THB_CK4(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_THB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_THB_CK3(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_THB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_THB_CK2(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_THB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_THB_CK1(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_THB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_THB_RK5(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_THB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_THB_RK4(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_THB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_THB_RK3(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_THB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_THB_RK2(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_THB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TNA_THB_RK1(matA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TNA_THB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_THB_IK5(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_THB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_THB_IK4(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_THB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_THB_IK3(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_THB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_THB_IK2(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_THB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_THB_IK1(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_THB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_THB_CK5(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_THB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_THB_CK4(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_THB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_THB_CK3(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_THB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_THB_CK2(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_THB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_THB_CK1(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_THB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_THB_RK5(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_THB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_THB_RK4(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_THB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_THB_RK3(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_THB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_THB_RK2(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_THB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_TSA_THB_RK1(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_TSA_THB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_THB_IK5(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_THB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_THB_IK4(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_THB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_THB_IK3(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_THB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_THB_IK2(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_THB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_THB_IK1(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_THB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IKC)            , intent(in)    , optional      :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(:,:)
        integer(IKC)            , intent(in)    , contiguous    :: matB(:,:)
        integer(IKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_THB_CK5(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_THB_CK5
#endif
        use pm_kind, only: CKC => CK5
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_THB_CK4(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_THB_CK4
#endif
        use pm_kind, only: CKC => CK4
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_THB_CK3(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_THB_CK3
#endif
        use pm_kind, only: CKC => CK3
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_THB_CK2(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_THB_CK2
#endif
        use pm_kind, only: CKC => CK2
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_THB_CK1(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_THB_CK1
#endif
        use pm_kind, only: CKC => CK1
        complex(CKC)            , intent(in)    , optional      :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(:,:)
        complex(CKC)            , intent(in)    , contiguous    :: matB(:,:)
        complex(CKC)            , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_THB_RK5(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_THB_RK5
#endif
        use pm_kind, only: RKC => RK5
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_THB_RK4(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_THB_RK4
#endif
        use pm_kind, only: RKC => RK4
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_THB_RK3(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_THB_RK3
#endif
        use pm_kind, only: RKC => RK3
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_THB_RK2(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_THB_RK2
#endif
        use pm_kind, only: RKC => RK2
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemm_ASS_SFA_SFB_THA_THB_RK1(matA, operationA, matB, operationB, matC, alpha, beta)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_ASS_SFA_SFB_THA_THB_RK1
#endif
        use pm_kind, only: RKC => RK1
        real(RKC)               , intent(in)    , optional      :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(:,:)
        real(RKC)               , intent(in)    , contiguous    :: matB(:,:)
        real(RKC)               , intent(inout) , contiguous    :: matC(:,:)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
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
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TNB_IK5(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TNB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TNB_IK4(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TNB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TNB_IK3(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TNB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TNB_IK2(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TNB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TNB_IK1(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TNB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TNB_CK5(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TNB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TNB_CK4(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TNB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TNB_CK3(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TNB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TNB_CK2(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TNB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TNB_CK1(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TNB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TNB_RK5(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TNB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TNB_RK4(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TNB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TNB_RK3(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TNB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TNB_RK2(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TNB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TNB_RK1(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TNB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TNB_IK5(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TNB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TNB_IK4(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TNB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TNB_IK3(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TNB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TNB_IK2(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TNB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TNB_IK1(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TNB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TNB_CK5(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TNB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TNB_CK4(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TNB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TNB_CK3(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TNB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TNB_CK2(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TNB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TNB_CK1(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TNB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TNB_RK5(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TNB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TNB_RK4(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TNB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TNB_RK3(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TNB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TNB_RK2(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TNB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TNB_RK1(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TNB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TNB_IK5(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TNB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TNB_IK4(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TNB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TNB_IK3(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TNB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TNB_IK2(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TNB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TNB_IK1(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TNB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TNB_CK5(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TNB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TNB_CK4(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TNB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TNB_CK3(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TNB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TNB_CK2(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TNB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TNB_CK1(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TNB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TNB_RK5(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TNB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TNB_RK4(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TNB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TNB_RK3(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TNB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TNB_RK2(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TNB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TNB_RK1(matA, operationA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TNB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TSB_IK5(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TSB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TSB_IK4(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TSB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TSB_IK3(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TSB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TSB_IK2(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TSB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TSB_IK1(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TSB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TSB_CK5(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TSB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TSB_CK4(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TSB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TSB_CK3(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TSB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TSB_CK2(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TSB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TSB_CK1(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TSB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TSB_RK5(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TSB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TSB_RK4(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TSB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TSB_RK3(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TSB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TSB_RK2(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TSB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_TSB_RK1(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_TSB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TSB_IK5(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TSB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TSB_IK4(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TSB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TSB_IK3(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TSB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TSB_IK2(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TSB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TSB_IK1(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TSB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TSB_CK5(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TSB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TSB_CK4(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TSB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TSB_CK3(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TSB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TSB_CK2(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TSB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TSB_CK1(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TSB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TSB_RK5(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TSB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TSB_RK4(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TSB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TSB_RK3(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TSB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TSB_RK2(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TSB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_TSB_RK1(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_TSB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TSB_IK5(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TSB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TSB_IK4(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TSB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TSB_IK3(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TSB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TSB_IK2(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TSB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TSB_IK1(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TSB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TSB_CK5(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TSB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TSB_CK4(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TSB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TSB_CK3(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TSB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TSB_CK2(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TSB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TSB_CK1(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TSB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TSB_RK5(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TSB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TSB_RK4(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TSB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TSB_RK3(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TSB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TSB_RK2(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TSB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_TSB_RK1(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_TSB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transSymm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_THB_IK5(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_THB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_THB_IK4(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_THB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_THB_IK3(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_THB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_THB_IK2(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_THB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_THB_IK1(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_THB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_THB_CK5(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_THB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_THB_CK4(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_THB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_THB_CK3(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_THB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_THB_CK2(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_THB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_THB_CK1(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_THB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_THB_RK5(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_THB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_THB_RK4(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_THB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_THB_RK3(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_THB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_THB_RK2(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_THB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TNA_THB_RK1(matA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TNA_THB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_THB_IK5(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_THB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_THB_IK4(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_THB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_THB_IK3(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_THB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_THB_IK2(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_THB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_THB_IK1(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_THB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_THB_CK5(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_THB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_THB_CK4(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_THB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_THB_CK3(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_THB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_THB_CK2(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_THB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_THB_CK1(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_THB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_THB_RK5(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_THB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_THB_RK4(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_THB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_THB_RK3(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_THB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_THB_RK2(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_THB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_TSA_THB_RK1(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_TSA_THB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transSymm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_THB_IK5(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_THB_IK5
#endif
        use pm_kind, only: IKC => IK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_THB_IK4(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_THB_IK4
#endif
        use pm_kind, only: IKC => IK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_THB_IK3(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_THB_IK3
#endif
        use pm_kind, only: IKC => IK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_THB_IK2(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_THB_IK2
#endif
        use pm_kind, only: IKC => IK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if IK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_THB_IK1(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_THB_IK1
#endif
        use pm_kind, only: IKC => IK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        integer(IKC)            , intent(in)                    :: alpha, beta
        integer(IKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        integer(IKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        integer(IKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_THB_CK5(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_THB_CK5
#endif
        use pm_kind, only: CKC => CK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_THB_CK4(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_THB_CK4
#endif
        use pm_kind, only: CKC => CK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_THB_CK3(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_THB_CK3
#endif
        use pm_kind, only: CKC => CK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_THB_CK2(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_THB_CK2
#endif
        use pm_kind, only: CKC => CK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_THB_CK1(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_THB_CK1
#endif
        use pm_kind, only: CKC => CK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        complex(CKC)            , intent(in)                    :: alpha, beta
        complex(CKC)            , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        complex(CKC)            , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        complex(CKC)            , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_THB_RK5(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_THB_RK5
#endif
        use pm_kind, only: RKC => RK5
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_THB_RK4(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_THB_RK4
#endif
        use pm_kind, only: RKC => RK4
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_THB_RK3(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_THB_RK3
#endif
        use pm_kind, only: RKC => RK3
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_THB_RK2(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_THB_RK2
#endif
        use pm_kind, only: RKC => RK2
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine gemm_EXP_SFA_SFB_THA_THB_RK1(matA, operationA, matB, operationB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: gemm_EXP_SFA_SFB_THA_THB_RK1
#endif
        use pm_kind, only: RKC => RK1
        integer(IK)             , intent(in)                    :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC
        real(RKC)               , intent(in)                    :: alpha, beta
        real(RKC)               , intent(in)    , contiguous    :: matA(1 - roffA :, 1 - coffA :)
        real(RKC)               , intent(in)    , contiguous    :: matB(1 - roffB :, 1 - coffB :)
        real(RKC)               , intent(inout) , contiguous    :: matC(1 - roffC :, 1 - coffC :)
        type(transHerm_type)    , intent(in)                    :: operationA
        type(transHerm_type)    , intent(in)                    :: operationB
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_matrixMulAdd ! LCOV_EXCL_LINE