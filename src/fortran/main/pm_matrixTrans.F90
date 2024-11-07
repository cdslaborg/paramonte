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
!>  This module contains abstract and concrete derived types and procedures related to various common matrix
!>  transposition operations for which there is a corresponding matrix class defined in [pm_matrixClass](@ref pm_matrixClass).<br>
!>
!>  \details
!>  There are a myriad of operations that can be applied to modify matrices, for example,
!>  matrix addition, scalar multiplication, transposition, matrix multiplication, row operations, and submatrix.<br>
!>  The entities and derived types of this module currently focus on operations that are required for compile-time
!>  resolution of procedures within the generic interfaces of the ParaMonte library for Linear Algebra operations.<br>
!>  Such procedures frequently need to work on either an Identical, Inverse, Symmetric transpose, Hermitian transpose,
!>  Orthogonal Transpose, Unitary transpose, or other forms of some of their input matrix arguments.<br>
!>
!>  Transposition operation
!>  =======================
!>
!>  In linear algebra, the **transpose of a matrix** is an operator which flips a matrix over its diagonal.<br>
!>  It switches the row and column indices of a given matrix \f$A\f$ by producing another (transposed) matrix.<br>
!>  The transpose of a matrix was introduced in **1858** by the British mathematician Arthur Cayley.<br>
!>
!>  **Notation**<br>
!>  The transpose of a given matrix \f$A\f$ is frequently denoted by \f$A^T\f$.<br>
!>  In the case of square matrices, \f$A^T\f$ may also denote the \f$T\f$th power of the matrix \f$A\f$.<br>
!>  To avoid a possible confusion, some use left upperscripts to denote transpose, that is, \f${}^TA\f$.<br>
!>  An advantage of this notation is that no parentheses are needed when exponents are involved as
!>  \f$({}^TA)^n = {}^T(A^n)\f$, notation \f${}^TA^n\f$ is not ambiguous.<br>
!>
!>  **Definition**<br>
!>  The transpose of a matrix \f$A\f$ may be constructed by any one of the following methods,
!>  <ol>
!>      <li>    Reflect \f$A\f$ over its main diagonal (which runs from top-left to bottom-right) to obtain \f$A^T\f$.
!>      <li>    Write the rows of \f$A\f$ as the columns of \f$A^T\f$.
!>      <li>    Write the columns of \f$A\f$ as the rows of \f$A^T\f$.
!>  </ol>
!>  Formally, the \f$i\f$-th row, \f$j\f$-th column element of \f$A^T\f$ is the \f$j\f$-th row, \f$i\f$-th column element of \f$A\f$,
!>  \f{equation}{
!>      \left[\mathbf{A}^{\up{T}}\right]_{ij} = \left[\mathbf{A} \right]_{ji} ~.
!>  \f}
!>  If \f$A\f$ is an \f$m \times n\f$ matrix, then \f$A^T\f$ is an \f$n \times m\f$ matrix.<br>
!>
!>  **Matrix transposition types and correspondence to [matrix classes](@ref pm_matrixClass):**<br>
!>  <ol>
!>      <li>    **[Symmetric transposition](@ref pm_matrixTrans::transSymm_type)**<br>
!>              A square matrix whose transpose is equal to itself is called a Symmetric matrix.<br>
!>              In other words, \f$A\f$ is Symmetric if \f$\mathbf{A}^{\up{T}} = \mathbf{A}\f$.<br>
!>              The corresponding transposition is called [Symmetric](@ref pm_matrixTrans::transSymm_type)
!>              denoted by the operator \f$\cdot^{\up{T}}\f$.<br>
!>      <li>    **[Skew-Symmetric transposition](@ref pm_matrixTrans::transSymmSkew_type)**<br>
!>              A square matrix whose transpose is equal to its negative is called a Skew-Symmetric matrix
!>              In other words, \f$A\f$ is Skew-Symmetric if \f$\mathbf{A}^{\up{T}} = -\mathbf{A}\f$.<br>
!>              The corresponding transposition is called [Skew-Symmetric](@ref pm_matrixTrans::transSymmSkew_type)
!>              denoted by the operator \f$-\cdot^{\up{T}}\f$.<br>
!>      <li>    **[Hermitian transposition](@ref pm_matrixTrans::transHerm_type)**<br>
!>              A square **complex matrix** whose transpose equals the matrix with every entry replaced
!>              by its complex conjugate (\f$\overline{\cdot}\f$) is called a Hermitian matrix.<br>
!>              In other words, \f$A\f$ is Hermitian if \f$\mathbf{A}^{\up{T}} = \overline{\mathbf{A}}\f$.<br>
!>              The corresponding transposition is called [Hermitian (Conjugate)](@ref pm_matrixTrans::transHerm_type)
!>              denoted by the operator \f$\cdot^{\up{H}}\f$.<br>
!>      <li>    **[Skew-Hermitian transposition](@ref pm_matrixTrans::transHermSkew_type)**<br>
!>              A square **complex matrix** whose transpose is equal to the negation of its complex conjugate is called a Skew-Hermitian matrix.<br>
!>              In other words, \f$A\f$ is Skew-Hermitian if \f$\mathbf{A}^{\up{T}} = -{\overline{\mathbf{A}}}\f$.<br>
!>              The corresponding transposition is called [Skew-Hermitian (Skew-Conjugate)](@ref pm_matrixTrans::transHermSkew_type)
!>              denoted by the operator \f$-\cdot^{\up{H}}\f$.<br>
!>  </ol>
!>
!>  Inversion-Transposition operation
!>  =================================
!>
!>  There are also special matrix operations that mix **inversion** with Symmetric
!>  and Hermitian each having a corresponding [matrix classes](@ref pm_matrixClass):**<br>
!>  <ol>
!>      <li>    **[Orthogonal Transposition](@ref pm_matrixTrans::transOrth_type)**<br>
!>              A square matrix whose transpose is equal to its inverse is called an Orthogonal matrix.<br>
!>              In other words, \f$A\f$ is Orthogonal if \f$\mathbf{A}^{\up{T}} = \mathbf{A}^{-1}\f$.<br>
!>              The corresponding transposition is called [Orthogonal](@ref pm_matrixTrans::transOrth_type) denoted by the operator \f$\cdot^{\up{-T}}\f$.<br>
!>      <li>    **[Unitary transposition](@ref pm_matrixTrans::transUnit_type)**<br>
!>              A square complex matrix whose transpose is equal to its conjugate inverse is called a Unitary matrix.<br>
!>              In other words, \f$A\f$ is Unitary if \f$\mathbf{A}^{\up{H}} = {\overline{\mathbf{A}^{-1}}}\f$.<br>
!>              The corresponding transposition is called [Unitary](@ref pm_matrixTrans::transUnit_type) denoted by the operator \f$\cdot^{\up{-H}}\f$.<br>
!>  </ol>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_matrixTrans

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

    character(*,SK), parameter :: MODULE_NAME = "@pm_matrixTrans"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require different forms of transposition (Symmetric, Hermitian, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [trans](@ref pm_matrixTrans::trans)<br>
    !>  [transSymm](@ref pm_matrixTrans::transSymm)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [transOrth](@ref pm_matrixTrans::transOrth)<br>
    !>  [transUnit](@ref pm_matrixTrans::transUnit)<br>
    !>  [transSymmSkew](@ref pm_matrixTrans::transSymmSkew)<br>
    !>  [transHermSkew](@ref pm_matrixTrans::transHermSkew)<br>
    !>  [transSymm_type](@ref pm_matrixTrans::transSymm_type)<br>
    !>  [transHerm_type](@ref pm_matrixTrans::transHerm_type)<br>
    !>  [transOrth_type](@ref pm_matrixTrans::transOrth_type)<br>
    !>  [transUnit_type](@ref pm_matrixTrans::transUnit_type)<br>
    !>  [transSymmSkew_type](@ref pm_matrixTrans::transSymmSkew_type)<br>
    !>  [transHermSkew_type](@ref pm_matrixTrans::transHermSkew_type)<br>
    !>  [trans_type](@ref pm_matrixTrans::trans_type)<br>
    !>
    !>  \final{trans_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: trans_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [trans_type](@ref pm_matrixTrans::trans_type) that is exclusively used
    !>  to request no transpose of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [trans](@ref pm_matrixTrans::trans)<br>
    !>  [transSymm](@ref pm_matrixTrans::transSymm)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [transOrth](@ref pm_matrixTrans::transOrth)<br>
    !>  [transUnit](@ref pm_matrixTrans::transUnit)<br>
    !>  [transSymmSkew](@ref pm_matrixTrans::transSymmSkew)<br>
    !>  [transHermSkew](@ref pm_matrixTrans::transHermSkew)<br>
    !>  [transSymm_type](@ref pm_matrixTrans::transSymm_type)<br>
    !>  [transHerm_type](@ref pm_matrixTrans::transHerm_type)<br>
    !>  [transOrth_type](@ref pm_matrixTrans::transOrth_type)<br>
    !>  [transUnit_type](@ref pm_matrixTrans::transUnit_type)<br>
    !>  [transSymmSkew_type](@ref pm_matrixTrans::transSymmSkew_type)<br>
    !>  [transHermSkew_type](@ref pm_matrixTrans::transHermSkew_type)<br>
    !>  [trans_type](@ref pm_matrixTrans::trans_type)<br>
    !>
    !>  \final{trans}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(trans_type), parameter :: trans = trans_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: trans
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request Symmetric transpose (\f$\cdot^T\f$) of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [transSymm](@ref pm_matrixTrans::transSymm)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [trans](@ref pm_matrixTrans::trans)<br>
    !>  [transSymm](@ref pm_matrixTrans::transSymm)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [transOrth](@ref pm_matrixTrans::transOrth)<br>
    !>  [transUnit](@ref pm_matrixTrans::transUnit)<br>
    !>  [transSymmSkew](@ref pm_matrixTrans::transSymmSkew)<br>
    !>  [transHermSkew](@ref pm_matrixTrans::transHermSkew)<br>
    !>  [transSymm_type](@ref pm_matrixTrans::transSymm_type)<br>
    !>  [transHerm_type](@ref pm_matrixTrans::transHerm_type)<br>
    !>  [transOrth_type](@ref pm_matrixTrans::transOrth_type)<br>
    !>  [transUnit_type](@ref pm_matrixTrans::transUnit_type)<br>
    !>  [transSymmSkew_type](@ref pm_matrixTrans::transSymmSkew_type)<br>
    !>  [transHermSkew_type](@ref pm_matrixTrans::transHermSkew_type)<br>
    !>  [trans_type](@ref pm_matrixTrans::trans_type)<br>
    !>
    !>  \final{transSymm_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(trans_type) :: transSymm_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [transSymm_type](@ref pm_matrixTrans::transSymm_type) that is exclusively used
    !>  to request Symmetric transpose (\f$\cdot^T\f$) of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [trans](@ref pm_matrixTrans::trans)<br>
    !>  [transSymm](@ref pm_matrixTrans::transSymm)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [transOrth](@ref pm_matrixTrans::transOrth)<br>
    !>  [transUnit](@ref pm_matrixTrans::transUnit)<br>
    !>  [transSymmSkew](@ref pm_matrixTrans::transSymmSkew)<br>
    !>  [transHermSkew](@ref pm_matrixTrans::transHermSkew)<br>
    !>  [transSymm_type](@ref pm_matrixTrans::transSymm_type)<br>
    !>  [transHerm_type](@ref pm_matrixTrans::transHerm_type)<br>
    !>  [transOrth_type](@ref pm_matrixTrans::transOrth_type)<br>
    !>  [transUnit_type](@ref pm_matrixTrans::transUnit_type)<br>
    !>  [transSymmSkew_type](@ref pm_matrixTrans::transSymmSkew_type)<br>
    !>  [transHermSkew_type](@ref pm_matrixTrans::transHermSkew_type)<br>
    !>  [trans_type](@ref pm_matrixTrans::trans_type)<br>
    !>
    !>  \final{transSymm}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(transSymm_type), parameter :: transSymm = transSymm_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: transSymm
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request Hermitian (conjugate) transpose (\f$\cdot^H\f$) of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [transHerm](@ref pm_matrixTrans::transHerm)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [trans](@ref pm_matrixTrans::trans)<br>
    !>  [transSymm](@ref pm_matrixTrans::transSymm)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [transOrth](@ref pm_matrixTrans::transOrth)<br>
    !>  [transUnit](@ref pm_matrixTrans::transUnit)<br>
    !>  [transSymmSkew](@ref pm_matrixTrans::transSymmSkew)<br>
    !>  [transHermSkew](@ref pm_matrixTrans::transHermSkew)<br>
    !>  [transSymm_type](@ref pm_matrixTrans::transSymm_type)<br>
    !>  [transHerm_type](@ref pm_matrixTrans::transHerm_type)<br>
    !>  [transOrth_type](@ref pm_matrixTrans::transOrth_type)<br>
    !>  [transUnit_type](@ref pm_matrixTrans::transUnit_type)<br>
    !>  [transSymmSkew_type](@ref pm_matrixTrans::transSymmSkew_type)<br>
    !>  [transHermSkew_type](@ref pm_matrixTrans::transHermSkew_type)<br>
    !>  [trans_type](@ref pm_matrixTrans::trans_type)<br>
    !>
    !>  \final{transHerm_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(trans_type) :: transHerm_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [transHerm_type](@ref pm_matrixTrans::transHerm_type) that is exclusively used
    !>  to request Hermitian (conjugate) transpose (\f$\cdot^T\f$) of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [trans](@ref pm_matrixTrans::trans)<br>
    !>  [transSymm](@ref pm_matrixTrans::transSymm)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [transOrth](@ref pm_matrixTrans::transOrth)<br>
    !>  [transUnit](@ref pm_matrixTrans::transUnit)<br>
    !>  [transSymmSkew](@ref pm_matrixTrans::transSymmSkew)<br>
    !>  [transHermSkew](@ref pm_matrixTrans::transHermSkew)<br>
    !>  [transSymm_type](@ref pm_matrixTrans::transSymm_type)<br>
    !>  [transHerm_type](@ref pm_matrixTrans::transHerm_type)<br>
    !>  [transOrth_type](@ref pm_matrixTrans::transOrth_type)<br>
    !>  [transUnit_type](@ref pm_matrixTrans::transUnit_type)<br>
    !>  [transSymmSkew_type](@ref pm_matrixTrans::transSymmSkew_type)<br>
    !>  [transHermSkew_type](@ref pm_matrixTrans::transHermSkew_type)<br>
    !>  [trans_type](@ref pm_matrixTrans::trans_type)<br>
    !>  [nothing](@ref pm_array::nothing)<br>
    !>
    !>  \final{transHerm}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(transHerm_type), parameter :: transHerm = transHerm_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: transHerm
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request Orthogonal Transpose (\f$\cdot^-T\f$) of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [transOrth](@ref pm_matrixTrans::transOrth)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [trans](@ref pm_matrixTrans::trans)<br>
    !>  [transSymm](@ref pm_matrixTrans::transSymm)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [transOrth](@ref pm_matrixTrans::transOrth)<br>
    !>  [transUnit](@ref pm_matrixTrans::transUnit)<br>
    !>  [transSymmSkew](@ref pm_matrixTrans::transSymmSkew)<br>
    !>  [transHermSkew](@ref pm_matrixTrans::transHermSkew)<br>
    !>  [transSymm_type](@ref pm_matrixTrans::transSymm_type)<br>
    !>  [transHerm_type](@ref pm_matrixTrans::transHerm_type)<br>
    !>  [transOrth_type](@ref pm_matrixTrans::transOrth_type)<br>
    !>  [transUnit_type](@ref pm_matrixTrans::transUnit_type)<br>
    !>  [transSymmSkew_type](@ref pm_matrixTrans::transSymmSkew_type)<br>
    !>  [transHermSkew_type](@ref pm_matrixTrans::transHermSkew_type)<br>
    !>  [trans_type](@ref pm_matrixTrans::trans_type)<br>
    !>
    !>  \final{transOrth_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(trans_type) :: transOrth_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [transOrth_type](@ref pm_matrixTrans::transOrth_type) that is exclusively used
    !>  to request Orthogonal Transpose (\f$\cdot^-T\f$) of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [trans](@ref pm_matrixTrans::trans)<br>
    !>  [transSymm](@ref pm_matrixTrans::transSymm)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [transOrth](@ref pm_matrixTrans::transOrth)<br>
    !>  [transUnit](@ref pm_matrixTrans::transUnit)<br>
    !>  [transSymmSkew](@ref pm_matrixTrans::transSymmSkew)<br>
    !>  [transHermSkew](@ref pm_matrixTrans::transHermSkew)<br>
    !>  [transSymm_type](@ref pm_matrixTrans::transSymm_type)<br>
    !>  [transHerm_type](@ref pm_matrixTrans::transHerm_type)<br>
    !>  [transOrth_type](@ref pm_matrixTrans::transOrth_type)<br>
    !>  [transUnit_type](@ref pm_matrixTrans::transUnit_type)<br>
    !>  [transSymmSkew_type](@ref pm_matrixTrans::transSymmSkew_type)<br>
    !>  [transHermSkew_type](@ref pm_matrixTrans::transHermSkew_type)<br>
    !>  [trans_type](@ref pm_matrixTrans::trans_type)<br>
    !>
    !>  \final{transOrth}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(transOrth_type), parameter :: transOrth = transOrth_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: transOrth
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request Unitary Transpose (\f$\cdot^{-H}\f$) of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [transUnit](@ref pm_matrixTrans::transUnit)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [trans](@ref pm_matrixTrans::trans)<br>
    !>  [transSymm](@ref pm_matrixTrans::transSymm)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [transOrth](@ref pm_matrixTrans::transOrth)<br>
    !>  [transUnit](@ref pm_matrixTrans::transUnit)<br>
    !>  [transSymmSkew](@ref pm_matrixTrans::transSymmSkew)<br>
    !>  [transHermSkew](@ref pm_matrixTrans::transHermSkew)<br>
    !>  [transSymm_type](@ref pm_matrixTrans::transSymm_type)<br>
    !>  [transHerm_type](@ref pm_matrixTrans::transHerm_type)<br>
    !>  [transOrth_type](@ref pm_matrixTrans::transOrth_type)<br>
    !>  [transUnit_type](@ref pm_matrixTrans::transUnit_type)<br>
    !>  [transSymmSkew_type](@ref pm_matrixTrans::transSymmSkew_type)<br>
    !>  [transHermSkew_type](@ref pm_matrixTrans::transHermSkew_type)<br>
    !>  [trans_type](@ref pm_matrixTrans::trans_type)<br>
    !>
    !>  \final{transUnit_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(trans_type) :: transUnit_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [transUnit_type](@ref pm_matrixTrans::transUnit_type) that is exclusively used
    !>  to request Unitary Transpose (\f$\cdot^{-H}\f$) of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [trans](@ref pm_matrixTrans::trans)<br>
    !>  [transSymm](@ref pm_matrixTrans::transSymm)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [transOrth](@ref pm_matrixTrans::transOrth)<br>
    !>  [transUnit](@ref pm_matrixTrans::transUnit)<br>
    !>  [transSymmSkew](@ref pm_matrixTrans::transSymmSkew)<br>
    !>  [transHermSkew](@ref pm_matrixTrans::transHermSkew)<br>
    !>  [transSymm_type](@ref pm_matrixTrans::transSymm_type)<br>
    !>  [transHerm_type](@ref pm_matrixTrans::transHerm_type)<br>
    !>  [transOrth_type](@ref pm_matrixTrans::transOrth_type)<br>
    !>  [transUnit_type](@ref pm_matrixTrans::transUnit_type)<br>
    !>  [transSymmSkew_type](@ref pm_matrixTrans::transSymmSkew_type)<br>
    !>  [transHermSkew_type](@ref pm_matrixTrans::transHermSkew_type)<br>
    !>  [trans_type](@ref pm_matrixTrans::trans_type)<br>
    !>
    !>  \final{transUnit}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(transUnit_type), parameter :: transUnit = transUnit_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: transUnit
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request Skew-Symmetric transpose (\f$-\cdot^T\f$) of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [transSymmSkew](@ref pm_matrixTrans::transSymmSkew)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [trans](@ref pm_matrixTrans::trans)<br>
    !>  [transSymm](@ref pm_matrixTrans::transSymm)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [transOrth](@ref pm_matrixTrans::transOrth)<br>
    !>  [transUnit](@ref pm_matrixTrans::transUnit)<br>
    !>  [transSymmSkew](@ref pm_matrixTrans::transSymmSkew)<br>
    !>  [transHermSkew](@ref pm_matrixTrans::transHermSkew)<br>
    !>  [transSymm_type](@ref pm_matrixTrans::transSymm_type)<br>
    !>  [transHerm_type](@ref pm_matrixTrans::transHerm_type)<br>
    !>  [transOrth_type](@ref pm_matrixTrans::transOrth_type)<br>
    !>  [transUnit_type](@ref pm_matrixTrans::transUnit_type)<br>
    !>  [transSymmSkew_type](@ref pm_matrixTrans::transSymmSkew_type)<br>
    !>  [transHermSkew_type](@ref pm_matrixTrans::transHermSkew_type)<br>
    !>  [trans_type](@ref pm_matrixTrans::trans_type)<br>
    !>
    !>  \final{transSymmSkew_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(trans_type) :: transSymmSkew_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [transSymmSkew_type](@ref pm_matrixTrans::transSymmSkew_type) that is exclusively used
    !>  to request Skew-Symmetric transpose (\f$\cdot^T\f$) of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [trans](@ref pm_matrixTrans::trans)<br>
    !>  [transSymm](@ref pm_matrixTrans::transSymm)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [transOrth](@ref pm_matrixTrans::transOrth)<br>
    !>  [transUnit](@ref pm_matrixTrans::transUnit)<br>
    !>  [transSymmSkew](@ref pm_matrixTrans::transSymmSkew)<br>
    !>  [transHermSkew](@ref pm_matrixTrans::transHermSkew)<br>
    !>  [transSymm_type](@ref pm_matrixTrans::transSymm_type)<br>
    !>  [transHerm_type](@ref pm_matrixTrans::transHerm_type)<br>
    !>  [transOrth_type](@ref pm_matrixTrans::transOrth_type)<br>
    !>  [transUnit_type](@ref pm_matrixTrans::transUnit_type)<br>
    !>  [transSymmSkew_type](@ref pm_matrixTrans::transSymmSkew_type)<br>
    !>  [transHermSkew_type](@ref pm_matrixTrans::transHermSkew_type)<br>
    !>  [trans_type](@ref pm_matrixTrans::trans_type)<br>
    !>
    !>  \final{transSymmSkew}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(transSymmSkew_type), parameter :: transSymmSkew = transSymmSkew_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: transSymmSkew
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to request Skew-Hermitian transpose (\f$-\cdot^H\f$) of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [transHermSkew](@ref pm_matrixTrans::transHermSkew)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [trans](@ref pm_matrixTrans::trans)<br>
    !>  [transSymm](@ref pm_matrixTrans::transSymm)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [transOrth](@ref pm_matrixTrans::transOrth)<br>
    !>  [transUnit](@ref pm_matrixTrans::transUnit)<br>
    !>  [transSymmSkew](@ref pm_matrixTrans::transSymmSkew)<br>
    !>  [transHermSkew](@ref pm_matrixTrans::transHermSkew)<br>
    !>  [transSymm_type](@ref pm_matrixTrans::transSymm_type)<br>
    !>  [transHerm_type](@ref pm_matrixTrans::transHerm_type)<br>
    !>  [transOrth_type](@ref pm_matrixTrans::transOrth_type)<br>
    !>  [transUnit_type](@ref pm_matrixTrans::transUnit_type)<br>
    !>  [transSymmSkew_type](@ref pm_matrixTrans::transSymmSkew_type)<br>
    !>  [transHermSkew_type](@ref pm_matrixTrans::transHermSkew_type)<br>
    !>  [trans_type](@ref pm_matrixTrans::trans_type)<br>
    !>
    !>  \final{transHermSkew_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(trans_type) :: transHermSkew_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [transHermSkew_type](@ref pm_matrixTrans::transHermSkew_type) that is exclusively used
    !>  to request Skew-Hermitian transpose (\f$\cdot^T\f$) of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [trans](@ref pm_matrixTrans::trans)<br>
    !>  [transSymm](@ref pm_matrixTrans::transSymm)<br>
    !>  [transHerm](@ref pm_matrixTrans::transHerm)<br>
    !>  [transOrth](@ref pm_matrixTrans::transOrth)<br>
    !>  [transUnit](@ref pm_matrixTrans::transUnit)<br>
    !>  [transSymmSkew](@ref pm_matrixTrans::transSymmSkew)<br>
    !>  [transHermSkew](@ref pm_matrixTrans::transHermSkew)<br>
    !>  [transSymm_type](@ref pm_matrixTrans::transSymm_type)<br>
    !>  [transHerm_type](@ref pm_matrixTrans::transHerm_type)<br>
    !>  [transOrth_type](@ref pm_matrixTrans::transOrth_type)<br>
    !>  [transUnit_type](@ref pm_matrixTrans::transUnit_type)<br>
    !>  [transSymmSkew_type](@ref pm_matrixTrans::transSymmSkew_type)<br>
    !>  [transHermSkew_type](@ref pm_matrixTrans::transHermSkew_type)<br>
    !>  [trans_type](@ref pm_matrixTrans::trans_type)<br>
    !>
    !>  \final{transHermSkew}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(transHermSkew_type), parameter :: transHermSkew = transHermSkew_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: transHermSkew
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the transpose of the input matrix of arbitrary type and kind using a cache-oblivious approach.<br>
    !>
    !>  \details
    !>  In computing, a **cache-oblivious** (or **cache-transcendent**) algorithm is a method designed to take advantage
    !>  of a [processor cache](https://en.wikipedia.org/wiki/CPU_cache) without having the size of the cache
    !>  (or the length of the cache lines, etc.) as an explicit parameter.<br>
    !>  An optimal cache-oblivious algorithm is a cache-oblivious algorithm that uses the cache optimally.<br>
    !>  Thus, a cache-oblivious algorithm is designed to perform well, without modification,
    !>  on multiple machines with different cache sizes, or for a memory hierarchy
    !>  with different levels of cache having different sizes.<br>
    !>  Cache-oblivious algorithms are contrasted with **explicit loop tiling**,
    !>  which explicitly breaks a problem into blocks that are optimally sized for a given cache.<br>
    !>
    !>  Typically, a cache-oblivious algorithm works by a recursive divide-and-conquer algorithm,
    !>  where the problem is divided into smaller and smaller subproblems.<br>
    !>  Eventually, one reaches a subproblem size that fits into the cache, regardless of the cache size.<br>
    !>  For example, an optimal cache-oblivious matrix multiplication is obtained by recursively dividing
    !>  each matrix into four sub-matrices to be multiplied, multiplying the submatrices in a depth-first fashion.<br>
    !>  In tuning for a specific machine, one may use a hybrid algorithm which uses loop tiling tuned for
    !>  the specific cache sizes at the bottom level but otherwise uses the cache-oblivious algorithm.<br>
    !>
    !>  \param[inout]   source      :   The input/output matrix (of rank `2`) of either<br>
    !>                                  <ol>
    !>                                      <li>    type `character` of kind \SKALL or, <br>
    !>                                      <li>    type `integer` of kind \IKALL or, <br>
    !>                                      <li>    type `logical` of kind \LKALL or, <br>
    !>                                      <li>    type `complex` of kind \CKALL or, <br>
    !>                                      <li>    type `real` of kind \RKALL or, <br>
    !>                                  </ol>
    !>                                  whose contents will be Symmetric or Hermitian transposed.<br>
    !>                                  <ol>
    !>                                      <li>    If the output matrix argument `destin` is missing, the result of transposition will be written to `source`.<br>
    !>                                              This is possible only if the input `source` is a square matrix.<br>
    !>                                      <li>    If the output matrix argument `destin` is present, the result of transposition will be written to `destin`.<br>
    !>                                              As such, the input `source` has `intent(in)` will not be modified by the algorithm.<br>
    !>                                  </ol>
    !>  \param[out]     destin      :   The output matrix of the same type and kind, but transposed shape of `source` containing the transposition.
    !>                                  (**optional**. If missing, the transposition result will be written to the input `source`, in which case, `source` must be square.)
    !>  \param[in]      bsize       :   The input positive scalar integer of default kind \IK representing the minimum submatrix size.<br>
    !>                                  Any input `source` or subset of it whose size along both dimensions is below `bsize`
    !>                                  will be transposed via the default Fortran `transpose()` procedure.<br>
    !>                                  (**optional**. default = `32`)
    !>  \param[in]      operation   :   The input scalar that can be,
    !>                                  <ol>
    !>                                      <li>    the constant [transHerm](@ref pm_matrixTrans::transHerm) exclusively when `source` is of type `complex` of kind \CKALL.
    !>                                              implying that a Hermitian transpose of the specified subset of `source` is to be computed and stored.
    !>                                  </ol>
    !>                                  This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.<br>
    !>                                  (**optional**. If missing, the Symmetric transposition will be returned for complex matrices.)
    !>
    !>  \interface{setMatTrans}
    !>  \code{.F90}
    !>
    !>      use pm_matrixTrans, only: setMatTrans, transHerm
    !>
    !>      call setMatTrans(source(1:ndim,1:ndim))
    !>      call setMatTrans(source(1:ndim,1:ndim), bsize)
    !>      call setMatTrans(source(1:ndim,1:ndim), operation)
    !>      call setMatTrans(source(1:ndim,1:ndim), operation, bsize)
    !>      call setMatTrans(source(1:nrow,1:ncol), destin(1:ncol,1:nrow))
    !>      call setMatTrans(source(1:nrow,1:ncol), destin(1:ncol,1:nrow), bsize)
    !>      call setMatTrans(source(1:nrow,1:ncol), destin(1:ncol,1:nrow), operation)
    !>      call setMatTrans(source(1:nrow,1:ncol), destin(1:ncol,1:nrow), operation, bsize)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `0 < bsize` must hold for the corresponding input arguments.<br>
    !>  The condition `size(source, 1) == size(source, 2)` must hold when the output argument `destin` is missing.<br>
    !>  The condition `size(source, 1) == size(destin, 2) .and. size(source, 2) == size(destin, 1)` must hold for the corresponding input arguments.<br>
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  Based on some relevant benchmarks performed, the `contiguous` attribute for the input arguments does not appear
    !>  to have any noticeable impact on the performance of the algorithm in the release optimized compilation modes.<br>
    !>
    !>  \see
    !>  [pm_matrixCopy](@ref pm_matrixCopy)<br>
    !>
    !>  \example{setMatTrans}
    !>  \include{lineno} example/pm_matrixTrans/setMatTrans/main.F90
    !>  \compilef{setMatTrans}
    !>  \output{setMatTrans}
    !>  \include{lineno} example/pm_matrixTrans/setMatTrans/main.out.F90
    !>
    !>  \benchmarks
    !>
    !>  \benchmark{setMatTrans_vs_transpose, The runtime performance of [setMatTrans](@ref pm_matrixTrans::setMatTrans) vs. Fortran intrinsic `transpose()`.}
    !>  \include{lineno} benchmark/pm_matrixTrans/setMatTrans_vs_transpose/main.F90
    !>  \compilefb{setMatTrans_vs_transpose}
    !>  \postprocb{setMatTrans_vs_transpose}
    !>  \include{lineno} benchmark/pm_matrixTrans/setMatTrans_vs_transpose/main.py
    !>  \visb{setMatTrans_vs_transpose}
    !>  \image html benchmark/pm_matrixTrans/setMatTrans_vs_transpose/benchmark.setMatTrans_vs_transpose.runtime.png width=1000
    !>  \image html benchmark/pm_matrixTrans/setMatTrans_vs_transpose/benchmark.setMatTrans_vs_transpose.runtime.ratio.png width=1000
    !>  \moralb{setMatTrans_vs_transpose}
    !>      -#  The procedures under the generic interface [setMatTrans](@ref pm_matrixTrans::setMatTrans)
    !>          use a cache-oblivious approach to matrix Symmetric transposition.<br>
    !>          As such, they are particularly efficient and cache-friendly for large matrices.<br>
    !>          As such, the generic interface [setMatTrans](@ref pm_matrixTrans::setMatTrans) can be significantly
    !>          faster than the Fortran intrinsic `transpose()`, depending on the Fortran compiler used.<br>
    !>          This is particularly true for large-order matrices.
    !>
    !>  \benchmark{setMatTrans_vs_transHerm, The runtime performance of [setMatTrans](@ref pm_matrixTrans::setMatTrans) vs. Fortran intrinsic `transpose()`.}
    !>  \include{lineno} benchmark/pm_matrixTrans/setMatTrans_vs_transHerm/main.F90
    !>  \compilefb{setMatTrans_vs_transHerm}
    !>  \postprocb{setMatTrans_vs_transHerm}
    !>  \include{lineno} benchmark/pm_matrixTrans/setMatTrans_vs_transHerm/main.py
    !>  \visb{setMatTrans_vs_transHerm}
    !>  \image html benchmark/pm_matrixTrans/setMatTrans_vs_transHerm/benchmark.setMatTrans_vs_transHerm.runtime.png width=1000
    !>  \image html benchmark/pm_matrixTrans/setMatTrans_vs_transHerm/benchmark.setMatTrans_vs_transHerm.runtime.ratio.png width=1000
    !>  \moralb{setMatTrans_vs_transHerm}
    !>      -#  The procedures under the generic interface [setMatTrans](@ref pm_matrixTrans::setMatTrans)
    !>          use a cache-oblivious approach to matrix Hermitian transposition.<br>
    !>          As such, they are particularly efficient and cache-friendly for large matrices.<br>
    !>          As such, the generic interface [setMatTrans](@ref pm_matrixTrans::setMatTrans) can be significantly
    !>          faster than the Fortran intrinsic `transpose()`, depending on the Fortran compiler used.<br>
    !>          This is particularly true for large-order matrices.
    !>
    !>  \benchmark{setMatTransBlock, The runtime performance of [setMatTrans](@ref pm_matrixTrans::setMatTrans) vs. Fortran intrinsic `transpose()`.}
    !>  \include{lineno} benchmark/pm_matrixTrans/setMatTransBlock/main.F90
    !>  \compilefb{setMatTransBlock}
    !>  \postprocb{setMatTransBlock}
    !>  \include{lineno} benchmark/pm_matrixTrans/setMatTransBlock/main.py
    !>  \visb{setMatTransBlock}
    !>  \image html benchmark/pm_matrixTrans/setMatTransBlock/benchmark.setMatTransBlock.runtime.png width=1000
    !>  \image html benchmark/pm_matrixTrans/setMatTransBlock/benchmark.setMatTransBlock.runtime.ratio.png width=1000
    !>  \moralb{setMatTransBlock}
    !>      -#  The procedures under the generic interface [setMatTrans](@ref pm_matrixTrans::setMatTrans)
    !>          use a cache-oblivious approach to matrix Hermitian transposition.<br>
    !>          As such, they are particularly efficient and cache-friendly for large matrices.<br>
    !>          However, despite its name and goals, the cache-oblivious algorithm is not entirely
    !>          independent of the cache size (and hence the minimum block size) as evidenced here.<br>
    !>
    !>  \test
    !>  [test_pm_matrixTrans](@ref test_pm_matrixTrans)
    !>
    !>  \final{setMatTrans}
    !>
    !>  \todo
    !>  \pmed
    !>  The performance of this algorithm could be possibly improved by converting the recursive procedure calls within the implementation to do-loops.<br>
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    !   Fixed block size Symmetric transpose.
    interface setMatTrans

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_SK5(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_SK5
#endif
        use pm_kind, only: LK, SKG => SK5
        character(*,SKG)    , intent(inout) :: source(:,:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_SK4(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_SK4
#endif
        use pm_kind, only: LK, SKG => SK4
        character(*,SKG)    , intent(inout) :: source(:,:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_SK3(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_SK3
#endif
        use pm_kind, only: LK, SKG => SK3
        character(*,SKG)    , intent(inout) :: source(:,:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_SK2(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_SK2
#endif
        use pm_kind, only: LK, SKG => SK2
        character(*,SKG)    , intent(inout) :: source(:,:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_SK1(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_SK1
#endif
        use pm_kind, only: LK, SKG => SK1
        character(*,SKG)    , intent(inout) :: source(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_IK5(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_IK5
#endif
        use pm_kind, only: LK, IKG => IK5
        integer(IKG)        , intent(inout) :: source(:,:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_IK4(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_IK4
#endif
        use pm_kind, only: LK, IKG => IK4
        integer(IKG)        , intent(inout) :: source(:,:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_IK3(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_IK3
#endif
        use pm_kind, only: LK, IKG => IK3
        integer(IKG)        , intent(inout) :: source(:,:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_IK2(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_IK2
#endif
        use pm_kind, only: LK, IKG => IK2
        integer(IKG)        , intent(inout) :: source(:,:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_IK1(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_IK1
#endif
        use pm_kind, only: LK, IKG => IK1
        integer(IKG)        , intent(inout) :: source(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_LK5(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_LK5
#endif
        use pm_kind, only: LK, LKG => LK5
        logical(LKG)        , intent(inout) :: source(:,:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_LK4(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_LK4
#endif
        use pm_kind, only: LK, LKG => LK4
        logical(LKG)        , intent(inout) :: source(:,:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_LK3(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_LK3
#endif
        use pm_kind, only: LK, LKG => LK3
        logical(LKG)        , intent(inout) :: source(:,:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_LK2(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_LK2
#endif
        use pm_kind, only: LK, LKG => LK2
        logical(LKG)        , intent(inout) :: source(:,:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_LK1(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_LK1
#endif
        use pm_kind, only: LK, LKG => LK1
        logical(LKG)        , intent(inout) :: source(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_CK5(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_CK5
#endif
        use pm_kind, only: LK, CKG => CK5
        complex(CKG)        , intent(inout) :: source(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_CK4(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_CK4
#endif
        use pm_kind, only: LK, CKG => CK4
        complex(CKG)        , intent(inout) :: source(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_CK3(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_CK3
#endif
        use pm_kind, only: LK, CKG => CK3
        complex(CKG)        , intent(inout) :: source(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_CK2(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_CK2
#endif
        use pm_kind, only: LK, CKG => CK2
        complex(CKG)        , intent(inout) :: source(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_CK1(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_CK1
#endif
        use pm_kind, only: LK, CKG => CK1
        complex(CKG)        , intent(inout) :: source(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_RK5(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_RK5
#endif
        use pm_kind, only: LK, RKG => RK5
        real(RKG)           , intent(inout) :: source(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_RK4(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_RK4
#endif
        use pm_kind, only: LK, RKG => RK4
        real(RKG)           , intent(inout) :: source(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_RK3(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_RK3
#endif
        use pm_kind, only: LK, RKG => RK3
        real(RKG)           , intent(inout) :: source(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_RK2(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_RK2
#endif
        use pm_kind, only: LK, RKG => RK2
        real(RKG)           , intent(inout) :: source(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setMatTransSymmOldFix_RK1(source)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldFix_RK1
#endif
        use pm_kind, only: LK, RKG => RK1
        real(RKG)           , intent(inout) :: source(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_SK5(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_SK5
#endif
        use pm_kind, only: LK, SKG => SK5
        character(*,SKG)    , intent(in)    :: source(:,:)
        character(*,SKG)    , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if SK4_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_SK4(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_SK4
#endif
        use pm_kind, only: LK, SKG => SK4
        character(*,SKG)    , intent(in)    :: source(:,:)
        character(*,SKG)    , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if SK3_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_SK3(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_SK3
#endif
        use pm_kind, only: LK, SKG => SK3
        character(*,SKG)    , intent(in)    :: source(:,:)
        character(*,SKG)    , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if SK2_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_SK2(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_SK2
#endif
        use pm_kind, only: LK, SKG => SK2
        character(*,SKG)    , intent(in)    :: source(:,:)
        character(*,SKG)    , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if SK1_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_SK1(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_SK1
#endif
        use pm_kind, only: LK, SKG => SK1
        character(*,SKG)    , intent(in)    :: source(:,:)
        character(*,SKG)    , intent(out)   :: destin(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_IK5(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_IK5
#endif
        use pm_kind, only: LK, IKG => IK5
        integer(IKG)        , intent(in)    :: source(:,:)
        integer(IKG)        , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if IK4_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_IK4(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_IK4
#endif
        use pm_kind, only: LK, IKG => IK4
        integer(IKG)        , intent(in)    :: source(:,:)
        integer(IKG)        , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if IK3_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_IK3(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_IK3
#endif
        use pm_kind, only: LK, IKG => IK3
        integer(IKG)        , intent(in)    :: source(:,:)
        integer(IKG)        , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if IK2_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_IK2(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_IK2
#endif
        use pm_kind, only: LK, IKG => IK2
        integer(IKG)        , intent(in)    :: source(:,:)
        integer(IKG)        , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if IK1_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_IK1(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_IK1
#endif
        use pm_kind, only: LK, IKG => IK1
        integer(IKG)        , intent(in)    :: source(:,:)
        integer(IKG)        , intent(out)   :: destin(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_LK5(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_LK5
#endif
        use pm_kind, only: LK, LKG => LK5
        logical(LKG)        , intent(in)    :: source(:,:)
        logical(LKG)        , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if LK4_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_LK4(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_LK4
#endif
        use pm_kind, only: LK, LKG => LK4
        logical(LKG)        , intent(in)    :: source(:,:)
        logical(LKG)        , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if LK3_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_LK3(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_LK3
#endif
        use pm_kind, only: LK, LKG => LK3
        logical(LKG)        , intent(in)    :: source(:,:)
        logical(LKG)        , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if LK2_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_LK2(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_LK2
#endif
        use pm_kind, only: LK, LKG => LK2
        logical(LKG)        , intent(in)    :: source(:,:)
        logical(LKG)        , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if LK1_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_LK1(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_LK1
#endif
        use pm_kind, only: LK, LKG => LK1
        logical(LKG)        , intent(in)    :: source(:,:)
        logical(LKG)        , intent(out)   :: destin(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_CK5(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_CK5
#endif
        use pm_kind, only: LK, CKG => CK5
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_CK4(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_CK4
#endif
        use pm_kind, only: LK, CKG => CK4
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_CK3(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_CK3
#endif
        use pm_kind, only: LK, CKG => CK3
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_CK2(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_CK2
#endif
        use pm_kind, only: LK, CKG => CK2
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_CK1(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_CK1
#endif
        use pm_kind, only: LK, CKG => CK1
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_RK5(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_RK5
#endif
        use pm_kind, only: LK, RKG => RK5
        real(RKG)           , intent(in)    :: source(:,:)
        real(RKG)           , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_RK4(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_RK4
#endif
        use pm_kind, only: LK, RKG => RK4
        real(RKG)           , intent(in)    :: source(:,:)
        real(RKG)           , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_RK3(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_RK3
#endif
        use pm_kind, only: LK, RKG => RK3
        real(RKG)           , intent(in)    :: source(:,:)
        real(RKG)           , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_RK2(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_RK2
#endif
        use pm_kind, only: LK, RKG => RK2
        real(RKG)           , intent(in)    :: source(:,:)
        real(RKG)           , intent(out)   :: destin(:,:)
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setMatTransSymmNewFix_RK1(source, destin)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewFix_RK1
#endif
        use pm_kind, only: LK, RKG => RK1
        real(RKG)           , intent(in)    :: source(:,:)
        real(RKG)           , intent(out)   :: destin(:,:)
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    !   Fixed block size Hermitian transpose.
    !cond excluded
    interface setMatTrans

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setMatTransHermOldFix_CK5(source, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermOldFix_CK5
#endif
        use pm_kind, only: LK, CKG => CK5
        complex(CKG)        , intent(inout) :: source(:,:)
        type(transHerm_type), intent(in)    :: operation
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setMatTransHermOldFix_CK4(source, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermOldFix_CK4
#endif
        use pm_kind, only: LK, CKG => CK4
        complex(CKG)        , intent(inout) :: source(:,:)
        type(transHerm_type), intent(in)    :: operation
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setMatTransHermOldFix_CK3(source, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermOldFix_CK3
#endif
        use pm_kind, only: LK, CKG => CK3
        complex(CKG)        , intent(inout) :: source(:,:)
        type(transHerm_type), intent(in)    :: operation
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setMatTransHermOldFix_CK2(source, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermOldFix_CK2
#endif
        use pm_kind, only: LK, CKG => CK2
        complex(CKG)        , intent(inout) :: source(:,:)
        type(transHerm_type), intent(in)    :: operation
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setMatTransHermOldFix_CK1(source, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermOldFix_CK1
#endif
        use pm_kind, only: LK, CKG => CK1
        complex(CKG)        , intent(inout) :: source(:,:)
        type(transHerm_type), intent(in)    :: operation
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setMatTransHermNewFix_CK5(source, destin, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermNewFix_CK5
#endif
        use pm_kind, only: LK, CKG => CK5
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
        type(transHerm_type), intent(in)    :: operation
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setMatTransHermNewFix_CK4(source, destin, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermNewFix_CK4
#endif
        use pm_kind, only: LK, CKG => CK4
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
        type(transHerm_type), intent(in)    :: operation
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setMatTransHermNewFix_CK3(source, destin, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermNewFix_CK3
#endif
        use pm_kind, only: LK, CKG => CK3
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
        type(transHerm_type), intent(in)    :: operation
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setMatTransHermNewFix_CK2(source, destin, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermNewFix_CK2
#endif
        use pm_kind, only: LK, CKG => CK2
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
        type(transHerm_type), intent(in)    :: operation
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setMatTransHermNewFix_CK1(source, destin, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermNewFix_CK1
#endif
        use pm_kind, only: LK, CKG => CK1
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
        type(transHerm_type), intent(in)    :: operation
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface
    !endcond excluded

    !   Arbitrary block size Symmetric transpose.
    !cond excluded
    interface setMatTrans

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_SK5(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_SK5
#endif
        use pm_kind, only: LK, SKG => SK5
        character(*,SKG)    , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if SK4_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_SK4(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_SK4
#endif
        use pm_kind, only: LK, SKG => SK4
        character(*,SKG)    , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if SK3_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_SK3(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_SK3
#endif
        use pm_kind, only: LK, SKG => SK3
        character(*,SKG)    , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if SK2_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_SK2(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_SK2
#endif
        use pm_kind, only: LK, SKG => SK2
        character(*,SKG)    , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if SK1_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_SK1(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_SK1
#endif
        use pm_kind, only: LK, SKG => SK1
        character(*,SKG)    , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_IK5(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_IK5
#endif
        use pm_kind, only: LK, IKG => IK5
        integer(IKG)        , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if IK4_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_IK4(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_IK4
#endif
        use pm_kind, only: LK, IKG => IK4
        integer(IKG)        , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if IK3_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_IK3(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_IK3
#endif
        use pm_kind, only: LK, IKG => IK3
        integer(IKG)        , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if IK2_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_IK2(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_IK2
#endif
        use pm_kind, only: LK, IKG => IK2
        integer(IKG)        , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if IK1_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_IK1(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_IK1
#endif
        use pm_kind, only: LK, IKG => IK1
        integer(IKG)        , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_LK5(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_LK5
#endif
        use pm_kind, only: LK, LKG => LK5
        logical(LKG)        , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if LK4_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_LK4(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_LK4
#endif
        use pm_kind, only: LK, LKG => LK4
        logical(LKG)        , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if LK3_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_LK3(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_LK3
#endif
        use pm_kind, only: LK, LKG => LK3
        logical(LKG)        , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if LK2_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_LK2(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_LK2
#endif
        use pm_kind, only: LK, LKG => LK2
        logical(LKG)        , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if LK1_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_LK1(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_LK1
#endif
        use pm_kind, only: LK, LKG => LK1
        logical(LKG)        , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_CK5(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_CK5
#endif
        use pm_kind, only: LK, CKG => CK5
        complex(CKG)        , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_CK4(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_CK4
#endif
        use pm_kind, only: LK, CKG => CK4
        complex(CKG)        , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_CK3(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_CK3
#endif
        use pm_kind, only: LK, CKG => CK3
        complex(CKG)        , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_CK2(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_CK2
#endif
        use pm_kind, only: LK, CKG => CK2
        complex(CKG)        , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_CK1(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_CK1
#endif
        use pm_kind, only: LK, CKG => CK1
        complex(CKG)        , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_RK5(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_RK5
#endif
        use pm_kind, only: LK, RKG => RK5
        real(RKG)           , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_RK4(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_RK4
#endif
        use pm_kind, only: LK, RKG => RK4
        real(RKG)           , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_RK3(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_RK3
#endif
        use pm_kind, only: LK, RKG => RK3
        real(RKG)           , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_RK2(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_RK2
#endif
        use pm_kind, only: LK, RKG => RK2
        real(RKG)           , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setMatTransSymmOldArb_RK1(source, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmOldArb_RK1
#endif
        use pm_kind, only: LK, RKG => RK1
        real(RKG)           , intent(inout) :: source(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_SK5(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_SK5
#endif
        use pm_kind, only: LK, SKG => SK5
        character(*,SKG)    , intent(in)    :: source(:,:)
        character(*,SKG)    , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if SK4_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_SK4(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_SK4
#endif
        use pm_kind, only: LK, SKG => SK4
        character(*,SKG)    , intent(in)    :: source(:,:)
        character(*,SKG)    , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if SK3_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_SK3(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_SK3
#endif
        use pm_kind, only: LK, SKG => SK3
        character(*,SKG)    , intent(in)    :: source(:,:)
        character(*,SKG)    , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if SK2_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_SK2(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_SK2
#endif
        use pm_kind, only: LK, SKG => SK2
        character(*,SKG)    , intent(in)    :: source(:,:)
        character(*,SKG)    , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if SK1_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_SK1(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_SK1
#endif
        use pm_kind, only: LK, SKG => SK1
        character(*,SKG)    , intent(in)    :: source(:,:)
        character(*,SKG)    , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_IK5(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_IK5
#endif
        use pm_kind, only: LK, IKG => IK5
        integer(IKG)        , intent(in)    :: source(:,:)
        integer(IKG)        , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if IK4_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_IK4(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_IK4
#endif
        use pm_kind, only: LK, IKG => IK4
        integer(IKG)        , intent(in)    :: source(:,:)
        integer(IKG)        , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if IK3_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_IK3(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_IK3
#endif
        use pm_kind, only: LK, IKG => IK3
        integer(IKG)        , intent(in)    :: source(:,:)
        integer(IKG)        , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if IK2_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_IK2(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_IK2
#endif
        use pm_kind, only: LK, IKG => IK2
        integer(IKG)        , intent(in)    :: source(:,:)
        integer(IKG)        , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if IK1_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_IK1(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_IK1
#endif
        use pm_kind, only: LK, IKG => IK1
        integer(IKG)        , intent(in)    :: source(:,:)
        integer(IKG)        , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_LK5(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_LK5
#endif
        use pm_kind, only: LK, LKG => LK5
        logical(LKG)        , intent(in)    :: source(:,:)
        logical(LKG)        , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if LK4_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_LK4(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_LK4
#endif
        use pm_kind, only: LK, LKG => LK4
        logical(LKG)        , intent(in)    :: source(:,:)
        logical(LKG)        , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if LK3_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_LK3(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_LK3
#endif
        use pm_kind, only: LK, LKG => LK3
        logical(LKG)        , intent(in)    :: source(:,:)
        logical(LKG)        , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if LK2_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_LK2(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_LK2
#endif
        use pm_kind, only: LK, LKG => LK2
        logical(LKG)        , intent(in)    :: source(:,:)
        logical(LKG)        , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if LK1_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_LK1(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_LK1
#endif
        use pm_kind, only: LK, LKG => LK1
        logical(LKG)        , intent(in)    :: source(:,:)
        logical(LKG)        , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_CK5(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_CK5
#endif
        use pm_kind, only: LK, CKG => CK5
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_CK4(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_CK4
#endif
        use pm_kind, only: LK, CKG => CK4
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_CK3(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_CK3
#endif
        use pm_kind, only: LK, CKG => CK3
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_CK2(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_CK2
#endif
        use pm_kind, only: LK, CKG => CK2
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_CK1(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_CK1
#endif
        use pm_kind, only: LK, CKG => CK1
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_RK5(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_RK5
#endif
        use pm_kind, only: LK, RKG => RK5
        real(RKG)           , intent(in)    :: source(:,:)
        real(RKG)           , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if RK4_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_RK4(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_RK4
#endif
        use pm_kind, only: LK, RKG => RK4
        real(RKG)           , intent(in)    :: source(:,:)
        real(RKG)           , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if RK3_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_RK3(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_RK3
#endif
        use pm_kind, only: LK, RKG => RK3
        real(RKG)           , intent(in)    :: source(:,:)
        real(RKG)           , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if RK2_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_RK2(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_RK2
#endif
        use pm_kind, only: LK, RKG => RK2
        real(RKG)           , intent(in)    :: source(:,:)
        real(RKG)           , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if RK1_ENABLED
    PURE recursive module subroutine setMatTransSymmNewArb_RK1(source, destin, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransSymmNewArb_RK1
#endif
        use pm_kind, only: LK, RKG => RK1
        real(RKG)           , intent(in)    :: source(:,:)
        real(RKG)           , intent(out)   :: destin(:,:)
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface
    !endcond excluded

    !   Arbitrary block size Hermitian transpose.
    !cond excluded
    interface setMatTrans

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setMatTransHermOldArb_CK5(source, operation, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermOldArb_CK5
#endif
        use pm_kind, only: LK, CKG => CK5
        complex(CKG)        , intent(inout) :: source(:,:)
        type(transHerm_type), intent(in)    :: operation
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setMatTransHermOldArb_CK4(source, operation, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermOldArb_CK4
#endif
        use pm_kind, only: LK, CKG => CK4
        complex(CKG)        , intent(inout) :: source(:,:)
        type(transHerm_type), intent(in)    :: operation
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setMatTransHermOldArb_CK3(source, operation, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermOldArb_CK3
#endif
        use pm_kind, only: LK, CKG => CK3
        complex(CKG)        , intent(inout) :: source(:,:)
        type(transHerm_type), intent(in)    :: operation
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setMatTransHermOldArb_CK2(source, operation, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermOldArb_CK2
#endif
        use pm_kind, only: LK, CKG => CK2
        complex(CKG)        , intent(inout) :: source(:,:)
        type(transHerm_type), intent(in)    :: operation
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setMatTransHermOldArb_CK1(source, operation, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermOldArb_CK1
#endif
        use pm_kind, only: LK, CKG => CK1
        complex(CKG)        , intent(inout) :: source(:,:)
        type(transHerm_type), intent(in)    :: operation
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE recursive module subroutine setMatTransHermNewArb_CK5(source, destin, operation, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermNewArb_CK5
#endif
        use pm_kind, only: LK, CKG => CK5
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
        type(transHerm_type), intent(in)    :: operation
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if CK4_ENABLED
    PURE recursive module subroutine setMatTransHermNewArb_CK4(source, destin, operation, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermNewArb_CK4
#endif
        use pm_kind, only: LK, CKG => CK4
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
        type(transHerm_type), intent(in)    :: operation
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if CK3_ENABLED
    PURE recursive module subroutine setMatTransHermNewArb_CK3(source, destin, operation, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermNewArb_CK3
#endif
        use pm_kind, only: LK, CKG => CK3
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
        type(transHerm_type), intent(in)    :: operation
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if CK2_ENABLED
    PURE recursive module subroutine setMatTransHermNewArb_CK2(source, destin, operation, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermNewArb_CK2
#endif
        use pm_kind, only: LK, CKG => CK2
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
        type(transHerm_type), intent(in)    :: operation
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

#if CK1_ENABLED
    PURE recursive module subroutine setMatTransHermNewArb_CK1(source, destin, operation, bsize)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatTransHermNewArb_CK1
#endif
        use pm_kind, only: LK, CKG => CK1
        complex(CKG)        , intent(in)    :: source(:,:)
        complex(CKG)        , intent(out)   :: destin(:,:)
        type(transHerm_type), intent(in)    :: operation
        integer(IK)         , intent(in)    :: bsize
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface
    !endcond excluded

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_matrixTrans ! LCOV_EXCL_LINE