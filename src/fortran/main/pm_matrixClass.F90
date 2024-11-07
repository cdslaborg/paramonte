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
!>  This module contains abstract and concrete derived types that are required for compile-time
!>  resolution of procedures within the generic interfaces of the ParaMonte library for Linear Algebra operations.<br>
!>  Such procedures frequently need to work on different classes of matrices (e.g., Symmetric, Hermitian, Unitary, Orthogonal, ...) as their input matrix arguments.<br>
!>
!>  \details
!>  A large important group of matrices are defined based on their transposition properties.
!>  In linear algebra, the **transpose of a matrix** is an operator which flips a matrix over its diagonal.<br>
!>  It switches the row and column indices of a given matrix \f$A\f$ by producing another (transposed) matrix.<br>
!>  The transpose of a matrix was introduced in 1858 by the British mathematician Arthur Cayley.<br>
!>
!>  **Notation**<br>
!>  The transpose of a given matrix \f$A\f$ is frequently denoted by \f$A^T\f$.<br>
!>  In the case of square matrices, \f$A^T\f$ may also denote the \f$T\f$th power of the matrix \f$A\f$.<br>
!>  To avoid a possible confusion, some use left upperscripts to denote transpose, that is, \f${}^TA\f$.<br>
!>  An advantage of this notation is that no parentheses are needed when exponents are involved as \f$({}^TA)^n = {}^T(A^n)\f$, notation \f${}^TA^n\f$ is not ambiguous.<br>
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
!>  **Matrix definitions involving transposition**<br>
!>  <ol>
!>      <li>    **Symmetric Matrix**<br>
!>              A square matrix whose transpose is equal to itself is called a symmetric matrix.<br>
!>              In other words, \f$A\f$ is symmetric if \f$\mathbf{A}^{\up{T}} = \mathbf{A}\f$.<br>
!>      <li>    **Skew-Symmetric Matrix**<br>
!>              A square matrix whose transpose is equal to its negative is called a skew-symmetric matrix
!>              In other words, \f$A\f$ is skew-symmetric if \f$\mathbf{A}^{\up{T}} = -\mathbf{A}\f$.<br>
!>      <li>    **Hermitian Matrix**<br>
!>              A square **complex matrix** whose transpose equals the matrix with every entry replaced by its complex conjugate (\f$\overline{\cdot}\f$) is called a Hermitian matrix.<br>
!>              In other words, \f$A\f$ is Hermitian if \f$\mathbf{A}^{\up{T}} = \overline{\mathbf{A}}\f$.<br>
!>              This special kind of transpose is also called a **Hermitian transpose** or **Conjugate transpose**.<br>
!>      <li>    **Skew-Hermitian Matrix**<br>
!>              A square **complex matrix** whose transpose is equal to the negation of its complex conjugate is called a skew-Hermitian matrix.<br>
!>              In other words, \f$A\f$ is skew-Hermitian if \f$\mathbf{A}^{\up{T}} = -{\overline{\mathbf{A}}}\f$.<br>
!>      <li>    **Unitary Matrix**<br>
!>              A square complex matrix whose transpose is equal to its conjugate inverse is called a unitary matrix.<br>
!>              In other words, \f$A\f$ is unitary if \f$\mathbf{A}^{\up{T}} = {\overline{\mathbf{A}^{-1}}}\f$.<br>
!>      <li>    **Orthogonal Matrix**<br>
!>              A square matrix whose transpose is equal to its inverse is called an orthogonal matrix.<br>
!>              In other words, \f$A\f$ is orthogonal if \f$\mathbf{A}^{\up{T}} = \mathbf{A}^{-1}\f$.<br>
!>  </ol>
!>
!>  \todo
!>  \phigh
!>  This module is a work in progress.<br>
!>  Various missing classifications of matrices should be added in the future and the existing classifications be refined.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_matrixClass

    use pm_kind, only: SK, IK, LK
    use pm_matrixSubset, only: uppDia, uppDia_type
    use pm_matrixSubset, only: lowDia, lowDia_type
    use pm_matrixPack, only: rdpack, rdpack_type
    use pm_matrixPack, only: rfpack, rfpack_type

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_matrixClass"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an `abstract` derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require different matrix classes (e.g., Symmetric, Hermitian, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>
    !>  \final{matrix_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, abstract :: matrix_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the general rectangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [genrecmat](@ref pm_matrixClass::genrecmat)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{genrecmat_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type :: genrecmat_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [genrecmat_type](@ref pm_matrixClass::genrecmat_type) that is exclusively used
    !>  to signify the general rectangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{genrecmat}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(genrecmat_type), parameter :: genrecmat = genrecmat_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: genrecmat
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the Square class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [square](@ref pm_matrixClass::square)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{square_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(matrix_type) :: square_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [square_type](@ref pm_matrixClass::square_type) that is exclusively used
    !>  to signify the Square class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{square}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(square_type), parameter :: square = square_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: square
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the Invertible class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [invertible](@ref pm_matrixClass::invertible)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{invertible_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(matrix_type) :: invertible_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [invertible_type](@ref pm_matrixClass::invertible_type) that is exclusively used
    !>  to signify the Invertible class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{invertible}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(invertible_type), parameter :: invertible = invertible_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: invertible
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the Factorization class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [factoring](@ref pm_matrixClass::factoring)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{factoring_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(matrix_type) :: factoring_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [factoring_type](@ref pm_matrixClass::factoring_type) that is exclusively used
    !>  to signify the Factorization class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{factoring}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(factoring_type), parameter :: factoring = factoring_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: factoring
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the Cholesky Factorization class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [cholesky](@ref pm_matrixClass::cholesky)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{cholesky_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(matrix_type) :: cholesky_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [cholesky_type](@ref pm_matrixClass::cholesky_type) that is exclusively used
    !>  to signify the Cholesky Factorization class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{cholesky}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(cholesky_type), parameter :: cholesky = cholesky_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: cholesky
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the lower-triangle Cholesky Factorization class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [choLow](@ref pm_matrixClass::choLow)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{choLow_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(cholesky_type) :: choLow_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [cholesky_type](@ref pm_matrixClass::cholesky_type) that is exclusively used
    !>  to signify lower-triangle Cholesky Factorization class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{choLow}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(choLow_type), parameter :: choLow = choLow_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: choLow
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the upper-triangle Cholesky Factorization class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [choUpp](@ref pm_matrixClass::choUpp)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{choUpp_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(cholesky_type) :: choUpp_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [cholesky_type](@ref pm_matrixClass::cholesky_type) that is exclusively used
    !>  to signify upper-triangle Cholesky Factorization class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{choUpp}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(choUpp_type), parameter :: choUpp = choUpp_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: choUpp
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the LUP Factorization class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [lup](@ref pm_matrixClass::lup)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{lup_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(matrix_type) :: lup_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [lup_type](@ref pm_matrixClass::lup_type) that is exclusively used
    !>  to signify the LUP Factorization class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{lup}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(lup_type), parameter :: lup = lup_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: lup
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the Symmetric class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [symmetric](@ref pm_matrixClass::symmetric)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{symmetric_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(matrix_type) :: symmetric_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [symmetric_type](@ref pm_matrixClass::symmetric_type) that is exclusively used
    !>  to signify the Symmetric class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{symmetric}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(symmetric_type), parameter :: symmetric = symmetric_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: symmetric
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the Hermitian class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [hermitian](@ref pm_matrixClass::hermitian)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{hermitian_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(matrix_type) :: hermitian_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [hermitian_type](@ref pm_matrixClass::hermitian_type) that is exclusively used
    !>  to signify the Hermitian class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{hermitian}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(hermitian_type), parameter :: hermitian = hermitian_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: hermitian
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify the Hermitian Positive-Definite class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [posdefmat](@ref pm_matrixClass::posdefmat)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{posdefmat_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(hermitian_type) :: posdefmat_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [hermitian_type](@ref pm_matrixClass::hermitian_type) that is exclusively used
    !>  to signify the Hermitian class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{posdefmat}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(posdefmat_type), parameter :: posdefmat = posdefmat_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: posdefmat
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the
    !>  triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  A matrix of the form,
    !>  \f{equation}{
    !>      L =
    !>      \begin{bmatrix}
    !>          \ell_{1,1}&&&&0 \\
    !>          \ell_{2,1}&\ell_{2,2}&&&\\
    !>          \ell_{3,1}&\ell_{3,2}&\ddots &&\\
    !>          \vdots &\vdots &\ddots &\ddots &\\
    !>          \ell_{n,1}&\ell_{n,2}&\ldots &\ell_{n,n-1}&\ell_{n,n}
    !>      \end{bmatrix}
    !>  \f}
    !>  is called a **lower triangular matrix** or **left triangular matrix**, and analogously a matrix of the form,
    !>  \f{equation}{
    !>      U =
    !>      \begin{bmatrix}
    !>          u_{1,1}&u_{1,2}&u_{1,3}&\ldots &u_{1,n}\\
    !>          &u_{2,2}&u_{2,3}&\ldots &u_{2,n}\\
    !>          &&\ddots &\ddots &\vdots\\
    !>          &&&\ddots &u_{n-1,n}\\
    !>          0&&&&u_{n,n}
    !>      \end{bmatrix}
    !>  \f}
    !>  is called an **upper triangular matrix** or **right triangular matrix**.<br>
    !>  A lower or left triangular matrix is commonly denoted with the variable **L**, and an upper or right triangular matrix is commonly denoted with the variable **U**.
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [triang](@ref pm_matrixClass::triang)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{triang_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(matrix_type) :: triang_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [triang_type](@ref pm_matrixClass::triang_type) that is exclusively used to signify the
    !>  triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{triang}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(triang_type), parameter :: triang = triang_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: triang
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the
    !>  upper-diagonal triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the documentation of [triang_type](@ref pm_matrixClass::triang_type) for more details on triangular matrices.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [upperDiag](@ref pm_matrixClass::upperDiag)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{upperDiag_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(triang_type) :: upperDiag_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [upperDiag_type](@ref pm_matrixClass::upperDiag_type) that is exclusively used to signify the
    !>  upper-diagonal triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{upperDiag}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(upperDiag_type), parameter :: upperDiag = upperDiag_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: upperDiag
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the
    !>  upper-unit-diagonal triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the documentation of [triang_type](@ref pm_matrixClass::triang_type) for more details on triangular matrices.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [upperUnit](@ref pm_matrixClass::upperUnit)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{upperUnit_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(upperDiag_type) :: upperUnit_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [upperUnit_type](@ref pm_matrixClass::upperUnit_type) that is exclusively used to signify the
    !>  upper-unit-diagonal triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{upperUnit}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(upperUnit_type), parameter :: upperUnit = upperUnit_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: upperUnit
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the
    !>  upper-zero-diagonal triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the documentation of [triang_type](@ref pm_matrixClass::triang_type) for more details on triangular matrices.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [upperZero](@ref pm_matrixClass::upperZero)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{upperZero_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(upperDiag_type) :: upperZero_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [upperUnit_type](@ref pm_matrixClass::upperZero_type) that is exclusively used to signify the
    !>  upper-zero-diagonal triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{upperZero}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(upperZero_type), parameter :: upperZero = upperZero_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: upperZero
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the
    !>  lower-diagonal triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the documentation of [triang_type](@ref pm_matrixClass::triang_type) for more details on triangular matrices.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [lowerDiag](@ref pm_matrixClass::lowerDiag)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{lowerDiag_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(triang_type) :: lowerDiag_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type) that is exclusively used to signify the
    !>  lower-diagonal triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{lowerDiag}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(lowerDiag_type), parameter :: lowerDiag = lowerDiag_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: lowerDiag
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the
    !>  lower-unit-diagonal triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the documentation of [triang_type](@ref pm_matrixClass::triang_type) for more details on triangular matrices.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [lowerUnit](@ref pm_matrixClass::lowerUnit)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{lowerUnit_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(lowerDiag_type) :: lowerUnit_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [upperUnit_type](@ref pm_matrixClass::upperUnit_type) that is exclusively used to signify the
    !>  lower-unit-diagonal triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{lowerUnit}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(lowerUnit_type), parameter :: lowerUnit = lowerUnit_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: lowerUnit
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the
    !>  lower-zero-diagonal triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  See the documentation of [triang_type](@ref pm_matrixClass::triang_type) for more details on triangular matrices.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [lowerZero](@ref pm_matrixClass::lowerZero)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{lowerZero_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(lowerDiag_type) :: lowerZero_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [lowerUnit_type](@ref pm_matrixClass::lowerZero_type) that is exclusively used to signify the
    !>  lower-zero-diagonal triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{lowerZero}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(lowerZero_type), parameter :: lowerZero = lowerZero_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: lowerZero
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the
    !>  unit-triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  An **atomic (upper or lower) triangular matrix** is a special form of [unitriangular matrix](@ref pm_matrixClass::unitTriang_type),
    !>  where all of the off-diagonal elements are zero, except for the entries in a single column or row.<br>
    !>  Such a matrix is also called a **Frobenius matrix** or **Gauss matrix**, or **Gauss transformation matrix**.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [unitTriang](@ref pm_matrixClass::unitTriang)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{unitTriang_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(triang_type) :: unitTriang_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [unitTriang_type](@ref pm_matrixClass::unitTriang_type) that is exclusively used to signify the
    !>  unit-triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{unitTriang}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(unitTriang_type), parameter :: unitTriang = unitTriang_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: unitTriang
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the
    !>  atomic-triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  An **atomic (upper or lower) triangular matrix** is a special form of [unitriangular matrix](@ref pm_matrixClass::unitTriang_type),
    !>  where all of the off-diagonal elements are zero, except for the entries in a single column or row.<br>
    !>  Such a matrix is also called a **Frobenius matrix** or **Gauss matrix**, or **Gauss transformation matrix**.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [atomicTriang](@ref pm_matrixClass::atomicTriang)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{atomicTriang_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(unitTriang_type) :: atomicTriang_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type) that is exclusively used to signify the
    !>  atomic-triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{atomicTriang}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(atomicTriang_type), parameter :: atomicTriang = atomicTriang_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: atomicTriang
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the
    !>  frobenius class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  A Frobenius matrix is a special kind of square matrix from numerical mathematics.<br>
    !>  A matrix is a Frobenius matrix if it has the following three properties:<br>
    !>  <ol>
    !>      <li>    All entries on the main diagonal are ones.<br>
    !>      <li>    The entries below the main diagonal of at most one column are arbitrary.<br>
    !>      <li>    Every other entry is zero.<br>
    !>  </ol>
    !>  The following matrix is an example.<br>
    !>  \f{equation}{
    !>      A =
    !>      \begin{pmatrix}
    !>          1&0&0&\cdots &0\\
    !>          0&1&0&\cdots &0\\
    !>          0&a_{32}&1&\cdots &0\\
    !>          \vdots &\vdots &\vdots &\ddots &\vdots \\
    !>          0&a_{n2}&0&\cdots &1
    !>      \end{pmatrix}
    !>  \f}
    !>  Frobenius matrices are invertible.<br>
    !>  The inverse of a Frobenius matrix is again a Frobenius matrix, equal to the original matrix with changed signs outside the main diagonal.<br>
    !>  The inverse of the example above is therefore:
    !>  \f{equation}{
    !>      A^{-1} =
    !>      \begin{pmatrix}
    !>          1&0&0&\cdots &0\\
    !>          0&1&0&\cdots &0\\
    !>          0&-a_{{32}}&1&\cdots &0\\
    !>          \vdots &\vdots &\vdots &\ddots &\vdots \\
    !>          0&-a_{{n2}}&0&\cdots &1
    !>      \end{pmatrix}
    !>  \f}
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [frobenius](@ref pm_matrixClass::frobenius)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{frobenius_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(atomicTriang_type) :: frobenius_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [frobenius_type](@ref pm_matrixClass::frobenius_type) that is exclusively used to signify the
    !>  atomic-triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{frobenius}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(frobenius_type), parameter :: frobenius = frobenius_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: frobenius
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used to signify the
    !>  Gauss-Transformation class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  A **Gauss matrix** is a special form of the [atomic triangular matrix](@ref pm_matrixClass::atomicTriang_type)
    !>  that differs from an Identity matrix only in the elements of a single row preceding the diagonal entry of that row
    !>  (as opposed to the [Frobenius matrix](@ref pm_matrixClass::frobenius_type) definition which has the matrix differing
    !>  from the identity matrix in a single column below the diagonal).<br>
    !>  The following example shows a 4-by-4 Gauss matrix with its 3rd row differing from the identity matrix.<br>
    !>  \f{equation}{
    !>      A =
    !>      \begin{pmatrix}
    !>          1&0&0&0\\
    !>          0&1&0&0\\
    !>          a_{31}&a_{32}&1&0\\
    !>          0&0&0&1
    !>      \end{pmatrix}
    !>  \f}
    !>  An alternative name for the Gauss matrix is **Gauss transformation matrix**, after Carl Friedrich Gauss.<br>
    !>  The Gauss matrix is used in the process of Gaussian elimination to represent the Gaussian transformations.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [gauss](@ref pm_matrixClass::gauss)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{gauss_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type, extends(atomicTriang_type) :: gauss_type
    end type

    !>  \brief
    !>  This is a scalar `parameter` object of type [gauss_type](@ref pm_matrixClass::gauss_type) that is exclusively used to signify the
    !>  atomic-triangular class of a given matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [genrecmat](@ref pm_matrixClass::genrecmat)<br>
    !>  [posdefmat](@ref pm_matrixClass::posdefmat)<br>
    !>  [symmetric](@ref pm_matrixClass::symmetric)<br>
    !>  [hermitian](@ref pm_matrixClass::hermitian)<br>
    !>  [triang_type](@ref pm_matrixClass::triang_type)<br>
    !>  [posdefmat_type](@ref pm_matrixClass::posdefmat_type)<br>
    !>  [symmetric_type](@ref pm_matrixClass::symmetric_type)<br>
    !>  [hermitian_type](@ref pm_matrixClass::hermitian_type)<br>
    !>  [upperDiag_type](@ref pm_matrixClass::upperDiag_type)<br>
    !>  [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type)<br>
    !>  [upperZero_type](@ref pm_matrixClass::upperZero_type)<br>
    !>  [lowerZero_type](@ref pm_matrixClass::lowerZero_type)<br>
    !>  [upperUnit_type](@ref pm_matrixClass::upperUnit_type)<br>
    !>  [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type)<br>
    !>  [unitTriang_type](@ref pm_matrixClass::unitTriang_type)<br>
    !>  [atomicTriang_type](@ref pm_matrixClass::atomicTriang_type)<br>
    !>  [frobenius_type](@ref pm_matrixClass::frobenius_type)<br>
    !>  [genrecmat_type](@ref pm_matrixClass::genrecmat_type)<br>
    !>  [matrix_type](@ref pm_matrixClass::matrix_type)<br>
    !>  [gauss_type](@ref pm_matrixClass::gauss_type)<br>
    !>
    !>  \final{gauss}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    type(gauss_type), parameter :: gauss = gauss_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: gauss
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if and only if the input matrix is of the specified input `class`.<br>
    !>
    !>  \details
    !>  See [pm_matrixClass](@ref pm_matrixClass) for the mathematical definitions of different matrix classes.<br>
    !>
    !>  \param[in]  mat     :   The input matrix of arbitrary shape `(:,:)` of,
    !>                          <ol>
    !>                              <li>    type `character` of kind \SKALL of arbitrary length type parameter, or <br>
    !>                              <li>    type `integer` of kind \IKALL, or <br>
    !>                              <li>    type `logical` of kind \LKALL, or <br>
    !>                              <li>    type `complex` of kind \CKALL, or <br>
    !>                              <li>    type `real` of kind \RKALL, <br>
    !>                          </ol>
    !>                          containing the matrix whose class membership is to be tested.<br>
    !>  \param[in]  class   :   The input scalar constant that can be one of the following:<br>
    !>                          <ol>
    !>                              <li>    The scalar constant [symmetric](@ref pm_matrixClass::symmetric) implying the matrix is square **Symmetric**.<br>
    !>                              <li>    The scalar constant [hermitian](@ref pm_matrixClass::hermitian) implying the matrix is square **Hermitian**.<br>
    !>                              <li>    The scalar constant [posdefmat](@ref pm_matrixClass::posdefmat) implying the matrix is square Hermitian **Positive Definite**.<br>
    !>                          </ol>
    !>  \param[in]  subset  :   The input scalar constant that can be one of the following:<br>
    !>                          <ol>
    !>                              <li>    The scalar constant [uppDia](@ref pm_matrixSubset::uppDia) implying the upper-diagonal subset of the input matrix must be used for testing.<br>
    !>                              <li>    The scalar constant [lowDia](@ref pm_matrixSubset::lowDia) implying the lower-diagonal subset of the input matrix must be used for testing.<br>
    !>                          </ol>
    !>                          Specifying this argument leads to faster runtimes.<br>
    !>                          If missing, the full input matrix will be considered for positive-definiteness testing.<br>
    !>                          (**optional**. It can be present **only if** the input argument `class` is set to [posdefmat](@ref pm_matrixClass::posdefmat).)
    !>  \param[in]  pack    :   The input scalar constant that can be one of the following:<br>
    !>                          <ol>
    !>                              <li>    The scalar constant [rdpack](@ref pm_matrixPack::rdpack) implying the Rectangular Default Packing format of the input matrix.<br>
    !>                              <li>    The scalar constant [rfpack](@ref pm_matrixPack::rfpack) implying the Rectangular Full Packing format of the input matrix.<br>
    !>                          </ol>
    !>                          (**optional**. It must be present **if and only if** the input arguments `class = posdefmat` and the input argument `subset` is present.)
    !>
    !>  \return
    !>  `itis`              :   The output scalar `logical` of default kind \LK that is `.true.` if and only if the input matrix belongs to the specified matrix type `class`.<br>
    !>
    !>  \interface{isMatClass}
    !>  \code{.F90}
    !>
    !>      use pm_matrixClass, only: isMatClass
    !>
    !>      itis = isMatClass(mat(:,:), class)
    !>      itis = isMatClass(mat(:,:), class, subset, pack) ! class = posdefmat; subset = uppDia, lowDia; pack = rdpack, rfpack
    !>      !
    !>  \endcode
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [pm_matrixChol](@ref pm_matrixChol)<br>
    !>  [pm_matrixPack](@ref pm_matrixPack)<br>
    !>  [pm_matrixSubset](@ref pm_matrixSubset)<br>
    !>
    !>  \example
    !>  \include{lineno} example/pm_matrixClass/isMatClass/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_matrixClass/isMatClass/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixClass](@ref test_pm_matrixClass)
    !>
    !>  \todo
    !>  \phigh
    !>  This generic interface must be extended to all matrix classes documented in this module.<br>
    !>
    !>  \todo
    !>  \phigh
    !>  The implementation for [rfpack](@ref pm_matrixPack::rfpack) can be improved once the corresponding improvements to the auxiliary routines used are implemented.<br>
    !>  For example, the current implementation for positive-definiteness check makes a copy of the input array which can be avoided if
    !>  the corresponding [setMatChol](@ref pm_matrixChol::setMatChol) interface is implemented.<br>
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    interface isMatClass

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function isMatClassPosDefFulRDP_CK5(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefFulRDP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if CK4_ENABLED
    PURE module function isMatClassPosDefFulRDP_CK4(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefFulRDP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if CK3_ENABLED
    PURE module function isMatClassPosDefFulRDP_CK3(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefFulRDP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if CK2_ENABLED
    PURE module function isMatClassPosDefFulRDP_CK2(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefFulRDP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if CK1_ENABLED
    PURE module function isMatClassPosDefFulRDP_CK1(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefFulRDP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function isMatClassPosDefFulRDP_RK5(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefFulRDP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if RK4_ENABLED
    PURE module function isMatClassPosDefFulRDP_RK4(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefFulRDP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if RK3_ENABLED
    PURE module function isMatClassPosDefFulRDP_RK3(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefFulRDP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if RK2_ENABLED
    PURE module function isMatClassPosDefFulRDP_RK2(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefFulRDP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if RK1_ENABLED
    PURE module function isMatClassPosDefFulRDP_RK1(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefFulRDP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function isMatClassPosDefUppRDP_CK5(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRDP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if CK4_ENABLED
    PURE module function isMatClassPosDefUppRDP_CK4(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRDP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if CK3_ENABLED
    PURE module function isMatClassPosDefUppRDP_CK3(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRDP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if CK2_ENABLED
    PURE module function isMatClassPosDefUppRDP_CK2(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRDP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if CK1_ENABLED
    PURE module function isMatClassPosDefUppRDP_CK1(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRDP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function isMatClassPosDefUppRDP_RK5(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRDP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if RK4_ENABLED
    PURE module function isMatClassPosDefUppRDP_RK4(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRDP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if RK3_ENABLED
    PURE module function isMatClassPosDefUppRDP_RK3(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRDP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if RK2_ENABLED
    PURE module function isMatClassPosDefUppRDP_RK2(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRDP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if RK1_ENABLED
    PURE module function isMatClassPosDefUppRDP_RK1(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRDP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function isMatClassPosDefLowRDP_CK5(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRDP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if CK4_ENABLED
    PURE module function isMatClassPosDefLowRDP_CK4(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRDP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if CK3_ENABLED
    PURE module function isMatClassPosDefLowRDP_CK3(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRDP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if CK2_ENABLED
    PURE module function isMatClassPosDefLowRDP_CK2(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRDP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if CK1_ENABLED
    PURE module function isMatClassPosDefLowRDP_CK1(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRDP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function isMatClassPosDefLowRDP_RK5(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRDP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if RK4_ENABLED
    PURE module function isMatClassPosDefLowRDP_RK4(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRDP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if RK3_ENABLED
    PURE module function isMatClassPosDefLowRDP_RK3(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRDP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if RK2_ENABLED
    PURE module function isMatClassPosDefLowRDP_RK2(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRDP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if RK1_ENABLED
    PURE module function isMatClassPosDefLowRDP_RK1(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRDP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rdpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function isMatClassPosDefUppRFP_CK5(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRFP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if CK4_ENABLED
    PURE module function isMatClassPosDefUppRFP_CK4(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRFP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if CK3_ENABLED
    PURE module function isMatClassPosDefUppRFP_CK3(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRFP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if CK2_ENABLED
    PURE module function isMatClassPosDefUppRFP_CK2(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRFP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if CK1_ENABLED
    PURE module function isMatClassPosDefUppRFP_CK1(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRFP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function isMatClassPosDefUppRFP_RK5(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRFP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if RK4_ENABLED
    PURE module function isMatClassPosDefUppRFP_RK4(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRFP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if RK3_ENABLED
    PURE module function isMatClassPosDefUppRFP_RK3(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRFP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if RK2_ENABLED
    PURE module function isMatClassPosDefUppRFP_RK2(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRFP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if RK1_ENABLED
    PURE module function isMatClassPosDefUppRFP_RK1(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefUppRFP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(uppDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function isMatClassPosDefLowRFP_CK5(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRFP_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if CK4_ENABLED
    PURE module function isMatClassPosDefLowRFP_CK4(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRFP_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if CK3_ENABLED
    PURE module function isMatClassPosDefLowRFP_CK3(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRFP_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if CK2_ENABLED
    PURE module function isMatClassPosDefLowRFP_CK2(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRFP_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if CK1_ENABLED
    PURE module function isMatClassPosDefLowRFP_CK1(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRFP_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function isMatClassPosDefLowRFP_RK5(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRFP_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if RK4_ENABLED
    PURE module function isMatClassPosDefLowRFP_RK4(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRFP_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if RK3_ENABLED
    PURE module function isMatClassPosDefLowRFP_RK3(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRFP_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if RK2_ENABLED
    PURE module function isMatClassPosDefLowRFP_RK2(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRFP_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

#if RK1_ENABLED
    PURE module function isMatClassPosDefLowRFP_RK1(mat, class, subset, pack) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassPosDefLowRFP_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(posdefmat_type)    , intent(in)                    :: class
        type(lowDia_type)       , intent(in)                    :: subset
        type(rfpack_type)       , intent(in)                    :: pack
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function isMatClassSymm_SK5(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if SK4_ENABLED
    PURE module function isMatClassSymm_SK4(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if SK3_ENABLED
    PURE module function isMatClassSymm_SK3(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if SK2_ENABLED
    PURE module function isMatClassSymm_SK2(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if SK1_ENABLED
    PURE module function isMatClassSymm_SK1(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function isMatClassSymm_IK5(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if IK4_ENABLED
    PURE module function isMatClassSymm_IK4(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if IK3_ENABLED
    PURE module function isMatClassSymm_IK3(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if IK2_ENABLED
    PURE module function isMatClassSymm_IK2(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if IK1_ENABLED
    PURE module function isMatClassSymm_IK1(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function isMatClassSymm_LK5(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if LK4_ENABLED
    PURE module function isMatClassSymm_LK4(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if LK3_ENABLED
    PURE module function isMatClassSymm_LK3(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if LK2_ENABLED
    PURE module function isMatClassSymm_LK2(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if LK1_ENABLED
    PURE module function isMatClassSymm_LK1(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function isMatClassSymm_CK5(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if CK4_ENABLED
    PURE module function isMatClassSymm_CK4(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if CK3_ENABLED
    PURE module function isMatClassSymm_CK3(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if CK2_ENABLED
    PURE module function isMatClassSymm_CK2(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if CK1_ENABLED
    PURE module function isMatClassSymm_CK1(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function isMatClassSymm_RK5(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if RK4_ENABLED
    PURE module function isMatClassSymm_RK4(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if RK3_ENABLED
    PURE module function isMatClassSymm_RK3(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if RK2_ENABLED
    PURE module function isMatClassSymm_RK2(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if RK1_ENABLED
    PURE module function isMatClassSymm_RK1(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassSymm_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(symmetric_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if SK5_ENABLED
    PURE module function isMatClassHerm_SK5(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_SK5
#endif
        use pm_kind, only: SKG => SK5
        character(*,SKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if SK4_ENABLED
    PURE module function isMatClassHerm_SK4(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_SK4
#endif
        use pm_kind, only: SKG => SK4
        character(*,SKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if SK3_ENABLED
    PURE module function isMatClassHerm_SK3(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_SK3
#endif
        use pm_kind, only: SKG => SK3
        character(*,SKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if SK2_ENABLED
    PURE module function isMatClassHerm_SK2(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_SK2
#endif
        use pm_kind, only: SKG => SK2
        character(*,SKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if SK1_ENABLED
    PURE module function isMatClassHerm_SK1(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_SK1
#endif
        use pm_kind, only: SKG => SK1
        character(*,SKG)        , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if IK5_ENABLED
    PURE module function isMatClassHerm_IK5(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_IK5
#endif
        use pm_kind, only: IKG => IK5
        integer(IKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if IK4_ENABLED
    PURE module function isMatClassHerm_IK4(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_IK4
#endif
        use pm_kind, only: IKG => IK4
        integer(IKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if IK3_ENABLED
    PURE module function isMatClassHerm_IK3(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_IK3
#endif
        use pm_kind, only: IKG => IK3
        integer(IKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if IK2_ENABLED
    PURE module function isMatClassHerm_IK2(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_IK2
#endif
        use pm_kind, only: IKG => IK2
        integer(IKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if IK1_ENABLED
    PURE module function isMatClassHerm_IK1(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_IK1
#endif
        use pm_kind, only: IKG => IK1
        integer(IKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if LK5_ENABLED
    PURE module function isMatClassHerm_LK5(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_LK5
#endif
        use pm_kind, only: LKG => LK5
        logical(LKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if LK4_ENABLED
    PURE module function isMatClassHerm_LK4(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_LK4
#endif
        use pm_kind, only: LKG => LK4
        logical(LKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if LK3_ENABLED
    PURE module function isMatClassHerm_LK3(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_LK3
#endif
        use pm_kind, only: LKG => LK3
        logical(LKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if LK2_ENABLED
    PURE module function isMatClassHerm_LK2(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_LK2
#endif
        use pm_kind, only: LKG => LK2
        logical(LKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if LK1_ENABLED
    PURE module function isMatClassHerm_LK1(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_LK1
#endif
        use pm_kind, only: LKG => LK1
        logical(LKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function isMatClassHerm_CK5(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if CK4_ENABLED
    PURE module function isMatClassHerm_CK4(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if CK3_ENABLED
    PURE module function isMatClassHerm_CK3(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if CK2_ENABLED
    PURE module function isMatClassHerm_CK2(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if CK1_ENABLED
    PURE module function isMatClassHerm_CK1(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)            , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function isMatClassHerm_RK5(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if RK4_ENABLED
    PURE module function isMatClassHerm_RK4(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if RK3_ENABLED
    PURE module function isMatClassHerm_RK3(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if RK2_ENABLED
    PURE module function isMatClassHerm_RK2(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

#if RK1_ENABLED
    PURE module function isMatClassHerm_RK1(mat, class) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatClassHerm_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)               , intent(in)    , contiguous    :: mat(:,:)
        type(hermitian_type)    , intent(in)                    :: class
        logical(LK)                                             :: itis
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_matrixClass ! LCOV_EXCL_LINE