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
!>  Such procedures frequently need to work with different packaging of Symmetric/Hermitian/Band 
!>  or other special types of their input matrix arguments.<br>
!>
!>  \details
!>  This module follows the conventions of the LAPACK Linear Algebra Fortran software for matrix package schemes.<br>
!>  LAPACK routines use the following matrix package schemes:<br>
!>  <ol>
!>      <li>    <b>Linear Full Packing (**LFP**) for Symmetric, Hermitian triangular matrices ([lfpack](@ref pm_matrixPack::lfpack))</b><br>
!>              This vector packing format compactly stores matrix elements as a vector.<br>
!>              This is useful when only one part of the matrix, the upper or lower triangle, is necessary to determine all matrix elements.<br>
!>              This is the case when the matrix is Symmetric, Hermitian, or is of [upper/lower triangular class](@ref pm_matrixClass).<br>
!>              A **lower-diagonal** triangular or Symmetric/Hermitian matrix of the form,<br>
!>              \f{equation}{
!>                  L = \begin{bmatrix}
!>                          \ell_{1,1}  &               &           &                   &               \\
!>                          \ell_{2,1}  & \ell_{2,2}    &           &                   &               \\
!>                          \ell_{3,1}  & \ell_{3,2}    & \ddots    &                   &               \\
!>                          \vdots      & \vdots        & \ddots    & \ell_{n-1,n-1}    &               \\
!>                          \ell_{n,1}  & \ell_{n,2}    & \ldots    & \ell_{n,n-1}      & \ell_{n,n}
!>                      \end{bmatrix} ~,
!>              \f}
!>              translates to the packed (column-major) format,
!>              \f{equation}{
!>                  L_{\mathrm{packed}} = [ \ell_{1,1}, \ell_{2,1}, \ell_{3,1}, \ldots, \ell_{n,1}, \ell_{2,2}, \ell_{3,2}, \ldots, \ell_{n,2}, \ldots, \ell_{n-1,n-1}, \ell_{n,n-1}, \ell_{n,n} ] ~.
!>              \f}
!>              In general, the indices \f$(i, j), 1\leq j\leq i \leq n\f$ of a full **lower-diagonal** \f$n\times n\f$ matrix translate to packed vector index,
!>              <ol>
!>                  <li>    \f$k(i, j) = i - 1 + (j - 1) * (2 * n - j) / 2\f$ for a column-major matrix layout (Fortran-style).
!>                  <li>    \f$k(i, j) = i - 1 + j * (j - 1) / 2\f$ for a row-major matrix layout (C-style).
!>              </ol>
!>              Similarly, an **upper-diagonal** triangular or Symmetric/Hermitian matrix of the form,
!>              \f{equation}{
!>                  U = \begin{bmatrix}
!>                          u_{1,1} &   u_{1,2} &   u_{1,3} &   \ldots  &   u_{1,n}     \\
!>                                  &   u_{2,2} &   u_{2,3} &   \ldots  &   u_{2,n}     \\
!>                                  &           &   u_{3,3} &   \ddots  &   \vdots      \\
!>                                  &           &           &   \ddots  &   u_{n-1,n}   \\
!>                                  &           &           &           &   u_{n,n}
!>                      \end{bmatrix} ~,
!>              \f}
!>              translates to the packed (column-major) format,
!>              \f{equation}{
!>                  U_{\mathrm{packed}} = [ u_{1,1}, u_{1,2}, u_{1,2}, u_{1,3}, u_{2,3}, u_{3,3}, \ldots, u_{1,n}, u_{2,n}, u_{3,2}, \ldots, u_{n-1,n}, u_{n,n} ] ~.
!>              \f}
!>              In general, the indices \f$(i, j), 1\leq i\leq j \leq n\f$ of a full **upper-diagonal** \f$n\times n\f$ matrix translate to the packed vector index,
!>              <ol>
!>                  <li>    \f$k(i, j) = i - 1 + j * (j - 1) / 2\f$ for a column-major matrix layout (Fortran-style).
!>                  <li>    \f$k(i, j) = j - 1 + (i - 1) * (2 * n - i) / 2\f$ for a row-major matrix layout (C-style).
!>              </ol>
!>      <li>    <b>Rectangular Band Packing (RBP or [rbpack](@ref pm_matrixPack))</b> which is to be added to this module.<br>
!>      <li>    <b>Rectangular Full Packing (RFP or [rfpack](@ref pm_matrixPack::rfpack))</b> for symmetric, Hermitian, or triangular matrices.<br>
!>              The Rectangular Full Packing is a combination of the contiguous and triangular packings.<br>
!>              It can be used to pack the upper or lower triangle of a symmetric, Hermitian, or triangular matrix **contiguously**.<br>
!>              It offers the package savings of the triangular packing plus the efficiency of using contiguous-pack Level 3 BLAS and LAPACK routines.<br>
!>              The RFP scheme is typically defined by three parameters in major libraries such as Intel MKL:<br>
!>              <ol>
!>                  <li>    A **matrix layout** parameter, which specifies column major (with the value `LAPACK_COL_MAJOR`) or row major (with the value `LAPACK_ROW_MAJOR`) matrix layout.<br>
!>                          By default, the Fortran routines of the ParaMonte library assume column-major layout, unless the procedure is to accessed from C-style programming languages.<br>
!>                  <li>    A [subset](@ref pm_matrixSubset::subset_type) parameter (represented by `uplo` in the LAPACK library),
!>                          which specifies which [upper](@ref pm_matrixSubset::uppDia)/[lower](@ref pm_matrixSubset::lowDia) triangle of the matrix is packed in RFP format.<br>
!>                          The upper/lower triangles are represented with the string values `"U"` and `"L"` respectively.<br>
!>                          <ol>
!>                              <li>    The following figure illustrates the RFP layout of an **even-order** matrix with [upper-diagonal triangle storage](@ref pm_matrixSubset::uppDia).<br>
!>                                      \image html pm_matrixPack@EURFP.png width=700
!>                              <li>    The following figure illustrates the RFP layout of an  **odd-order** matrix with [upper-diagonal triangle storage](@ref pm_matrixSubset::uppDia).<br>
!>                                      \image html pm_matrixPack@OURFP.png width=700
!>                              <li>    The following figure illustrates the RFP layout of an **even-order** matrix with [lower-diagonal triangle storage](@ref pm_matrixSubset::lowDia).<br>
!>                                      \image html pm_matrixPack@ELRFP.png width=700
!>                              <li>    The following figure illustrates the RFP layout of an  **odd-order** matrix with [lower-diagonal triangle storage](@ref pm_matrixSubset::lowDia).<br>
!>                                      \image html pm_matrixPack@OLRFP.png width=700
!>                          </ol>
!>                  <li>    A [transposition operation](@ref pm_matrixTrans) parameter (represented frequently by `transA` or `transB` in the LAPACK library),
!>                          which specifies one of the following three transposition operations on the matrix,
!>                          <ol>
!>                              <li>    [nothing](@ref pm_array::nothing) indicating No transposition (**as is**, with the LAPACK string value `"N"`).
!>                              <li>    [Symmetric transpose](@ref pm_matrixTrans::transSymm) (with the value `"T"`).
!>                              <li>    [Conjugate (Hermitian) transpose](@ref pm_matrixTrans::transHerm) (with the value `"C"`).
!>                          </ol>
!>              </ol>
!>              \note
!>              A matrix in RFP format, will always have an odd number of rows, unless it is transposed.<br>
!>              An odd-order  matrix has \f$\ms{ncol} = \ms{nrow} / 2 + 1\f$ columns.<br>
!>              An even-order matrix has \f$\ms{ncol} = \ms{nrow} / 2\f$ columns.<br>
!>      <li>    <b>Rectangular Default (i.e., standard, no specific) packing (RDP or [rdpack](@ref pm_matrixPack::rdpack))</b><br>
!>              This is the regular matrix package.<br>
!>              An \f$m\times n\f$ full matrix \f$A\f$ has the form,
!>              \f{equation}{
!>                  \large
!>                  \mathbf{A} = \begin{bmatrix}
!>                                  a_{11} & a_{12} & \cdots & a_{1n} \\
!>                                  a_{21} & a_{22} & \cdots & a_{2n} \\
!>                                  \vdots & \vdots & \ddots & \vdots \\
!>                                  a_{m1} & a_{m2} & \cdots & a_{mn}
!>                              \end{bmatrix} ~.
!>              \f}
!>
!>  </ol>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_matrixPack

    use pm_kind, only: SK, IK, LK
    use pm_matrixSubset, only: uppDia, uppDia_type

    implicit none

    character(*,SK), parameter :: MODULE_NAME = "@pm_matrixPack"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is an `abstract` derived type for constructing concrete derived types to
    !>  distinguish various procedure signatures that require different forms of matrix packing (triangular, Band, ...).<br>
    !>
    !>  \details
    !>  This `abstract` derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users must use `parameter` objects instantiated from the concrete subclasses of this parent `abstract` derived type.<br>
    !>
    !>  \see
    !>  [lfpack](@ref pm_matrixPack::lfpack)<br>
    !>  [rfpack](@ref pm_matrixPack::rfpack)<br>
    !>  [rdpack](@ref pm_matrixPack::rdpack)<br>
    !>  [lfpack_type](@ref pm_matrixPack::lfpack_type)<br>
    !>  [rfpack_type](@ref pm_matrixPack::rfpack_type)<br>
    !>  [lcpack_type](@ref pm_matrixPack::lcpack_type)<br>
    !>  [rcpack_type](@ref pm_matrixPack::rcpack_type)<br>
    !>  [ldpack_type](@ref pm_matrixPack::ldpack_type)<br>
    !>  [rdpack_type](@ref pm_matrixPack::rdpack_type)<br>
    !>  [package_type](@ref pm_matrixPack::package_type)<br>
    !>
    !>  \final{package_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, abstract :: package_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify <b>L</b>inear <b>S</b>parse (or <b>S</b>tandard) <b>Pack</b>ing format of a
    !>  (triangular) subset of a matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  or its subtypes as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [lfpack](@ref pm_matrixPack::lfpack)<br>
    !>  [rfpack](@ref pm_matrixPack::rfpack)<br>
    !>  [rdpack](@ref pm_matrixPack::rdpack)<br>
    !>  [lfpack_type](@ref pm_matrixPack::lfpack_type)<br>
    !>  [rfpack_type](@ref pm_matrixPack::rfpack_type)<br>
    !>  [lcpack_type](@ref pm_matrixPack::lcpack_type)<br>
    !>  [rcpack_type](@ref pm_matrixPack::rcpack_type)<br>
    !>  [ldpack_type](@ref pm_matrixPack::ldpack_type)<br>
    !>  [rdpack_type](@ref pm_matrixPack::rdpack_type)<br>
    !>  [package_type](@ref pm_matrixPack::package_type)<br>
    !>
    !>  \final{ldpack_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(package_type) :: ldpack_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify <b>R</b>ectangular <b>S</b>parse (or <b>S</b>tandard) <b>Pack</b>ing format of a
    !>  (triangular) subset of a matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [rdpack](@ref pm_matrixPack::rdpack)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [lfpack](@ref pm_matrixPack::lfpack)<br>
    !>  [rfpack](@ref pm_matrixPack::rfpack)<br>
    !>  [rdpack](@ref pm_matrixPack::rdpack)<br>
    !>  [lfpack_type](@ref pm_matrixPack::lfpack_type)<br>
    !>  [rfpack_type](@ref pm_matrixPack::rfpack_type)<br>
    !>  [lcpack_type](@ref pm_matrixPack::lcpack_type)<br>
    !>  [rcpack_type](@ref pm_matrixPack::rcpack_type)<br>
    !>  [ldpack_type](@ref pm_matrixPack::ldpack_type)<br>
    !>  [rdpack_type](@ref pm_matrixPack::rdpack_type)<br>
    !>  [package_type](@ref pm_matrixPack::package_type)<br>
    !>
    !>  \final{rdpack_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(package_type) :: rdpack_type
    end type

    !>  \brief
    !>  This is an object instance of class [rdpack_type](@ref pm_matrixPack::rdpack_type) that is exclusively used
    !>  to signify <b>R</b>ectangular <b>S</b>parse (or <b>S</b>tandard) <b>Pack</b>ing format of a given
    !>  (triangular) subset of a matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [lfpack](@ref pm_matrixPack::lfpack)<br>
    !>  [rfpack](@ref pm_matrixPack::rfpack)<br>
    !>  [rdpack](@ref pm_matrixPack::rdpack)<br>
    !>  [lfpack_type](@ref pm_matrixPack::lfpack_type)<br>
    !>  [rfpack_type](@ref pm_matrixPack::rfpack_type)<br>
    !>  [lcpack_type](@ref pm_matrixPack::lcpack_type)<br>
    !>  [rcpack_type](@ref pm_matrixPack::rcpack_type)<br>
    !>  [ldpack_type](@ref pm_matrixPack::ldpack_type)<br>
    !>  [rdpack_type](@ref pm_matrixPack::rdpack_type)<br>
    !>  [package_type](@ref pm_matrixPack::package_type)<br>
    !>
    !>  \final{rdpack}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(rdpack_type), parameter :: rdpack = rdpack_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: rdpack
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify <b>L</b>inear <b>C</b>ontiguous <b>Pack</b>ing format of a subset of a matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  sub-types of this parent type as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [lfpack](@ref pm_matrixPack::lfpack)<br>
    !>  [rfpack](@ref pm_matrixPack::rfpack)<br>
    !>  [rdpack](@ref pm_matrixPack::rdpack)<br>
    !>  [lfpack_type](@ref pm_matrixPack::lfpack_type)<br>
    !>  [rfpack_type](@ref pm_matrixPack::rfpack_type)<br>
    !>  [lcpack_type](@ref pm_matrixPack::lcpack_type)<br>
    !>  [rcpack_type](@ref pm_matrixPack::rcpack_type)<br>
    !>  [ldpack_type](@ref pm_matrixPack::ldpack_type)<br>
    !>  [rdpack_type](@ref pm_matrixPack::rdpack_type)<br>
    !>  [package_type](@ref pm_matrixPack::package_type)<br>
    !>
    !>  \final{lcpack_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(ldpack_type) :: lcpack_type
    end type

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify <b>R</b>ectangular <b>F</b>ull contiguous <b>Pack</b>ing format of a
    !>  (triangular) subset of a matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  sub-types of this parent type as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [lfpack](@ref pm_matrixPack::lfpack)<br>
    !>  [rfpack](@ref pm_matrixPack::rfpack)<br>
    !>  [rdpack](@ref pm_matrixPack::rdpack)<br>
    !>  [lfpack_type](@ref pm_matrixPack::lfpack_type)<br>
    !>  [rfpack_type](@ref pm_matrixPack::rfpack_type)<br>
    !>  [lcpack_type](@ref pm_matrixPack::lcpack_type)<br>
    !>  [rcpack_type](@ref pm_matrixPack::rcpack_type)<br>
    !>  [ldpack_type](@ref pm_matrixPack::ldpack_type)<br>
    !>  [rdpack_type](@ref pm_matrixPack::rdpack_type)<br>
    !>  [package_type](@ref pm_matrixPack::package_type)<br>
    !>
    !>  \final{rcpack_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(rdpack_type) :: rcpack_type
    end type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify <b>L</b>inear <b>F</b>ull contiguous <b>Pack</b>ing format of a
    !>  (triangular) subset of a matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [lfpack](@ref pm_matrixPack::lfpack)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [lfpack](@ref pm_matrixPack::lfpack)<br>
    !>  [rfpack](@ref pm_matrixPack::rfpack)<br>
    !>  [rdpack](@ref pm_matrixPack::rdpack)<br>
    !>  [lfpack_type](@ref pm_matrixPack::lfpack_type)<br>
    !>  [rfpack_type](@ref pm_matrixPack::rfpack_type)<br>
    !>  [lcpack_type](@ref pm_matrixPack::lcpack_type)<br>
    !>  [rcpack_type](@ref pm_matrixPack::rcpack_type)<br>
    !>  [ldpack_type](@ref pm_matrixPack::ldpack_type)<br>
    !>  [rdpack_type](@ref pm_matrixPack::rdpack_type)<br>
    !>  [package_type](@ref pm_matrixPack::package_type)<br>
    !>
    !>  \final{lfpack_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(lcpack_type) :: lfpack_type
    end type

    !>  \brief
    !>  This is an object instance of class [lfpack_type](@ref pm_matrixPack::lfpack_type) that is exclusively used
    !>  to signify <b>L</b>inear <b>F</b>ull contiguous <b>Pack</b>ing format of a
    !>  (triangular) subset of a matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [lfpack](@ref pm_matrixPack::lfpack)<br>
    !>  [rfpack](@ref pm_matrixPack::rfpack)<br>
    !>  [rdpack](@ref pm_matrixPack::rdpack)<br>
    !>  [lfpack_type](@ref pm_matrixPack::lfpack_type)<br>
    !>  [rfpack_type](@ref pm_matrixPack::rfpack_type)<br>
    !>  [lcpack_type](@ref pm_matrixPack::lcpack_type)<br>
    !>  [rcpack_type](@ref pm_matrixPack::rcpack_type)<br>
    !>  [ldpack_type](@ref pm_matrixPack::ldpack_type)<br>
    !>  [rdpack_type](@ref pm_matrixPack::rdpack_type)<br>
    !>  [package_type](@ref pm_matrixPack::package_type)<br>
    !>
    !>  \final{lfpack}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(lfpack_type), parameter :: lfpack = lfpack_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: lfpack
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  This is a concrete derived type whose instances are exclusively used
    !>  to signify <b>R</b>ectangular <b>F</b>ull contiguous <b>Pack</b>ing format of a
    !>  (triangular) subset of a matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  Objects instantiated from this derived type are exclusively used to differentiate
    !>  the procedures within the various generic interfaces of the ParaMonte library.<br>
    !>  As such, this concrete derived type does not contain any attributes.<br>
    !>
    !>  \note
    !>  This concrete derived type is not meant to be directly accessed by the end users.<br>
    !>  Instead, the end users should use the specific object parameter instance of this derived type
    !>  (e.g., [rfpack](@ref pm_matrixPack::rfpack)) as directed by the documentation of the specific procedure they intend to use.<br>
    !>
    !>  \see
    !>  [lfpack](@ref pm_matrixPack::lfpack)<br>
    !>  [rfpack](@ref pm_matrixPack::rfpack)<br>
    !>  [rdpack](@ref pm_matrixPack::rdpack)<br>
    !>  [lfpack_type](@ref pm_matrixPack::lfpack_type)<br>
    !>  [rfpack_type](@ref pm_matrixPack::rfpack_type)<br>
    !>  [lcpack_type](@ref pm_matrixPack::lcpack_type)<br>
    !>  [rcpack_type](@ref pm_matrixPack::rcpack_type)<br>
    !>  [ldpack_type](@ref pm_matrixPack::ldpack_type)<br>
    !>  [rdpack_type](@ref pm_matrixPack::rdpack_type)<br>
    !>  [package_type](@ref pm_matrixPack::package_type)<br>
    !>
    !>  \final{rfpack_type}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type, extends(rcpack_type) :: rfpack_type
    end type

    !>  \brief
    !>  This is an object instance of class [rfpack_type](@ref pm_matrixPack::rfpack_type) that is exclusively used
    !>  to signify <b>R</b>ectangular <b>F</b>ull contiguous <b>Pack</b>ing format of a
    !>  (triangular) subset of a matrix within an interface of a procedure of the ParaMonte library.<br>
    !>
    !>  \details
    !>  For example usage, see the documentation of the target procedure requiring this object.<br>
    !>
    !>  \see
    !>  [lfpack](@ref pm_matrixPack::lfpack)<br>
    !>  [rfpack](@ref pm_matrixPack::rfpack)<br>
    !>  [rdpack](@ref pm_matrixPack::rdpack)<br>
    !>  [lfpack_type](@ref pm_matrixPack::lfpack_type)<br>
    !>  [rfpack_type](@ref pm_matrixPack::rfpack_type)<br>
    !>  [lcpack_type](@ref pm_matrixPack::lcpack_type)<br>
    !>  [rcpack_type](@ref pm_matrixPack::rcpack_type)<br>
    !>  [ldpack_type](@ref pm_matrixPack::ldpack_type)<br>
    !>  [rdpack_type](@ref pm_matrixPack::rdpack_type)<br>
    !>  [package_type](@ref pm_matrixPack::package_type)<br>
    !>
    !>  \final{rfpack}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 12:00 AM, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin
    type(rfpack_type), parameter :: rfpack = rfpack_type()
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
    !DIR$ ATTRIBUTES DLLEXPORT :: rfpack
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return `.true.` if and only if a desired matrix with the specified input `shape` is of the specified input packing `pack`.<br>
    !>
    !>  \details
    !>  See [pm_matrixPack](@ref pm_matrixPack) for the mathematical definitions of different matrix packings.<br>
    !>  Note that a matrix in RFP format can be either in normal (long) or transpose (wide) format.<br>
    !>  This generic interface checks for both conditions to hold.<br>
    !>
    !>  \param[in]  pack    :   The input scalar constant that can be one of the following:<br>
    !>                          <ol>
    !>                              <li>    The scalar constant [lfpack](@ref pm_matrixPack::lfpack) implying the matrix has the Rectangular Full Packing (RFP) format.<br>
    !>                              <li>    The scalar constant [rfpack](@ref pm_matrixPack::rfpack) implying the matrix has the Rectangular Full Packing (RFP) format.<br>
    !>                          </ol>
    !>  \param[in]  shape   :   The input vector of length `2` of type `integer` of default kind \IK,
    !>                          representing the shape of the desired matrix whose specific packing is to be confirmed.<br>
    !>                          This input argument can be et to the direct output of the intrinsic Fortran function `shape()`.<br>
    !>
    !>  \return
    !>  `itis`              :   The output scalar `logical` of default kind \LK that is `.true.` if and only if
    !>                          the desired matrix with the specified input `shape` has the specified input packing `pack`.<br>
    !>
    !>  \interface{isMatPack}
    !>  \code{.F90}
    !>
    !>      use pm_matrixPack, only: isMatPack
    !>
    !>      itis = isMatPack(pack, shape)
    !>
    !>  \endcode
    !>
    !>  \pure
    !>
    !>  \example
    !>  \include{lineno} example/pm_matrixPack/isMatPack/main.F90
    !>  \compilef
    !>  \output
    !>  \include{lineno} example/pm_matrixPack/isMatPack/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixPack](@ref test_pm_matrixPack)
    !>
    !>  \pvhigh
    !>  This generic interface must be extended to `lfpack` format.
    !>
    !>  \final
    !>
    !>  \author
    !>  \AmirShahmoradi, Monday March 6, 2017, 3:22 pm, Institute for Computational Engineering and Sciences (ICES), The University of Texas at Austin.<br>
    interface isMatPack

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure module function isMatPackRFPD(pack, shape) result(itis)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: isMatPackRFPD
#endif
        type(rfpack_type)       , intent(in)                    :: pack
        integer(IK)             , intent(in)                    :: shape(2)
        logical(LK)                                             :: itis
    end function

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_matrixPack ! LCOV_EXCL_LINE