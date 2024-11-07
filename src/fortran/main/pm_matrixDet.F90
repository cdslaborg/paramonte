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
!>  This module contains procedures and generic interfaces relevant to the computation of the determinants of square matrices.<br>
!>
!>  \details
!>  The determinant is a scalar value that is a function of the entries of a square matrix.<br>
!>  It characterizes some properties of the matrix and the linear map represented by the matrix.<br>
!>  In particular, the determinant is nonzero if and only if the matrix is invertible and the linear map represented by the matrix is an isomorphism.<br>
!>  The determinant of a product of matrices is the product of their determinants (the preceding property is a corollary of this one).<br>
!>  The determinant of a matrix \f$A\f$ is denoted \f$\det(A)\f$, or \f$|A|\f$.<br>
!>  The determinant of a \f$2 \times 2\f$ matrix is,
!>  \f{equation}{
!>      \begin{vmatrix}a&b\\c&d\end{vmatrix} = ad - bc ~,
!>  \f}
!>  and the determinant of a \f$3 \times 3\f$ matrix is,
!>  \f{equation}{
!>      \begin{vmatrix}a&b&c\\d&e&f\\g&h&i\end{vmatrix} = aei + bfg + cdh - ceg - bdi - afh ~.
!>  \f}
!>  The determinant of an \f$n \times n\f$ matrix can be defined in several equivalent ways.<br>
!>  The **Leibniz formula** expresses the determinant as a sum of signed products of matrix entries such that each summand is the product of \f$n\f$ different entries,
!>  and the number of these summands is \f$n!\f$, the factorial of \f$n\f$ (the product of the \f$n\f$ first positive integers).<br>
!>  The Laplace expansion expresses the determinant of an \f$n\times n\f$ matrix as a linear combination of determinants of \f$(n-1)\times (n-1)\f$ submatrices.<br>
!>  Gaussian elimination expresses the determinant as the product of the diagonal entries of a diagonal matrix that is obtained by a succession of elementary row operations.<br>
!>  Determinants are used for defining the characteristic polynomial of a matrix, whose roots are the eigenvalues.<br>
!>  In geometry, the signed \f$n\f$-dimensional volume of a \f$n\f$-dimensional parallelepiped is expressed by a determinant.<br>
!>  This is used in calculus with exterior differential forms and the Jacobian determinant,
!>  in particular for changes of variables in multiple integrals.<br>
!>
!>  Properties of the determinant
!>  -----------------------------
!>
!>  The determinant is the **unique function** defined on the \f$n\times n\f$ matrices that has the four following properties:<br>
!>  <ol>
!>      <li>    The determinant of the identity matrix is 1.<br>
!>      <li>    The determinant is a homogeneous function, i.e.,<br>
!>              \f{equation}{
!>                  \det(cA) = c^{n}\det(A) ~,~ \text{for an} ~n\times n~ \text{matrix} ~A ~.
!>              \f}
!>      <li>    The exchange of two rows (or of two columns) multiplies the determinant by −1.<br>
!>      <li>    Multiplying a row (or a column) by a number multiplies the determinant by this number.<br>
!>      <li>    Adding to a row (or a column) a multiple of another row (or column) does not change the determinant.<br>
!>      <li>    If some column can be expressed as a linear combination of the other columns (i.e. the columns of the matrix form a linearly dependent set), the determinant is \f$0\f$.<br>
!>              As a special case, this includes: if some column is such that all its entries are zero, then the determinant of that matrix is \f$0\f$.<br>
!>      <li>    **Adding a scalar multiple of one column to another column does not change the value of the determinant**.<br>
!>              This is a consequence of multilinearity and being alternative:<br>
!>              By multilinearity the determinant changes by a multiple of the determinant of a matrix with two equal columns, which determinant is \f$0\f$, since the determinant is alternating.<br>
!>              If \f$A\f$ is a triangular matrix, i.e., \f$a_{ij} = 0\f$, whenever \f$i > j\f$ or, alternatively, whenever \f$i < j\f$,
!>              then its determinant equals the product of the diagonal entries:<br>
!>              \f{equation}{
!>                  \det(A) = a_{11}a_{22}\cdots a_{nn} = \prod _{i=1}^{n}a_{ii} ~.
!>              \f}
!>              Indeed, such a matrix can be reduced, by appropriately adding multiples of the columns with fewer nonzero entries to those with more entries,
!>              to a **diagonal matrix** (without changing the determinant).<br>
!>      <li>    The determinant of a general square matrix `mat` whose **LUP** factorization `lup(:,:)` is returned by [setMatLUP](@ref pm_matrixLUP::setMatLUP)
!>              can be computed as `(-1)**nswap * getMatMulTrace(lup)`, where `nswaps` is the number of swaps made to rows of the original `mat` compute its `lup`.<br>
!>              In other words, `nswap = count([(rperm(i) /= i, i = 1, size(rperm))]` where `rperm` is the row permutation vector returned by [setMatLUP](@ref pm_matrixLUP::setMatLUP).<br>
!>      <li>    The determinant of a positive-definite square matrix `mat` whose upper or lower Cholesky factorization `chol(:,:)` is returned by [setMatChol](@ref pm_matrixChol::setMatChol)
!>              is the square of the product of the diagonal elements of its Cholesky factorization and, by definition, is always a positive real number.<br>
!>  </ol>
!>
!>  \see
!>  [pm_matrixLUP](@ref pm_matrixLUP)<br>
!>  [pm_matrixChol](@ref pm_matrixChol)<br>
!>  [pm_matrixTrace](@ref pm_matrixTrace)<br>
!>  [pm_matrixTrans](@ref pm_matrixTrans)<br>
!>  [pm_matrixInv](@ref pm_matrixInv)<br>
!>  [Wikipedia Determinant](https://en.wikipedia.org/wiki/Determinant)<br>
!>  [Wolfram Determinant](https://mathworld.wolfram.com/Determinant.html)<br>
!>
!>  \test
!>  [test_pm_matrixDet](@ref test_pm_matrixDet)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Oct 18, 2009, 4:10 PM, Houghton, Michigan

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_matrixDet

    use pm_kind, only: SK, IK, LK
    use pm_matrixSubset, only: subset_type
    use pm_matrixSubset, only: uppDia, uppDia_type
    use pm_matrixSubset, only: lowDia, lowDia_type
    use pm_matrixTrans, only: transHerm, transHerm_type
    use pm_array, only: nothing, nothing_type

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_matrixDet"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the determinant of the input general **square** matrix.<br>
    !>
    !>  \details
    !>  This generic interface uses the [pivoted LU (LUP) factorization](@ref pm_matrixLUP) to compute the determinant.<br>
    !>  If the matrix is already known to be positive-definite, use the other more appropriate faster generic
    !>  interfaces of this module (e.g., [getMatDetSqrt](@ref pm_matrixDet::getMatDetSqrt)).<br>
    !>
    !>  \param[in]  mat :   The input `contiguous` square matrix of shape `(ndim,ndim)` of type `real` of kind \RKALL.
    !>
    !>  \return
    !>  `det`           :   The output scalar of the same type and kind as the input `mat`
    !>                      containing the determinant of the input square matrix.<br>
    !>
    !>  \interface{getMatDet}
    !>  \code{.F90}
    !>
    !>      use pm_matrixDet, only: getMatDet
    !>
    !>      det = getMatDet(mat)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(mat, 1) == size(mat, 2)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warning
    !>  The procedures under this generic interface call `error stop` if the Pivoted LU factorization of `mat` fails.<br>
    !>  See [setMatDet](@ref pm_matrixDet::setMatDet) for the fault-tolerant potentially faster subroutine version of this interface.<br>
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The procedures under this generic interface are computationally demanding since the LU factorization requires
    !>  creating a copy of the input `mat` to avoid the direct manipulation of the input `intent(in) :: mat`.<br>
    !>  See [setMatDet](@ref pm_matrixDet::setMatDet) for the performant subroutine versions of these procedures
    !>  where `mat` has `intent(inout)`.<br>
    !>
    !>  \see
    !>  [getMatDet](@ref pm_matrixDet::getMatDet)<br>
    !>  [setMatDet](@ref pm_matrixDet::setMatDet)<br>
    !>  [getMatDetSqrtLog](@ref pm_matrixDet::getMatDetSqrtLog)<br>
    !>  [setMatDetSqrtLog](@ref pm_matrixDet::setMatDetSqrtLog)<br>
    !>  [getMatDetSqrt](@ref pm_matrixDet::getMatDetSqrt)<br>
    !>  [setMatDetSqrt](@ref pm_matrixDet::setMatDetSqrt)<br>
    !>
    !>  \example{getMatDet}
    !>  \include{lineno} example/pm_matrixDet/getMatDet/main.F90
    !>  \compilef{getMatDet}
    !>  \output{getMatDet}
    !>  \include{lineno} example/pm_matrixDet/getMatDet/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixDet](@ref test_pm_matrixDet)
    !>
    !>  \bug
    !>  \status \unresolved
    !>  \source \gfortran{11.2}
    !>  \desc
    !>  The \gfortran{11.2} cannot compile the following code,
    !>  \code{.F90}
    !>
    !>      interface
    !>      function test(mat) result(res)
    !>          real, value :: mat(:,:)
    !>          real :: res
    !>      end function
    !>      end interface
    !>      end
    !>
    !>  \endcode
    !>  yielding the following error message,
    !>  \code{bash}
    !>      Error: VALUE attribute conflicts with DIMENSION attribute at (1)
    !>  \endcode
    !>  This error appears to conform with the Fortran 2003 standard but not with the newer standards.<br>
    !>  The \ifort has no problem compiling this code.<br>
    !>  \remedy
    !>  Avoid such interfaces until the \gfortran bug is resolved.
    !>
    !>  \todo
    !>  \pvhigh
    !>  As soon as the Gfortran bug described above is resolved, the `mat` argument of the
    !>  procedures under this generic interface can be changed to have the `value` attribute instead of the `intent(in)`.<br>
    !>  This will obviate the need for creating an additional copy of the array inside the routine, which has to be also modified subsequently.<br>
    !>
    !>  \todo
    !>  \plow
    !>  This generic interface could be extended to accept different classes of matrices , for example,
    !>  <ol>
    !>      <li>    [posdefmat_type](@ref pm_matrixClass::posdefmat_type),
    !>      <li>    [upperDiag_type](@ref pm_matrixClass::upperDiag_type),
    !>      <li>    [lowerDiag_type](@ref pm_matrixClass::lowerDiag_type),
    !>      <li>    [upperUnit_type](@ref pm_matrixClass::upperUnit_type),
    !>      <li>    [lowerUnit_type](@ref pm_matrixClass::lowerUnit_type),
    !>  </ol>
    !>  and other types to facilitate runtime dispatch to the appropriate routines.<br>
    !>  However, class-specific routines were deemed unnecessary of this version of the library as
    !>  the task of computing the determinant of such matrices is either trivially explained in the
    !>  [module documentation](@ref pm_matrixDet) or already implemented by other generic interfaces of this module.<br>
    !>
    !>  \final{getMatDet}
    !>
    !>  \author
    !>  \AmirShahmoradi, Apr 21, 2017, 1:43 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface getMatDet

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatDet_CK5(mat) result(det)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDet_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(in), contiguous    :: mat(:,:)
        complex(CKG)                                    :: det
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatDet_CK4(mat) result(det)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDet_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(in), contiguous    :: mat(:,:)
        complex(CKG)                                    :: det
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatDet_CK3(mat) result(det)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDet_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(in), contiguous    :: mat(:,:)
        complex(CKG)                                    :: det
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatDet_CK2(mat) result(det)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDet_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(in), contiguous    :: mat(:,:)
        complex(CKG)                                    :: det
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatDet_CK1(mat) result(det)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDet_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(in), contiguous    :: mat(:,:)
        complex(CKG)                                    :: det
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatDet_RK5(mat) result(det)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDet_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(in), contiguous    :: mat(:,:)
        real(RKG)                                       :: det
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatDet_RK4(mat) result(det)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDet_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(in), contiguous    :: mat(:,:)
        real(RKG)                                       :: det
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatDet_RK3(mat) result(det)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDet_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(in), contiguous    :: mat(:,:)
        real(RKG)                                       :: det
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatDet_RK2(mat) result(det)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDet_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(in), contiguous    :: mat(:,:)
        real(RKG)                                       :: det
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatDet_RK1(mat) result(det)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDet_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(in), contiguous    :: mat(:,:)
        real(RKG)                                       :: det
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the determinant of the input **square** matrix.
    !>
    !>  \details
    !>  This generic interface uses the [pivoted LU (LUP) factorization](@ref pm_matrixLUP) to compute the determinant.<br>
    !>  If the matrix is already known to be positive-definite, use the other more appropriate faster generic
    !>  interfaces of this module (e.g., [setMatDetSqrt](@ref pm_matrixDet::setMatDetSqrt)).<br>
    !>
    !>  \param[inout]   mat     :   The input/output `contiguous` square matrix of shape `(1:ndim, 1:ndim)` of,
    !>                              <ol>
    !>                                  <li>    type `complex` of kind \CKALL,
    !>                                  <li>    type `real` of kind \RKALL,
    !>                              </ol>
    !>                              On input, it represents the general matrix whose determinant must be computed.<br>
    !>                              On output, it is completely overwritten by its [LUP factorization](@ref pm_matrixLUP).<br>
    !>  \param[out]     det     :   The output scalar of the same type and kind as the input `mat` representing its determinant.
    !>  \param[out]     info    :   The output scalar of type `integer` of default kind \IK that is non-zero **if and only if**
    !>                              the Pivoted LU-factorization of the input `mat` fails.<br>
    !>                              See the corresponding `info` argument of [setMatLUP](@ref pm_matrixLUP::setMatLUP)
    !>                              for possible meanings of the output non-zero error code.<br>
    !>
    !>  \interface{setMatDet}
    !>  \code{.F90}
    !>
    !>      use pm_matrixDet, only: setMatDet
    !>
    !>      call setMatDet(mat, det, info)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(mat, 1) == size(mat, 2)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getMatDet](@ref pm_matrixDet::getMatDet)<br>
    !>  [setMatDet](@ref pm_matrixDet::setMatDet)<br>
    !>  [getMatDetSqrtLog](@ref pm_matrixDet::getMatDetSqrtLog)<br>
    !>  [setMatDetSqrtLog](@ref pm_matrixDet::setMatDetSqrtLog)<br>
    !>  [getMatDetSqrt](@ref pm_matrixDet::getMatDetSqrt)<br>
    !>  [setMatDetSqrt](@ref pm_matrixDet::setMatDetSqrt)<br>
    !>
    !>  \example{setMatDet}
    !>  \include{lineno} example/pm_matrixDet/setMatDet/main.F90
    !>  \compilef{setMatDet}
    !>  \output{setMatDet}
    !>  \include{lineno} example/pm_matrixDet/setMatDet/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixDet](@ref test_pm_matrixDet)
    !>
    !>  \final{setMatDet}
    !>
    !>  \author
    !>  \AmirShahmoradi, Apr 21, 2017, 1:43 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    interface setMatDet

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatDet_CK5(mat, det, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDet_CK5
#endif
        use pm_kind, only: CKG => CK5
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        complex(CKG)        , intent(out)                   :: det
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatDet_CK4(mat, det, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDet_CK4
#endif
        use pm_kind, only: CKG => CK4
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        complex(CKG)        , intent(out)                   :: det
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatDet_CK3(mat, det, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDet_CK3
#endif
        use pm_kind, only: CKG => CK3
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        complex(CKG)        , intent(out)                   :: det
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatDet_CK2(mat, det, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDet_CK2
#endif
        use pm_kind, only: CKG => CK2
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        complex(CKG)        , intent(out)                   :: det
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatDet_CK1(mat, det, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDet_CK1
#endif
        use pm_kind, only: CKG => CK1
        complex(CKG)        , intent(inout) , contiguous    :: mat(:,:)
        complex(CKG)        , intent(out)                   :: det
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatDet_RK5(mat, det, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDet_RK5
#endif
        use pm_kind, only: RKG => RK5
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        real(RKG)           , intent(out)                   :: det
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatDet_RK4(mat, det, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDet_RK4
#endif
        use pm_kind, only: RKG => RK4
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        real(RKG)           , intent(out)                   :: det
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatDet_RK3(mat, det, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDet_RK3
#endif
        use pm_kind, only: RKG => RK3
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        real(RKG)           , intent(out)                   :: det
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatDet_RK2(mat, det, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDet_RK2
#endif
        use pm_kind, only: RKG => RK2
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        real(RKG)           , intent(out)                   :: det
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatDet_RK1(mat, det, info)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDet_RK1
#endif
        use pm_kind, only: RKG => RK1
        real(RKG)           , intent(inout) , contiguous    :: mat(:,:)
        real(RKG)           , intent(out)                   :: det
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the square-root of the determinant of the input **positive-definite** square matrix.
    !>
    !>  \details
    !>  The positive-definiteness guarantee by the user enables a fast method of computing the
    !>  determinant of the input positive-definite square matrix via its Cholesky Factorization.<br>
    !>  For high rank matrices, there is an overflow possibility for the output of this generic interface.<br>
    !>  The overflow risk can be avoided by calling [getMatDetSqrtLog](@ref pm_matrixDet::getMatDetSqrtLog) instead.<br>
    !>
    !>  \param[in]      mat     :   The input `contiguous` positive-definite square matrix of shape `(ndim,ndim)` of type `real` of kind \RKALL.
    !>                              Only the upper triangle and diagonals of the matrix are needed and used to compute the Cholesky Factorization.
    !>  \param[in]      subset  :   The input scalar that can be either,
    !>                              <ol>
    !>                                  <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) implying that the upper-diagonal triangular block
    !>                                          of `mat` should be used for building the **lower-diagonal** Cholesky factor in the output `chol`.<br>
    !>                                  <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) implying that the lower-diagonal triangular block
    !>                                          of `mat` should be used for building the **upper-diagonal** Cholesky factor in the output `chol`.<br>
    !>                              </ol>
    !>                              This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.<br>
    !>                              (**optional**, default = [uppDia](@ref pm_matrixSubset::uppDia))
    !>
    !>  \return
    !>  `detSqrtLog`            :   The output scalar of type `real` of the same kind as the input `mat` containing
    !>                              the square-root of the determinant of the input positive-definite square matrix.<br>
    !>                              If the input `mat` is not positive-definite, the algorithm will halt the program by calling `error stop`.<br>
    !>
    !>  \interface{getMatDetSqrtLog}
    !>  \code{.F90}
    !>
    !>      use pm_matrixDet, only: getMatDetSqrtLog
    !>
    !>      detSqrtLog = getMatDetSqrtLog(mat, subset = subset)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(mat, 1) == size(mat, 2)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warning
    !>  The input `mat` is assumed to be positive-definite, otherwise the Cholesky Factorization within the procedure will fail,
    !>  leading to the invocation of the `error stop` Fortran statement.<br>
    !>  See [setMatDetSqrtLog](@ref pm_matrixDet::setMatDetSqrtLog) for the faster fault-tolerant subroutine versions of these procedures.
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The procedures under this generic interface are potentially computationally demanding compared to [setMatDetSqrtLog](@ref pm_matrixDet::setMatDetSqrtLog).<br>
    !>  See [setMatDetSqrtLog](@ref pm_matrixDet::setMatDetSqrtLog) for a more performant subroutine version of this generic interface for repeated calls.<br>
    !>
    !>  \see
    !>  [getMatDet](@ref pm_matrixDet::getMatDet)<br>
    !>  [setMatDet](@ref pm_matrixDet::setMatDet)<br>
    !>  [getMatDetSqrt](@ref pm_matrixDet::getMatDetSqrt)<br>
    !>  [setMatDetSqrt](@ref pm_matrixDet::setMatDetSqrt)<br>
    !>  [getMatDetSqrtLog](@ref pm_matrixDet::getMatDetSqrtLog)<br>
    !>  [setMatDetSqrtLog](@ref pm_matrixDet::setMatDetSqrtLog)<br>
    !>  [getMatMulTrace](@ref pm_matrixTrace::getMatMulTrace)<br>
    !>
    !>  \example{getMatDetSqrtLog}
    !>  \include{lineno} example/pm_matrixDet/getMatDetSqrtLog/main.F90
    !>  \compilef{getMatDetSqrtLog}
    !>  \output{getMatDetSqrtLog}
    !>  \include{lineno} example/pm_matrixDet/getMatDetSqrtLog/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixDet](@ref test_pm_matrixDet)
    !>
    !>  \final{getMatDetSqrtLog}
    !>
    !>  \author
    !>  \AmirShahmoradi, Apr 21, 2017, 1:43 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    interface getMatDetSqrt

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatDetSqrt_CK5(mat, subset) result(detSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrt_CK5
#endif
        use pm_kind, only: CKG => CK5
        class(subset_type)  , intent(in), optional      :: subset
        complex(CKG)        , intent(in), contiguous    :: mat(:,:)
        real(CKG)                                       :: detSqrt
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatDetSqrt_CK4(mat, subset) result(detSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrt_CK4
#endif
        use pm_kind, only: CKG => CK4
        class(subset_type)  , intent(in), optional      :: subset
        complex(CKG)        , intent(in), contiguous    :: mat(:,:)
        real(CKG)                                       :: detSqrt
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatDetSqrt_CK3(mat, subset) result(detSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrt_CK3
#endif
        use pm_kind, only: CKG => CK3
        class(subset_type)  , intent(in), optional      :: subset
        complex(CKG)        , intent(in), contiguous    :: mat(:,:)
        real(CKG)                                       :: detSqrt
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatDetSqrt_CK2(mat, subset) result(detSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrt_CK2
#endif
        use pm_kind, only: CKG => CK2
        class(subset_type)  , intent(in), optional      :: subset
        complex(CKG)        , intent(in), contiguous    :: mat(:,:)
        real(CKG)                                       :: detSqrt
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatDetSqrt_CK1(mat, subset) result(detSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrt_CK1
#endif
        use pm_kind, only: CKG => CK1
        class(subset_type)  , intent(in), optional      :: subset
        complex(CKG)        , intent(in), contiguous    :: mat(:,:)
        real(CKG)                                       :: detSqrt
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatDetSqrt_RK5(mat, subset) result(detSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrt_RK5
#endif
        use pm_kind, only: RKG => RK5
        class(subset_type)  , intent(in), optional      :: subset
        real(RKG)           , intent(in), contiguous    :: mat(:,:)
        real(RKG)                                       :: detSqrt
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatDetSqrt_RK4(mat, subset) result(detSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrt_RK4
#endif
        use pm_kind, only: RKG => RK4
        class(subset_type)  , intent(in), optional      :: subset
        real(RKG)           , intent(in), contiguous    :: mat(:,:)
        real(RKG)                                       :: detSqrt
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatDetSqrt_RK3(mat, subset) result(detSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrt_RK3
#endif
        use pm_kind, only: RKG => RK3
        class(subset_type)  , intent(in), optional      :: subset
        real(RKG)           , intent(in), contiguous    :: mat(:,:)
        real(RKG)                                       :: detSqrt
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatDetSqrt_RK2(mat, subset) result(detSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrt_RK2
#endif
        use pm_kind, only: RKG => RK2
        class(subset_type)  , intent(in), optional      :: subset
        real(RKG)           , intent(in), contiguous    :: mat(:,:)
        real(RKG)                                       :: detSqrt
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatDetSqrt_RK1(mat, subset) result(detSqrt)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrt_RK1
#endif
        use pm_kind, only: RKG => RK1
        class(subset_type)  , intent(in), optional      :: subset
        real(RKG)           , intent(in), contiguous    :: mat(:,:)
        real(RKG)                                       :: detSqrt
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the determinant of the input **positive-definite** square matrix.
    !>
    !>  \details
    !>  The positive-definiteness guarantee by the user enables a fast method of computing the
    !>  square root of the determinant of the input positive-definite square matrix via its Cholesky Factorization.<br>
    !>  For high rank matrices, there is an overflow possibility for the output of this generic interface.<br>
    !>  The overflow risk can be avoided by calling [setMatDetSqrtLog](@ref pm_matrixDet::setMatDetSqrtLog) instead.<br>
    !>
    !>  \note
    !>  This generic interface is meant to facilitate simultaneous computation of
    !>  the Cholesky factorization and determinant of a input positive definite matrix.<br>
    !>  If the Cholesky factorization is already computed, simply compute the determinant
    !>  as the square of the [multiplicative trace](@ref pm_matrixTrace) of the Cholesky factor.<br>
    !>
    !>  \param[in]      mat         :   The input `contiguous` matrix of shape `(1:ndim, 1:ndim)` of
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing a `subset` of the positive-definite matrix whose `log(sqrt(determinant))` is to be computed.<br>
    !>  \param[in]      subset      :   The input scalar that can be either,
    !>                                  <ol>
    !>                                      <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) implying that the upper-diagonal triangular block
    !>                                              of `mat` should be used for building the **lower-diagonal** Cholesky factor in the output `chol`.<br>
    !>                                      <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) implying that the lower-diagonal triangular block
    !>                                              of `mat` should be used for building the **upper-diagonal** Cholesky factor in the output `chol`.<br>
    !>                                  </ol>
    !>                                  This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.
    !>  \param[out]     detSqrt     :   The output scalar of type `real` of the same kind as the input `mat` representing the square root of its determinant.<br>
    !>  \param[out]     info        :   The output scalar `integer` of default kind \IK that is `0` <b>if and only if</b> the Cholesky factorization succeeds.<br>
    !>                                  Otherwise, it is set to the order of the leading minor of the specified input subset of `mat` that is not positive definite,
    !>                                  indicating the occurrence of an error and that the factorization could not be completed.<br>
    !>                                  **A non-zero value implies failure of the Cholesky factorization.**<br>
    !>                                  A non-zero `info` implied the input `mat` is not positive definite.<br>
    !>  \param[inout]   chol        :   The input/output matrix of the same type, kind, and shape as the input `mat` containing
    !>                                  the Cholesky factorization in its triangular subset dictated by the input arguments `subset` and `operation`.<br>
    !>                                  See the documentation of the input argument `mat` above for more information.<br>
    !>  \param[in]      operation   :   The input scalar that can be either,
    !>                                  <ol>
    !>                                      <li>    the constant [nothing](@ref pm_array::nothing) implying that
    !>                                              the same subset of the output `chol` as the input `subset` must be populated with the Cholesky factorization.<br>
    !>                                      <li>    the constant [transHerm](@ref pm_matrixTrans::transHerm) implying that
    !>                                              the complex transpose of the computed Cholesky factorization must be written to the output `chol`.<br>
    !>                                  </ol>
    !>
    !>  \interface{setMatDetSqrtLog}
    !>  \code{.F90}
    !>
    !>      use pm_matrixDet, only: setMatDetSqrtLog
    !>
    !>      call setMatDetSqrtLog(mat, subset, detSqrt, info, chol, operation)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(mat, 1) == size(mat, 2)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warning
    !>  The input `mat` is assumed to be positive-definite, otherwise the Cholesky Factorization within the procedure will fail,
    !>  in which case, the procedure will return the control with `info /= 0`.<br>
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getMatDet](@ref pm_matrixDet::getMatDet)<br>
    !>  [setMatDet](@ref pm_matrixDet::setMatDet)<br>
    !>  [getMatDetSqrt](@ref pm_matrixDet::getMatDetSqrt)<br>
    !>  [setMatDetSqrt](@ref pm_matrixDet::setMatDetSqrt)<br>
    !>  [getMatDetSqrtLog](@ref pm_matrixDet::getMatDetSqrtLog)<br>
    !>  [setMatDetSqrtLog](@ref pm_matrixDet::setMatDetSqrtLog)<br>
    !>  [getMatMulTrace](@ref pm_matrixTrace::getMatMulTrace)<br>
    !>
    !>  \example{setMatDetSqrtLog}
    !>  \include{lineno} example/pm_matrixDet/setMatDetSqrtLog/main.F90
    !>  \compilef{setMatDetSqrtLog}
    !>  \output{setMatDetSqrtLog}
    !>  \include{lineno} example/pm_matrixDet/setMatDetSqrtLog/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixDet](@ref test_pm_matrixDet)
    !>
    !>  \final{setMatDetSqrtLog}
    !>
    !>  \author
    !>  \AmirShahmoradi, Apr 21, 2017, 1:43 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! nothing

    interface setMatDetSqrt

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_ONO_CK5(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_ONO_CK4(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_ONO_CK3(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_ONO_CK2(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_ONO_CK1(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_ONO_RK5(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_ONO_RK4(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_ONO_RK3(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_ONO_RK2(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_ONO_RK1(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_ONO_CK5(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_ONO_CK4(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_ONO_CK3(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_ONO_CK2(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_ONO_CK1(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_ONO_RK5(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_ONO_RK4(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_ONO_RK3(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_ONO_RK2(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_ONO_RK1(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! transHerm

    interface setMatDetSqrt

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_OTH_CK5(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_OTH_CK5
#endif
        use pm_kind, only: CKG => CK5
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_OTH_CK4(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_OTH_CK4
#endif
        use pm_kind, only: CKG => CK4
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_OTH_CK3(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_OTH_CK3
#endif
        use pm_kind, only: CKG => CK3
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_OTH_CK2(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_OTH_CK2
#endif
        use pm_kind, only: CKG => CK2
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_OTH_CK1(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_OTH_CK1
#endif
        use pm_kind, only: CKG => CK1
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_OTH_RK5(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_OTH_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_OTH_RK4(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_OTH_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_OTH_RK3(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_OTH_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_OTH_RK2(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_OTH_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatDetSqrt_UXD_OTH_RK1(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_UXD_OTH_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_OTH_CK5(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_OTH_CK5
#endif
        use pm_kind, only: CKG => CK5
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_OTH_CK4(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_OTH_CK4
#endif
        use pm_kind, only: CKG => CK4
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_OTH_CK3(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_OTH_CK3
#endif
        use pm_kind, only: CKG => CK3
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_OTH_CK2(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_OTH_CK2
#endif
        use pm_kind, only: CKG => CK2
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_OTH_CK1(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_OTH_CK1
#endif
        use pm_kind, only: CKG => CK1
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_OTH_RK5(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_OTH_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_OTH_RK4(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_OTH_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_OTH_RK3(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_OTH_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_OTH_RK2(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_OTH_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatDetSqrt_XLD_OTH_RK1(mat, subset, detSqrt, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrt_XLD_OTH_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrt
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Generate and return the natural logarithm of the square-root of the determinant of the input **positive-definite** square matrix.
    !>
    !>  \details
    !>  The positive-definiteness guarantee by the user enables a fast method of computing the
    !>  determinant of the input positive-definite square matrix via its Cholesky Factorization.
    !>
    !>  \param[in]      mat     :   The input `contiguous` positive-definite square matrix of shape `(ndim,ndim)` of type `real` of kind \RKALL.
    !>                              Only the upper triangle and diagonals of the matrix are needed and used to compute the Cholesky Factorization.
    !>  \param[in]      subset  :   The input scalar that can be either,
    !>                              <ol>
    !>                                  <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) implying that the upper-diagonal triangular block
    !>                                          of `mat` should be used for building the **lower-diagonal** Cholesky factor in the output `chol`.<br>
    !>                                  <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) implying that the lower-diagonal triangular block
    !>                                          of `mat` should be used for building the **upper-diagonal** Cholesky factor in the output `chol`.<br>
    !>                              </ol>
    !>                              This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.
    !>                              (**optional**, default = [uppDia](@ref pm_matrixSubset::uppDia))
    !>
    !>  \return
    !>  `detSqrtLog`            :   The output scalar of the type `real` of the same kind as the input `mat` containing the
    !>                              natural logarithm of the square-root of the determinant of the input positive-definite square matrix.<br>
    !>                              If the input `mat` is not positive-definite, the algorithm will halt the program by calling `error stop`.<br>
    !>
    !>  \interface{getMatDetSqrtLog}
    !>  \code{.F90}
    !>
    !>      use pm_matrixDet, only: getMatDetSqrtLog
    !>
    !>      detSqrtLog = getMatDetSqrtLog(mat, subset = subset)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(mat, 1) == size(mat, 2)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warning
    !>  The input `mat` is assumed to be positive-definite, otherwise the Cholesky Factorization within the procedure will fail,
    !>  leading to the invocation of the `error stop` Fortran statement.<br>
    !>  See [setMatDetSqrtLog](@ref pm_matrixDet::setMatDetSqrtLog) for the faster fault-tolerant subroutine versions of these procedures.
    !>
    !>  \warnpure
    !>
    !>  \remark
    !>  The procedures under this generic interface are potentially computationally demanding compared to [setMatDetSqrtLog](@ref pm_matrixDet::setMatDetSqrtLog).<br>
    !>  See [setMatDetSqrtLog](@ref pm_matrixDet::setMatDetSqrtLog) for a more performant subroutine version of this generic interface for repeated calls.<br>
    !>
    !>  \see
    !>  [getMatDet](@ref pm_matrixDet::getMatDet)<br>
    !>  [setMatDet](@ref pm_matrixDet::setMatDet)<br>
    !>  [getMatDetSqrt](@ref pm_matrixDet::getMatDetSqrt)<br>
    !>  [setMatDetSqrt](@ref pm_matrixDet::setMatDetSqrt)<br>
    !>  [getMatDetSqrtLog](@ref pm_matrixDet::getMatDetSqrtLog)<br>
    !>  [setMatDetSqrtLog](@ref pm_matrixDet::setMatDetSqrtLog)<br>
    !>  [getMatMulTrace](@ref pm_matrixTrace::getMatMulTrace)<br>
    !>
    !>  \example{getMatDetSqrtLog}
    !>  \include{lineno} example/pm_matrixDet/getMatDetSqrtLog/main.F90
    !>  \compilef{getMatDetSqrtLog}
    !>  \output{getMatDetSqrtLog}
    !>  \include{lineno} example/pm_matrixDet/getMatDetSqrtLog/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixDet](@ref test_pm_matrixDet)
    !>
    !>  \final{getMatDetSqrtLog}
    !>
    !>  \author
    !>  \AmirShahmoradi, Apr 21, 2017, 1:43 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    interface getMatDetSqrtLog

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module function getMatDetSqrtLog_CK5(mat, subset) result(detSqrtLog)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrtLog_CK5
#endif
        use pm_kind, only: CKG => CK5
        class(subset_type)  , intent(in), optional      :: subset
        complex(CKG)        , intent(in), contiguous    :: mat(:,:)
        real(CKG)                                       :: detSqrtLog
    end function
#endif

#if CK4_ENABLED
    PURE module function getMatDetSqrtLog_CK4(mat, subset) result(detSqrtLog)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrtLog_CK4
#endif
        use pm_kind, only: CKG => CK4
        class(subset_type)  , intent(in), optional      :: subset
        complex(CKG)        , intent(in), contiguous    :: mat(:,:)
        real(CKG)                                       :: detSqrtLog
    end function
#endif

#if CK3_ENABLED
    PURE module function getMatDetSqrtLog_CK3(mat, subset) result(detSqrtLog)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrtLog_CK3
#endif
        use pm_kind, only: CKG => CK3
        class(subset_type)  , intent(in), optional      :: subset
        complex(CKG)        , intent(in), contiguous    :: mat(:,:)
        real(CKG)                                       :: detSqrtLog
    end function
#endif

#if CK2_ENABLED
    PURE module function getMatDetSqrtLog_CK2(mat, subset) result(detSqrtLog)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrtLog_CK2
#endif
        use pm_kind, only: CKG => CK2
        class(subset_type)  , intent(in), optional      :: subset
        complex(CKG)        , intent(in), contiguous    :: mat(:,:)
        real(CKG)                                       :: detSqrtLog
    end function
#endif

#if CK1_ENABLED
    PURE module function getMatDetSqrtLog_CK1(mat, subset) result(detSqrtLog)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrtLog_CK1
#endif
        use pm_kind, only: CKG => CK1
        class(subset_type)  , intent(in), optional      :: subset
        complex(CKG)        , intent(in), contiguous    :: mat(:,:)
        real(CKG)                                       :: detSqrtLog
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module function getMatDetSqrtLog_RK5(mat, subset) result(detSqrtLog)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrtLog_RK5
#endif
        use pm_kind, only: RKG => RK5
        class(subset_type)  , intent(in), optional      :: subset
        real(RKG)           , intent(in), contiguous    :: mat(:,:)
        real(RKG)                                       :: detSqrtLog
    end function
#endif

#if RK4_ENABLED
    PURE module function getMatDetSqrtLog_RK4(mat, subset) result(detSqrtLog)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrtLog_RK4
#endif
        use pm_kind, only: RKG => RK4
        class(subset_type)  , intent(in), optional      :: subset
        real(RKG)           , intent(in), contiguous    :: mat(:,:)
        real(RKG)                                       :: detSqrtLog
    end function
#endif

#if RK3_ENABLED
    PURE module function getMatDetSqrtLog_RK3(mat, subset) result(detSqrtLog)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrtLog_RK3
#endif
        use pm_kind, only: RKG => RK3
        class(subset_type)  , intent(in), optional      :: subset
        real(RKG)           , intent(in), contiguous    :: mat(:,:)
        real(RKG)                                       :: detSqrtLog
    end function
#endif

#if RK2_ENABLED
    PURE module function getMatDetSqrtLog_RK2(mat, subset) result(detSqrtLog)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrtLog_RK2
#endif
        use pm_kind, only: RKG => RK2
        class(subset_type)  , intent(in), optional      :: subset
        real(RKG)           , intent(in), contiguous    :: mat(:,:)
        real(RKG)                                       :: detSqrtLog
    end function
#endif

#if RK1_ENABLED
    PURE module function getMatDetSqrtLog_RK1(mat, subset) result(detSqrtLog)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMatDetSqrtLog_RK1
#endif
        use pm_kind, only: RKG => RK1
        class(subset_type)  , intent(in), optional      :: subset
        real(RKG)           , intent(in), contiguous    :: mat(:,:)
        real(RKG)                                       :: detSqrtLog
    end function
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the natural logarithm of the square-root of the determinant of the input **positive-definite** square matrix.
    !>
    !>  \details
    !>  The positive-definiteness guarantee by the user enables a fast method of computing the
    !>  determinant of the input positive-definite square matrix via its Cholesky Factorization.<br>
    !>  Computing the logarithm of the square-root of the determinant is preferrable for high rank matrices where there is a possibility of overflow,
    !>  although it comes with the additional computational cost of computing `size(mat, 1)` extra logarithms.<br>
    !>
    !>  \note
    !>  This generic interface is meant to facilitate simultaneous computation of
    !>  the Cholesky factorization and determinant of a input positive definite matrix.<br>
    !>  If the Cholesky factorization is already computed, simply compute the determinant
    !>  as the square of the [multiplicative trace](@ref pm_matrixTrace) of the Cholesky factor.<br>
    !>
    !>  \param[in]      mat         :   The input `contiguous` matrix of shape `(1:ndim, 1:ndim)` of
    !>                                  <ol>
    !>                                      <li>    type `complex` of kind \CKALL,
    !>                                      <li>    type `real` of kind \RKALL,
    !>                                  </ol>
    !>                                  containing a `subset` of the positive-definite matrix whose `log(sqrt(determinant))` is to be computed.<br>
    !>  \param[in]      subset      :   The input scalar that can be either,
    !>                                  <ol>
    !>                                      <li>    the constant [uppDia](@ref pm_matrixSubset::uppDia) implying that the upper-diagonal triangular block
    !>                                              of `mat` should be used for building the **lower-diagonal** Cholesky factor in the output `chol`.<br>
    !>                                      <li>    the constant [lowDia](@ref pm_matrixSubset::lowDia) implying that the lower-diagonal triangular block
    !>                                              of `mat` should be used for building the **upper-diagonal** Cholesky factor in the output `chol`.<br>
    !>                                  </ol>
    !>                                  This argument is merely a convenience to differentiate the different procedure functionalities within this generic interface.
    !>  \param[out]     detSqrtLog  :   The output scalar of type `real of the same kind as the input `mat` representing the natural logarithm of the square root of its determinant.<br>
    !>  \param[out]     info        :   The output scalar `integer` of default kind \IK that is `0` <b>if and only if</b> the Cholesky factorization succeeds.<br>
    !>                                  Otherwise, it is set to the order of the leading minor of the specified input subset of `mat` that is not positive definite,
    !>                                  indicating the occurrence of an error and that the factorization could not be completed.<br>
    !>                                  **A non-zero value implies failure of the Cholesky factorization.**<br>
    !>                                  A non-zero `info` implied the input `mat` is not positive definite.<br>
    !>  \param[inout]   chol        :   The input/output matrix of the same type, kind, and shape as the input `mat` containing
    !>                                  the Cholesky factorization in its triangular subset dictated by the input arguments `subset` and `operation`.<br>
    !>                                  See the documentation of the input argument `mat` above for more information.<br>
    !>  \param[in]      operation   :   The input scalar that can be either,
    !>                                  <ol>
    !>                                      <li>    the constant [nothing](@ref pm_array::nothing) implying that
    !>                                              the same subset of the output `chol` as the input `subset` must be populated with the Cholesky factorization.<br>
    !>                                      <li>    the constant [transHerm](@ref pm_matrixTrans::transHerm) implying that
    !>                                              the complex transpose of the computed Cholesky factorization must be written to the output `chol`.<br>
    !>                                  </ol>
    !>
    !>  \interface{setMatDetSqrtLog}
    !>  \code{.F90}
    !>
    !>      use pm_matrixDet, only: setMatDetSqrtLog
    !>
    !>      call setMatDetSqrtLog(mat, subset, detSqrtLog, info, chol, operation)
    !>
    !>  \endcode
    !>
    !>  \warning
    !>  The condition `size(mat, 1) == size(mat, 2)` must hold for the corresponding input arguments.<br>
    !>  \vericon
    !>
    !>  \warning
    !>  The input `mat` is assumed to be positive-definite, otherwise the Cholesky Factorization within the procedure will fail,
    !>  in which case, the procedure will return the control with `info /= 0`.<br>
    !>
    !>  \warnpure
    !>
    !>  \see
    !>  [getMatDet](@ref pm_matrixDet::getMatDet)<br>
    !>  [setMatDet](@ref pm_matrixDet::setMatDet)<br>
    !>  [getMatDetSqrt](@ref pm_matrixDet::getMatDetSqrt)<br>
    !>  [setMatDetSqrt](@ref pm_matrixDet::setMatDetSqrt)<br>
    !>  [getMatDetSqrtLog](@ref pm_matrixDet::getMatDetSqrtLog)<br>
    !>  [setMatDetSqrtLog](@ref pm_matrixDet::setMatDetSqrtLog)<br>
    !>  [getMatMulTrace](@ref pm_matrixTrace::getMatMulTrace)<br>
    !>
    !>  \example{setMatDetSqrtLog}
    !>  \include{lineno} example/pm_matrixDet/setMatDetSqrtLog/main.F90
    !>  \compilef{setMatDetSqrtLog}
    !>  \output{setMatDetSqrtLog}
    !>  \include{lineno} example/pm_matrixDet/setMatDetSqrtLog/main.out.F90
    !>
    !>  \test
    !>  [test_pm_matrixDet](@ref test_pm_matrixDet)
    !>
    !>  \final{setMatDetSqrtLog}
    !>
    !>  \author
    !>  \AmirShahmoradi, Apr 21, 2017, 1:43 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

    ! nothing

    interface setMatDetSqrtLog

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_ONO_CK5(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_ONO_CK4(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_ONO_CK3(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_ONO_CK2(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_ONO_CK1(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_ONO_RK5(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_ONO_RK4(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_ONO_RK3(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_ONO_RK2(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_ONO_RK1(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(uppDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_ONO_CK5(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_ONO_CK5
#endif
        use pm_kind, only: CKG => CK5
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_ONO_CK4(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_ONO_CK4
#endif
        use pm_kind, only: CKG => CK4
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_ONO_CK3(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_ONO_CK3
#endif
        use pm_kind, only: CKG => CK3
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_ONO_CK2(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_ONO_CK2
#endif
        use pm_kind, only: CKG => CK2
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_ONO_CK1(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_ONO_CK1
#endif
        use pm_kind, only: CKG => CK1
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_ONO_RK5(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_ONO_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_ONO_RK4(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_ONO_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_ONO_RK3(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_ONO_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_ONO_RK2(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_ONO_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_ONO_RK1(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_ONO_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(lowDia_type)   , intent(in)                    :: subset
        type(nothing_type)  , intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

    ! transHerm

    interface setMatDetSqrtLog

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_OTH_CK5(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_OTH_CK5
#endif
        use pm_kind, only: CKG => CK5
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_OTH_CK4(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_OTH_CK4
#endif
        use pm_kind, only: CKG => CK4
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_OTH_CK3(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_OTH_CK3
#endif
        use pm_kind, only: CKG => CK3
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_OTH_CK2(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_OTH_CK2
#endif
        use pm_kind, only: CKG => CK2
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_OTH_CK1(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_OTH_CK1
#endif
        use pm_kind, only: CKG => CK1
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_OTH_RK5(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_OTH_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_OTH_RK4(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_OTH_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_OTH_RK3(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_OTH_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_OTH_RK2(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_OTH_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatDetSqrtLog_UXD_OTH_RK1(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_UXD_OTH_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(uppDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if CK5_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_OTH_CK5(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_OTH_CK5
#endif
        use pm_kind, only: CKG => CK5
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK4_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_OTH_CK4(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_OTH_CK4
#endif
        use pm_kind, only: CKG => CK4
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK3_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_OTH_CK3(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_OTH_CK3
#endif
        use pm_kind, only: CKG => CK3
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK2_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_OTH_CK2(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_OTH_CK2
#endif
        use pm_kind, only: CKG => CK2
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if CK1_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_OTH_CK1(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_OTH_CK1
#endif
        use pm_kind, only: CKG => CK1
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        complex(CKG)        , intent(in)    , contiguous    :: mat(:,:)
        complex(CKG)        , intent(inout) , contiguous    :: chol(:,:)
        real(CKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if RK5_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_OTH_RK5(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_OTH_RK5
#endif
        use pm_kind, only: RKG => RK5
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK4_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_OTH_RK4(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_OTH_RK4
#endif
        use pm_kind, only: RKG => RK4
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK3_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_OTH_RK3(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_OTH_RK3
#endif
        use pm_kind, only: RKG => RK3
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK2_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_OTH_RK2(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_OTH_RK2
#endif
        use pm_kind, only: RKG => RK2
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

#if RK1_ENABLED
    PURE module subroutine setMatDetSqrtLog_XLD_OTH_RK1(mat, subset, detSqrtLog, info, chol, operation)
#if __INTEL_COMPILER && DLL_ENABLED && (_WIN32 || _WIN64)
        !DEC$ ATTRIBUTES DLLEXPORT :: setMatDetSqrtLog_XLD_OTH_RK1
#endif
        use pm_kind, only: RKG => RK1
        type(lowDia_type)   , intent(in)                    :: subset
        type(transHerm_type), intent(in)                    :: operation
        real(RKG)           , intent(in)    , contiguous    :: mat(:,:)
        real(RKG)           , intent(inout) , contiguous    :: chol(:,:)
        real(RKG)           , intent(out)                   :: detSqrtLog
        integer(IK)         , intent(out)                   :: info
    end subroutine
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_matrixDet ! LCOV_EXCL_LINE