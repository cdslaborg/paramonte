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
!>  This module contains a set of generic interfaces to the [BLAS routines](https://www.netlib.org/blas/).<br>
!>
!>  \details
!>  The BLAS generic interfaces of this module facilitate runtime dispatch of matrix algebra
!>  to the appropriate BLAS routines of an optimized BLAS library linked to the ParaMonte library.<br>
!>  The runtime dispatch occurs only if there is a corresponding routine (single or double precision) in the optimized BLAS library.<br>
!>
!>  \see
!>  [pm_lapack](@ref pm_lapack)<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_blas

    use pm_kind, only: SK, IK, LK, RKS, RKD

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_blas"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Interchange the elements of vectors `x` and `y`.<br>
    !>
    !>  \details
    !>  This generic interface offers a compile-time resolution to the
    !>  BLAS `sswap`, `dswap`, `cswap`, `zswap`, swapping routines.<br>
    !>  See the documentation of reference BLAS library for the definition of the input arguments.<br>
    !>
    !>  \see
    !>  [pm_matrixMulAdd](@ref pm_matrixMulAdd)<br>
    !>
    !>  \final{blasGEMV}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface blasSWAP
        pure subroutine zswap(n, x, incx, y, incy)
            use pm_kind , only: IK, TKC => RKD
            integer(IK) , intent(in) :: n, incx, incy
            complex(TKC), intent(inout) :: x(*), y(*)
        end subroutine zswap
        pure subroutine cswap(n, x, incx, y, incy)
            use pm_kind , only: IK, TKC => RKS
            integer(IK) , intent(in) :: n, incx, incy
            complex(TKC), intent(inout) :: x(*), y(*)
        end subroutine cswap
        pure subroutine dswap(n, x, incx, y, incy)
            use pm_kind , only: IK, TKC => RKD
            integer(IK) , intent(in) :: n, incx, incy
            real(TKC)   , intent(inout) :: x(*), y(*)
        end subroutine dswap
        pure subroutine sswap(n, x, incx, y, incy)
            use pm_kind , only: IK, TKC => RKS
            integer(IK) , intent(in) :: n, incx, incy
            real(TKC)   , intent(inout) :: x(*), y(*)
        end subroutine sswap
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the result of the multiplication of a
    !>  General matrix `matA` in [Rectangular Default packing](@ref pm_matrixPack) format
    !>  with a vector/column-like matrix `matB`, added to a third vector/column-like matrix `matC`.<br>
    !>
    !>  \details
    !>  This generic interface offers a compile-time resolution to the
    !>  BLAS `sgemv`, `dgemv`, `cgemv`, `zgemv`, General matrix-vector multiplication routines.<br>
    !>  See the documentation of reference BLAS library for the definition of the input arguments.<br>
    !>
    !>  \see
    !>  [pm_matrixMulAdd](@ref pm_matrixMulAdd)<br>
    !>
    !>  \final{blasGEMV}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface blasGEMV
        pure subroutine zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            complex(TKC), intent(inout) :: y
            complex(TKC), intent(in)    :: alpha, beta, a, x
            integer(IK) , intent(in)    :: m, n, lda, incx, incy
            character   , intent(in)    :: trans
        end subroutine
        pure subroutine cgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            complex(TKC), intent(inout) :: y
            complex(TKC), intent(in)    :: alpha, beta, a, x
            integer(IK) , intent(in)    :: m, n, lda, incx, incy
            character   , intent(in)    :: trans
        end subroutine
        pure subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            real(TKC)   , intent(inout) :: y
            real(TKC)   , intent(in)    :: alpha, beta, a, x
            integer(IK) , intent(in)    :: m, n, lda, incx, incy
            character   , intent(in)    :: trans
        end subroutine
        pure subroutine sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            real(TKC)   , intent(inout) :: y
            real(TKC)   , intent(in)    :: alpha, beta, a, x
            integer(IK) , intent(in)    :: m, n, lda, incx, incy
            character   , intent(in)    :: trans
        end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the result of the multiplication of a
    !>  Hermitian matrix `matA` in [Linear Full Packed](@ref pm_matrixPack) format
    !>  with a vector/column-like matrix `matB`, added to a third vector/column-like matrix `matC`.<br>
    !>
    !>  \details
    !>  This generic interface offers a compile-time resolution to the
    !>  BLAS `chpmv`, `zhpmv` Hermitian matrix multiplication routines.<br>
    !>  See the documentation of reference BLAS library for the definition of the input arguments.<br>
    !>
    !>  \see
    !>  [pm_matrixMulAdd](@ref pm_matrixMulAdd)<br>
    !>
    !>  \final{blasHPMV}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface blasHPMV
        pure subroutine zhpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            complex(TKC), intent(inout) :: y
            complex(TKC), intent(in)    :: alpha, beta, ap, x
            integer(IK) , intent(in)    :: n, incx, incy
            character   , intent(in)    :: uplo
        end subroutine
        pure subroutine chpmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            complex(TKC), intent(inout) :: y
            complex(TKC), intent(in)    :: alpha, beta, ap, x
            integer(IK) , intent(in)    :: n, incx, incy
            character   , intent(in)    :: uplo
        end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the result of the multiplication of a
    !>  Symmetric matrix `matA` in [Linear Full Packed](@ref pm_matrixPack) format
    !>  with a vector/column-like matrix `matB`, added to a third vector/column-like matrix `matC`.<br>
    !>
    !>  \details
    !>  This generic interface offers a compile-time resolution to the
    !>  BLAS `sspmv`, `dspmv` Symmetric matrix multiplication routines.<br>
    !>  See the documentation of reference BLAS library for the definition of the input arguments.<br>
    !>
    !>  \see
    !>  [pm_matrixMulAdd](@ref pm_matrixMulAdd)<br>
    !>
    !>  \final{blasSPMV}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface blasSPMV
        pure subroutine dspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            real(TKC)   , intent(inout) :: y
            real(TKC)   , intent(in)    :: alpha, beta, ap, x
            integer(IK) , intent(in)    :: n, incx, incy
            character   , intent(in)    :: uplo
        end subroutine
        pure subroutine sspmv(uplo, n, alpha, ap, x, incx, beta, y, incy)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            real(TKC)   , intent(inout) :: y
            real(TKC)   , intent(in)    :: alpha, beta, ap, x
            integer(IK) , intent(in)    :: n, incx, incy
            character   , intent(in)    :: uplo
        end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the result of the multiplication of a
    !>  Hermitian matrix `matA` in [Rectangular Default](@ref pm_matrixPack) packing format
    !>  with a vector/column-like matrix `matB`, added to a third vector/column-like matrix `matC`.<br>
    !>
    !>  \details
    !>  This generic interface offers a compile-time resolution to the
    !>  BLAS `chemv`, `zhemv` Hermitian matrix multiplication routines.<br>
    !>  See the documentation of reference BLAS library for the definition of the input arguments.<br>
    !>
    !>  \see
    !>  [pm_matrixMulAdd](@ref pm_matrixMulAdd)<br>
    !>
    !>  \final{blasHEMV}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface blasHEMV
        pure subroutine zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            complex(TKC), intent(inout) :: y
            complex(TKC), intent(in)    :: alpha, beta, a, x
            integer(IK) , intent(in)    :: n, lda, incx, incy
            character   , intent(in)    :: uplo
        end subroutine
        pure subroutine chemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            complex(TKC), intent(inout) :: y
            complex(TKC), intent(in)    :: alpha, beta, a, x
            integer(IK) , intent(in)    :: n, lda, incx, incy
            character   , intent(in)    :: uplo
        end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the result of the multiplication of a
    !>  Symmetric matrix `matA` in [Rectangular Default](@ref pm_matrixPack) packing format
    !>  with a vector/column-like matrix `matB`, added to a third vector/column-like matrix `matC`.<br>
    !>
    !>  \details
    !>  This generic interface offers a compile-time resolution to the
    !>  BLAS `ssymv`, `dsymv` Symmetric matrix multiplication routines.<br>
    !>  See the documentation of reference BLAS library for the definition of the input arguments.<br>
    !>
    !>  \see
    !>  [pm_matrixMulAdd](@ref pm_matrixMulAdd)<br>
    !>
    !>  \final{blasSYMV}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface blasSYMV
        pure subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            real(TKC)   , intent(inout) :: y
            real(TKC)   , intent(in)    :: alpha, beta, a, x
            integer(IK) , intent(in)    :: n, lda, incx, incy
            character   , intent(in)    :: uplo
        end subroutine
        pure subroutine ssymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            real(TKC)   , intent(inout) :: y
            real(TKC)   , intent(in)    :: alpha, beta, a, x
            integer(IK) , intent(in)    :: n, lda, incx, incy
            character   , intent(in)    :: uplo
        end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the result of the multiplication of two input
    !>  Symmetric-General/General-Symmetric matrices `matA` and `matB` added to a third General matrix `matC`.
    !>
    !>  \details
    !>  This generic interface offers a compile-time resolution to the
    !>  BLAS `ssymm`, `dsymm`, `csymm`, `zsymm` Symmetric-General/General-Symmetric matrix multiplication routines.<br>
    !>  See the documentation of reference BLAS library for the definition of the input arguments.<br>
    !>
    !>  \see
    !>  [pm_matrixMulAdd](@ref pm_matrixMulAdd)<br>
    !>
    !>  \final{blasSYMM}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface blasSYMM
        pure subroutine zsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            complex(TKC), intent(inout) :: c
            complex(TKC), intent(in)    :: alpha, beta, a, b
            integer(IK) , intent(in)    :: m, n, lda, ldb, ldc
            character   , intent(in)    :: side, uplo
        end subroutine
        pure subroutine csymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            complex(TKC), intent(inout) :: c
            complex(TKC), intent(in)    :: alpha, beta, a, b
            integer(IK) , intent(in)    :: m, n, lda, ldb, ldc
            character   , intent(in)    :: side, uplo
        end subroutine
        pure subroutine dsymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            real(TKC)   , intent(inout) :: c
            real(TKC)   , intent(in)    :: alpha, beta, a, b
            integer(IK) , intent(in)    :: m, n, lda, ldb, ldc
            character   , intent(in)    :: side, uplo
        end subroutine
        pure subroutine ssymm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            real(TKC)   , intent(inout) :: c
            real(TKC)   , intent(in)    :: alpha, beta, a, b
            integer(IK) , intent(in)    :: m, n, lda, ldb, ldc
            character   , intent(in)    :: side, uplo
        end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the result of the multiplication of two input
    !>  Hermitian-General/General-Hermitian matrices `matA` and `matB` added to a third General matrix `matC`.
    !>
    !>  \details
    !>  This generic interface offers a compile-time resolution to the
    !>  BLAS `chemm`, `zhemm` Hermitian-General/General-Hermitian matrix multiplication routines.<br>
    !>  See the documentation of reference BLAS library for the definition of the input arguments.<br>
    !>
    !>  \see
    !>  [pm_matrixMulAdd](@ref pm_matrixMulAdd)<br>
    !>
    !>  \final{blasHEMM}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface blasHEMM
        pure subroutine zhemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            complex(TKC), intent(inout) :: c
            complex(TKC), intent(in)    :: alpha, beta, a, b
            integer(IK) , intent(in)    :: m, n, lda, ldb, ldc
            character   , intent(in)    :: side, uplo
        end subroutine zhemm
        pure subroutine chemm(side, uplo, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            complex(TKC), intent(inout) :: c
            complex(TKC), intent(in)    :: alpha, beta, a, b
            integer(IK) , intent(in)    :: m, n, lda, ldb, ldc
            character   , intent(in)    :: side, uplo
        end subroutine chemm
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the result of the multiplication of two input General matrices `matA` and `matB` added to a third General matrix `matC`.
    !>
    !>  \details
    !>  This generic interface offers a compile-time resolution to the
    !>  BLAS `sgemm`, `dgemm`, `cgemm`, `zgemm` General matrix multiplication routines.<br>
    !>  See the documentation of reference BLAS library for the definition of the input arguments.<br>
    !>
    !>  \see
    !>  [pm_matrixMulAdd](@ref pm_matrixMulAdd)<br>
    !>
    !>  \final{blasGEMM}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface blasGEMM
    pure subroutine zgemm(transa, transb, l, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
        use pm_kind , only: IK, TKC => RKD
        implicit none
        complex(TKC), intent(inout) :: c
        complex(TKC), intent(in)    :: alpha, beta, a, b
        integer(IK) , intent(in)    :: l, m, n, lda, ldb, ldc
        character   , intent(in)    :: transa, transb
    end subroutine
    pure subroutine cgemm(transa, transb, l, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
        use pm_kind , only: IK, TKC => RKS
        implicit none
        complex(TKC), intent(inout) :: c
        complex(TKC), intent(in)    :: alpha, beta, a, b
        integer(IK) , intent(in)    :: l, m, n, lda, ldb, ldc
        character   , intent(in)    :: transa, transb
    end subroutine
    pure subroutine dgemm(transa, transb, l, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
        use pm_kind , only: IK, TKC => RKD
        implicit none
        real(TKC)   , intent(inout) :: c
        real(TKC)   , intent(in)    :: alpha, beta, a, b
        integer(IK) , intent(in)    :: l, m, n, lda, ldb, ldc
        character   , intent(in)    :: transa, transb
    end subroutine
    pure subroutine sgemm(transa, transb, l, m, n, alpha, a, lda, b, ldb, beta, c, ldc)
        use pm_kind , only: IK, TKC => RKS
        implicit none
        real(TKC)   , intent(inout) :: c
        real(TKC)   , intent(in)    :: alpha, beta, a, b
        integer(IK) , intent(in)    :: l, m, n, lda, ldb, ldc
        character   , intent(in)    :: transa, transb
    end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the following matrix-vector products, using the vector `x` and triangular matrix `A` or its transpose:<br>
    !>  \f{equation}{
    !>      x\leftarrow A x
    !>      x\leftarrow A^T x
    !>      x\leftarrow A^H x
    !>  \f}
    !>
    !>  \details
    !>  This generic interface offers a compile-time resolution to the
    !>  BLAS `strmv`, `dtrmv`, `ctrmv`, `ztrmv` matrix multiplication routines.<br>
    !>  See the documentation of reference BLAS library for the definition of the input arguments.<br>
    !>
    !>  \see
    !>  [pm_matrixMulTri](@ref pm_matrixMulTri)<br>
    !>
    !>  \final{blasTRMV}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface blasTRMV
        pure subroutine ztrmv(uplo, trans, diag, n, a, lda, x, incx)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            character   , intent(in)    :: uplo, trans, diag
            integer(IK) , intent(in)    :: n, lda, incx
            complex(TKC), intent(inout) :: x
            complex(TKC), intent(in)    :: a
        end subroutine
        pure subroutine ctrmv(uplo, trans, diag, n, a, lda, x, incx)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            character   , intent(in)    :: uplo, trans, diag
            integer(IK) , intent(in)    :: n, lda, incx
            complex(TKC), intent(inout) :: x
            complex(TKC), intent(in)    :: a
        end subroutine
        pure subroutine dtrmv(uplo, trans, diag, n, a, lda, x, incx)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            character   , intent(in)    :: uplo, trans, diag
            integer(IK) , intent(in)    :: n, lda, incx
            real(TKC)   , intent(inout) :: x
            real(TKC)   , intent(in)    :: a
        end subroutine
        pure subroutine strmv(uplo, trans, diag, n, a, lda, x, incx)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            character   , intent(in)    :: uplo, trans, diag
            integer(IK) , intent(in)    :: n, lda, incx
            real(TKC)   , intent(inout) :: x
            real(TKC)   , intent(in)    :: a
        end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the following matrix-vector products, using the vector `x` and triangular matrix `A` or its transpose:<br>
    !>  \f{equation}{
    !>      x\leftarrow A-1 x
    !>      x\leftarrow A^-T x
    !>      x\leftarrow A^-H x
    !>  \f}
    !>
    !>  \details
    !>  This generic interface offers a compile-time resolution to the
    !>  BLAS `strsv`, `dtrsv`, `ctrsv`, `ztrsv` matrix multiplication routines.<br>
    !>  See the documentation of reference BLAS library for the definition of the input arguments.<br>
    !>
    !>  \see
    !>  [pm_matrixMulTri](@ref pm_matrixMulTri)<br>
    !>
    !>  \final{blasTRSV}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface blasTRSV
        pure subroutine ztrsv(uplo, trans, diag, n, a, lda, x, incx)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            character   , intent(in)    :: uplo, trans, diag
            integer(IK) , intent(in)    :: n, lda, incx
            complex(TKC), intent(inout) :: x
            complex(TKC), intent(in)    :: a
        end subroutine
        pure subroutine ctrsv(uplo, trans, diag, n, a, lda, x, incx)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            character   , intent(in)    :: uplo, trans, diag
            integer(IK) , intent(in)    :: n, lda, incx
            complex(TKC), intent(inout) :: x
            complex(TKC), intent(in)    :: a
        end subroutine
        pure subroutine dtrsv(uplo, trans, diag, n, a, lda, x, incx)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            character   , intent(in)    :: uplo, trans, diag
            integer(IK) , intent(in)    :: n, lda, incx
            real(TKC)   , intent(inout) :: x
            real(TKC)   , intent(in)    :: a
        end subroutine
        pure subroutine strsv(uplo, trans, diag, n, a, lda, x, incx)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            character   , intent(in)    :: uplo, trans, diag
            integer(IK) , intent(in)    :: n, lda, incx
            real(TKC)   , intent(inout) :: x
            real(TKC)   , intent(in)    :: a
        end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the result of the triangular/general matrix multiplication, using scalar `alpha`,
    !>  rectangular general matrix `B`, and triangular matrix `A` or its Symmetric/Hermitian transpose.<br>
    !>
    !>  \details
    !>  This generic interface offers a compile-time resolution to the
    !>  BLAS `strmm`, `dtrmm`, `ctrmm`, `ztrmm` matrix multiplication routines.<br>
    !>  See the documentation of reference BLAS library for the definition of the input arguments.<br>
    !>
    !>  \see
    !>  [pm_matrixMulTri](@ref pm_matrixMulTri)<br>
    !>
    !>  \final{blasTRMM}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface blasTRMM
        pure subroutine ztrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            integer(IK) , intent(in)    :: m, n, lda, ldb
            complex(TKC), intent(inout) :: b
            complex(TKC), intent(in)    :: alpha, a
            character   , intent(in)    :: side, uplo, transa, diag
        end subroutine
        pure subroutine ctrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            integer(IK) , intent(in)    :: m, n, lda, ldb
            complex(TKC), intent(inout) :: b
            complex(TKC), intent(in)    :: alpha, a
            character   , intent(in)    :: side, uplo, transa, diag
        end subroutine
        pure subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            integer(IK) , intent(in)    :: m, n, lda, ldb
            real(TKC)   , intent(inout) :: b
            real(TKC)   , intent(in)    :: alpha, a
            character   , intent(in)    :: side, uplo, transa, diag
        end subroutine
        pure subroutine strmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            integer(IK) , intent(in)    :: m, n, lda, ldb
            real(TKC)   , intent(inout) :: b
            real(TKC)   , intent(in)    :: alpha, a
            character   , intent(in)    :: side, uplo, transa, diag
        end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the solution to a triangular system of equations with multiple right-hand sides, using scalar `alpha`,
    !>  rectangular general matrix `B`, and triangular matrix `A` or its Symmetric/Hermitian transpose.<br>
    !>
    !>  \details
    !>  This generic interface offers a compile-time resolution to the
    !>  BLAS `strsm`, `dtrsm`, `ctrsm`, `ztrsm` matrix multiplication routines.<br>
    !>  See the documentation of reference BLAS library for the definition of the input arguments.<br>
    !>
    !>  \see
    !>  [pm_matrixMulTri](@ref pm_matrixMulTri)<br>
    !>
    !>  \final{blasTRSM}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface blasTRSM
        pure subroutine ztrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            integer(IK) , intent(in)    :: m, n, lda, ldb
            complex(TKC), intent(inout) :: b
            complex(TKC), intent(in)    :: alpha, a
            character   , intent(in)    :: side, uplo, transa, diag
        end subroutine
        pure subroutine ctrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            integer(IK) , intent(in)    :: m, n, lda, ldb
            complex(TKC), intent(inout) :: b
            complex(TKC), intent(in)    :: alpha, a
            character   , intent(in)    :: side, uplo, transa, diag
        end subroutine
        pure subroutine dtrsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            integer(IK) , intent(in)    :: m, n, lda, ldb
            real(TKC)   , intent(inout) :: b
            real(TKC)   , intent(in)    :: alpha, a
            character   , intent(in)    :: side, uplo, transa, diag
        end subroutine
        pure subroutine strsm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            integer(IK) , intent(in)    :: m, n, lda, ldb
            real(TKC)   , intent(inout) :: b
            real(TKC)   , intent(in)    :: alpha, a
            character   , intent(in)    :: side, uplo, transa, diag
        end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Rank-K Update of a Real or Complex Symmetric or a Complex Hermitian matrix.<br>
    !>  \f{equation}{
    !>      C = \alpha A A^T + \beta C,
    !>      C = \alpha A^T A + \beta C,
    !>  \f}
    !>
    !>  \details
    !>  This generic interface offers a compile-time resolution to the
    !>  BLAS `ssyrk`, `dsyrk`, `csyrk`, `zsyrk` matrix update routines.<br>
    !>  See the documentation of reference BLAS library for the definition of the input arguments.<br>
    !>
    !>  \see
    !>  [pm_matrixUpdate](@ref pm_matrixUpdate)<br>
    !>
    !>  \final{blasSYRK}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface blasSYRK
        pure subroutine zsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            character   , intent(in)    :: uplo, trans
            integer(IK) , intent(in)    :: n, k, lda, ldc
            complex(TKC), intent(in)    :: a, alpha, beta
            complex(TKC), intent(inout) :: c
        end subroutine
        pure subroutine csyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            character   , intent(in)    :: uplo, trans
            integer(IK) , intent(in)    :: n, k, lda, ldc
            complex(TKC), intent(in)    :: a, alpha, beta
            complex(TKC), intent(inout) :: c
        end subroutine
        pure subroutine dsyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            character   , intent(in)    :: uplo, trans
            integer(IK) , intent(in)    :: n, k, lda, ldc
            real(TKC)   , intent(in)    :: a, alpha, beta
            real(TKC)   , intent(inout) :: c
        end subroutine
        pure subroutine ssyrk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            character   , intent(in)    :: uplo, trans
            integer(IK) , intent(in)    :: n, k, lda, ldc
            real(TKC)   , intent(in)    :: a, alpha, beta
            real(TKC)   , intent(inout) :: c
        end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  Return the Rank-K Update of a Complex Hermitian matrix.<br>
    !>  \f{equation}{
    !>      C = \alpha A A^H + \beta C,
    !>      C = \alpha A^H A + \beta C,
    !>  \f}
    !>
    !>  \details
    !>  This generic interface offers a compile-time resolution to the
    !>  BLAS `cherk`, `zherk` matrix update routines.<br>
    !>  See the documentation of reference BLAS library for the definition of the input arguments.<br>
    !>
    !>  \see
    !>  [pm_matrixUpdate](@ref pm_matrixUpdate)<br>
    !>
    !>  \final{blasHERK}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface blasHERK
        pure subroutine zherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
            use pm_kind , only: IK, TKC => RKD
            implicit none
            character   , intent(in)    :: uplo, trans
            integer(IK) , intent(in)    :: n, k, lda, ldc
            real(TKC)   , intent(in)    :: alpha, beta
            complex(TKC), intent(in)    :: a
            complex(TKC), intent(inout) :: c
        end subroutine
        pure subroutine cherk(uplo, trans, n, k, alpha, a, lda, beta, c, ldc)
            use pm_kind , only: IK, TKC => RKS
            implicit none
            character   , intent(in)    :: uplo, trans
            integer(IK) , intent(in)    :: n, k, lda, ldc
            real(TKC)   , intent(in)    :: alpha, beta
            complex(TKC), intent(in)    :: a
            complex(TKC), intent(inout) :: c
        end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_blas ! LCOV_EXCL_LINE