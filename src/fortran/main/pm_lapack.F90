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
!>  This module contains a set of generic interfaces to the [LAPACK routines](https://www.netlib.org/lapack/).<br>
!>
!>  \details
!>  The LAPACK generic interfaces of this module facilitate runtime dispatch of matrix algebra
!>  to the appropriate LAPACK routines of an optimized LAPACK library linked to the ParaMonte library.<br>
!>  The runtime dispatch occurs only if there is a corresponding routine (single or double precision) in the optimized LAPACK library.<br>
!>
!>  \see
!>  [pm_blas](@ref pm_blas)<br>
!>
!>  \finmain
!>
!>  \author
!>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_lapack

    use pm_kind, only: SK, IK, LK, RKS, RKD

    implicit none

    character(*, SK), parameter :: MODULE_NAME = "@pm_lapack"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !>  \brief
    !>  computes an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges.<br>
    !>
    !>  \details
    !>  This generic interface offers a compile-time resolution to the
    !>  LAPACK `sgetrf`, `dgetrf`, `cgetrf`, `zgetrf`, LAPACK routines.<br>
    !>  See the documentation of reference LAPACK library for the definition of the input arguments.<br>
    !>
    !>  \see
    !>  [pm_matrixLUP](@ref pm_matrixLUP)<br>
    !>
    !>  \finmain{lapackGETRF}
    !>
    !>  \author
    !>  \AmirShahmoradi, Friday 1:54 AM, April 21, 2017, Institute for Computational Engineering and Sciences (ICES), The University of Texas, Austin, TX
    interface lapackGETRF
        pure subroutine zgetrf(m, n, a, lda, ipiv, info)
            use pm_kind , only: IK, TKC => CKD
            integer(IK) , intent(in) :: m, n, lda
            integer(IK) , intent(out) :: info, ipiv(*)
            complex(TKC), intent(inout) :: a
        end subroutine
        pure subroutine cgetrf(m, n, a, lda, ipiv, info)
            use pm_kind , only: IK, TKC => CKS
            integer(IK) , intent(in) :: m, n, lda
            integer(IK) , intent(out) :: info, ipiv(*)
            complex(TKC), intent(inout) :: a
        end subroutine
        pure subroutine dgetrf(m, n, a, lda, ipiv, info)
            use pm_kind , only: IK, TKC => RKD
            integer(IK) , intent(in) :: m, n, lda
            integer(IK) , intent(out) :: info, ipiv(*)
            real(TKC)   , intent(inout) :: a
        end subroutine
        pure subroutine sgetrf(m, n, a, lda, ipiv, info)
            use pm_kind , only: IK, TKC => RKS
            integer(IK) , intent(in) :: m, n, lda
            integer(IK) , intent(out) :: info, ipiv(*)
            real(TKC)   , intent(inout) :: a
        end subroutine
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_lapack ! LCOV_EXCL_LINE