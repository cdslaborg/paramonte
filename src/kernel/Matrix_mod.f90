!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!   MIT License
!!!!
!!!!   ParaMonte: plain powerful parallel Monte Carlo library.
!!!!
!!!!   Copyright (C) 2012-present, The Computational Data Science Lab
!!!!
!!!!   This file is part of the ParaMonte library.
!!!!
!!!!   Permission is hereby granted, free of charge, to any person obtaining a
!!!!   copy of this software and associated documentation files (the "Software"),
!!!!   to deal in the Software without restriction, including without limitation
!!!!   the rights to use, copy, modify, merge, publish, distribute, sublicense,
!!!!   and/or sell copies of the Software, and to permit persons to whom the
!!!!   Software is furnished to do so, subject to the following conditions:
!!!!
!!!!   The above copyright notice and this permission notice shall be
!!!!   included in all copies or substantial portions of the Software.
!!!!
!!!!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!!!!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!!!!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!!!!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
!!!!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
!!!!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
!!!!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!!!!
!!!!   ACKNOWLEDGMENT
!!!!
!!!!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!!!!   As per the ParaMonte library license agreement terms, if you use any parts of
!!!!   this library for any purposes, kindly acknowledge the use of ParaMonte in your
!!!!   work (education/research/industry/development/...) by citing the ParaMonte
!!!!   library as described on this page:
!!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This module contains mathematical procedures.
!>  \author Amir Shahmoradi

module Matrix_mod

    implicit none

    character(*), parameter :: MODULE_NAME = "@Matrix_mod"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the Cholesky factorization of the input positive-definite matrix.
    !>
    !> \param[in]       nd          :   The size of the input square matrix - `nd` by `nd`.
    !> \param[in,out]   PosDefMat   :   The input square matrix.
    !> \param[out]      Diagonal    :   The Diagonal elements of the Cholesky factorization.
    !>
    !> \remark
    !> If `nd = 1`, `PosDefMat` will not be touched, only `sqrt(PosDefMat)` will be output to `Diagonal`.
    !>
    !> \warning
    !> If Cholesky factorization fails, `Diagonal(1)` will be set to `-1` to indicate error on return.
    !>
    !> \details
    !> Returns in the lower triangle of `PosDefMat`, the Cholesky factorization L of \f$\texttt{PosDefMat} = L.L^T\f$.
    !> On input, the upper upper triangle of `PosDefMat` should be given, which remains intact on output.
    pure subroutine getCholeskyFactor(nd,PosDefMat,Diagonal)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getCholeskyFactor
#endif
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK), intent(in)    :: nd
        real(RK)   , intent(inout) :: PosDefMat(nd,nd) ! Upper triangle + diagonal is input matrix, lower is output.
        real(RK)   , intent(out)   :: Diagonal(nd)
        real(RK)                   :: summ
        integer(IK)                :: i
        do i=1,nd
            summ = PosDefMat(i,i) - dot_product(PosDefMat(i,1:i-1),PosDefMat(i,1:i-1))
            if (summ <= 0._RK) then
                Diagonal(1) = -1._RK
                return
            end if
            Diagonal(i) = sqrt(summ)
            PosDefMat(i+1:nd,i) = ( PosDefMat(i,i+1:nd) - matmul(PosDefMat(i+1:nd,1:i-1),PosDefMat(i,1:i-1)) ) / Diagonal(i)
        end do
    end subroutine getCholeskyFactor

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Solve the linear equation system of the form: \f$ \texttt{PosDefMat} \times \texttt{InputSolution} = \texttt{Intercept} \f$
    !>
    !> \param[in]       nd              :   The size of the input square matrix - `nd` by `nd`.
    !> \param[in]       PosDefMat       :   The input square matrix.
    !> \param[in]       Diagonal        :   The Diagonal elements of the Cholesky factorization.
    !> \param[in]       Intercept       :   The intercept.
    !> \param[in,out]   InputSolution   :   The input right-hand-side which becomes the solution on return.
    !>
    !> \remark
    !> `PosDefMat` and `Diagonal` are the output of subroutine [getCholeskyFactor](@ref getcholeskyfactor)
    !> (i.g., only the lower triangle of `PosDefMat` is used).
    subroutine solveLinearPosDefSystem(nd,PosDefMat,Diagonal,Intercept,InputSolution)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: solveLinearPosDefSystem
#endif
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK), intent(in)    :: nd
        real(RK)   , intent(in)    :: PosDefMat(nd,nd),Diagonal(nd),Intercept(nd)
        real(RK)   , intent(inout) :: InputSolution(nd)
        integer(IK)                :: i
        do i=1,nd
            InputSolution(i) = ( Intercept(i) - dot_product(PosDefMat(i,1:i-1),InputSolution(1:i-1)) ) / Diagonal(i)
        end do
        do i = nd,1,-1
            InputSolution(i) = ( InputSolution(i) - dot_product(PosDefMat(i+1:nd,i),InputSolution(i+1:nd)) ) / Diagonal(i)
        end do
    end subroutine solveLinearPosDefSystem

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the inverse matrix of a symmetric-positive-definite input matrix, which is given in the upper triangle of `MatInvMat`.
    !> On output `MatInvMat` is completely overwritten by in the inverse of the matrix. Also, return the square root of determinant
    !> of the inverse matrix.
    !>
    !> \param[in]       nd                  :   The size of the input square matrix - `nd` by `nd`.
    !> \param[in,out]   MatInvMat           :   The input square matrix. On input, the upper half must be covariance matrix.
    !>                                          On output, it is completely overwritten by the inverse matrix.
    !> \param[out]      sqrtDetInvPosDefMat :   The square root of the determinant of the inverse matrix.
    !>
    !> \warning
    !> Do not call this routine when `nd = 1`. There is no need to call this routine when `nd = 1`.
    !>
    !> \warning
    !> If the input matrix is not positive-definite, the output `sqrtDetInvPosDefMat`
    !> will be set to `-1` on return to indicate the occurrence of an error.
    !>
    !> \author
    !> Amir Shahmoradi, Apr 21, 2017, 1:54 AM, ICES, UT Austin
    pure subroutine getInvPosDefMatSqrtDet(nd,MatInvMat,sqrtDetInvPosDefMat)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getInvPosDefMatSqrtDet
#endif
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK), intent(in)    :: nd
        real(RK)   , intent(inout) :: MatInvMat(nd,nd)           ! input: upper half is covariance matrix, output: inverse matrix
        real(RK)   , intent(out)   :: sqrtDetInvPosDefMat        !
        real(RK)                   :: CholeskyLower(nd,nd)       ! Cholesky factor
        real(RK)                   :: Diagonal(nd)               ! diagonal terms of the inverse matrix
        real(RK)                   :: summ
        integer(IK)                :: i,j,k
        if (nd==1_IK) then
            MatInvMat = 1._RK / MatInvMat
            sqrtDetInvPosDefMat = MatInvMat(1,1)
            return
        end if
        do j=1,nd
            do i=1,j
                CholeskyLower(i,j) = MatInvMat(i,j)
            end do
        end do
        call getCholeskyFactor(nd,CholeskyLower,Diagonal)
        if (Diagonal(1)<0._RK) then
            sqrtDetInvPosDefMat = -1._RK
            return
        end if
        sqrtDetInvPosDefMat = 1._RK / product(Diagonal)
        do i = 1,nd
            CholeskyLower(i,i) = 1._RK / Diagonal(i)
            do j = i+1,nd
                summ = 0._RK
                do k = i,j-1
                    summ = summ - CholeskyLower(j,k) * CholeskyLower(k,i)
                end do
                CholeskyLower(j,i) = summ / Diagonal(j)
            end do
        end do
        do i = 1,nd
            do j = i,nd
                MatInvMat(j,i) = dot_product(CholeskyLower(j:nd,j),CholeskyLower(j:nd,i))
            end do
            MatInvMat(i,i:nd) = MatInvMat(i:nd,i)
        end do
  end subroutine getInvPosDefMatSqrtDet

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the inverse matrix of a symmetric-positive-definite matrix, whose Cholesky Lower triangle is given in the lower part
    !> of `CholeskyLower` and and its diagonal elements in `Diagonal`.
    !>
    !> \param[in]       nd              :   The size of the input square matrix - `nd` by `nd`.
    !> \param[in]       CholeskyLower   :   The Cholesky factorization of the matrix.
    !> \param[in]       Diagonal        :   The diagonal elements of the Cholesky factorization of the matrix.
    !>
    !> \return
    !> `InvMatFromCholFac` : The full inverse matrix.
    !>
    !> \warning
    !> Do not call this routine when `nd = 1`. For `nd = 1`: `InvMatFromCholFac = 1._RK / Diagonal(1)^2`.
    !>
    !> \author
    !> Amir Shahmoradi, Apr 21, 2017, 1:54 AM, ICES, UT Austin
    pure function getInvMatFromCholFac(nd,CholeskyLower,Diagonal) result(InvMatFromCholFac)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getInvMatFromCholFac
#endif
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: CholeskyLower(nd,nd),Diagonal(nd)
        real(RK)                :: InvMatFromCholFac(nd,nd)
        real(RK)                :: summ
        integer(IK)             :: i,j,k
        if (nd==1_IK) then
            InvMatFromCholFac(1,1) = 1._RK / Diagonal(1)**2
            return
        end if
        InvMatFromCholFac = 0._RK
        do j=1,nd-1
            do i=j+1,nd
                InvMatFromCholFac(i,j) = CholeskyLower(i,j)
            end do
        end do
        do i = 1,nd
            InvMatFromCholFac(i,i) = 1._RK / Diagonal(i)
            do j = i+1,nd
                summ = 0._RK
                do k = i,j-1
                    summ = summ - InvMatFromCholFac(j,k) * InvMatFromCholFac(k,i)
                end do
                InvMatFromCholFac(j,i) = summ / Diagonal(j)
            end do
        end do
        do i = 1,nd
            do j = i,nd
                InvMatFromCholFac(j,i) = dot_product(InvMatFromCholFac(j:nd,j),InvMatFromCholFac(j:nd,i))
                InvMatFromCholFac(i,j) = InvMatFromCholFac(j,i)
            end do
            !InvMatFromCholFac(i,i:nd) = InvMatFromCholFac(i:nd,i)
        end do
    end function getInvMatFromCholFac

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the inverse matrix of an input symmetric-positive-definite matrix `PosDefMat`.
    !>
    !> \param[in]   nd          :   The size of the input square matrix - `nd` by `nd`.
    !> \param[in]   PosDefMat   :   The input symmetric-positive-definite matrix.
    !>
    !> \return
    !> `InvPosDefMat` : The full inverse matrix.
    !>
    !> \warning
    !> If the input matrix is not positive definite, the function will return with `InvPosDefMat(1,1) = -1._RK` to indicate
    !> the occurrence of an error. This is done to preserve the purity of the function.
    !>
    !> \warning
    !> Do not call this routine when `nd = 1`. For `nd = 1`: `InvPosDefMat = 1._RK / InvPosDefMat`.
    !>
    !> \remark
    !> On input, only the upper half of `PosDefMat` needs to be given.
    !>
    !> \remark
    !> This routine has the same functionality as the subroutine [getInvPosDefMatSqrtDet()](@ref getinvposdefmatsqrtdet),
    !> but returns the result as a function output.
    !> According to the timing tests, `compare_InvMatRoutines_1()`, the function version appears to be 15-30% faster than
    !> the subroutine version [getInvPosDefMatSqrtDet()](@ref getinvposdefmatsqrtdet).
    !>
    !> \author
    !> Amir Shahmoradi, Apr 21, 2017, 1:54 AM, ICES, UT Austin
    pure function getInvPosDefMat(nd,PosDefMat) result(InvPosDefMat)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getInvPosDefMat
#endif
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: PosDefMat(nd,nd)
        real(RK)                :: InvPosDefMat(nd,nd), CholeskyLower(nd,nd)
        real(RK)                :: Diagonal(nd)
        real(RK)                :: summ
        integer(IK)             :: i,j,k
       !if (nd==1) then
       !    InvPosDefMat(1,1) = 1._RK / sqrt(PosDefMat(1,1))
       !    return
       !end if
        do j=1,nd
            do i=1,j
                CholeskyLower(i,j) = PosDefMat(i,j)
            end do
        end do
        call getCholeskyFactor(nd,CholeskyLower,Diagonal)
        if (Diagonal(1)<0._RK) then
            InvPosDefMat = -1._RK   ! error occurred: getCholeskyFactor() failed in getInvPosDefMat()
            return
        end if
        do i = 1,nd
            CholeskyLower(i,i) = 1._RK / Diagonal(i)
            do j = i+1,nd
                summ = 0._RK
                do k = i,j-1
                    summ = summ - CholeskyLower(j,k) * CholeskyLower(k,i)
                end do
                CholeskyLower(j,i) = summ / Diagonal(j)
            end do
        end do
        do i = 1,nd
            do j = i,nd
                InvPosDefMat(j,i) = dot_product(CholeskyLower(j:nd,j),CholeskyLower(j:nd,i))
                InvPosDefMat(i,j) = InvPosDefMat(j,i)
            end do
        end do
  end function getInvPosDefMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the inverse matrix `InverseMatrix` of a `nd*nd` input matrix `MatrixLU`, and its determinant, using `LU` decomposition.
    !>
    !> \param[in]       nd              : The size of the square matrix and its `LU` decomposition - `nd` by `nd`.
    !> \param[inout]    MatrixLU        : The target symmetric-positive-definite matrix. On input it is the matrix, on output it is the LU decomposition.
    !> \param[out]      InverseMatrix   : The input symmetric-positive-definite matrix.
    !> \param[out]      detInvMat       : The determinant of the inverse of the symmetric-positive-definite matrix.
    !>
    !> \warning
    !> Do not call this routine when `nd = 1`. For `nd = 1`: `InverseMatrix = detInvMat = 1._RK / InverseMatrix`.
    !>
    !> \remark
    !> This routine has the same functionality as the function [getInvMat()](@ref getinvmat),
    !> but returns the result as a function output.
    !> According to the timing tests, `compare_InvMatRoutines_1()`, the function version of this code [getInvMat()](@ref getinvmat)
    !> appears to be 15-30% faster than this subroutine version.
    !>
    !> \author
    !> Amir Shahmoradi, Oct 18, 2009, 1:54 AM, MTU
    subroutine getInvMatDet(nd,MatrixLU,InverseMatrix,detInvMat)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getInvMatDet
#endif
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK), intent(in)    :: nd
        real(RK)   , intent(inout) :: MatrixLU(nd,nd) ! On input it is the matrix, on output it is the LU decomposition.
        real(RK)   , intent(out)   :: InverseMatrix(nd,nd)
        real(RK)   , intent(out)   :: detInvMat   ! determinant of the inverse matrix
        integer(IK)                :: i,j,Permutation(nd)
        do i = 1,nd
            do j = 1,nd
                InverseMatrix(i,j) = 0._RK
            end do
            InverseMatrix(i,i) = 1._RK
        end do
        call getLU(nd,MatrixLU,Permutation,detInvMat)
        do j = 1,nd
            detInvMat = detInvMat * MatrixLU(j,j)
            call solveLinearSystem(nd,MatrixLU,Permutation,InverseMatrix(1:nd,j))
        end do
        detInvMat = 1._RK/detInvMat
  end subroutine getInvMatDet

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the inverse matrix `InverseMatrix` of a `nd*nd` input matrix `Matrix`, and its determinant, using `LU` decomposition.
    !>
    !> \param[in]       nd      :   The size of the square matrix and its `LU` decomposition - `nd` by `nd`.
    !> \param[inout]    Matrix  :   The target symmetric-positive-definite matrix.
    !>
    !> \return
    !> `InverseMatrix` : The full inverse matrix.
    !>
    !> \remark
    !> This routine has the same functionality as the subroutine [getInvMatDet()](@ref getinvmatdet),
    !> but returns the result as a function output.
    !> According to the timing tests, `compare_InvMatRoutines_1()`, this function version appears
    !> to be 15-30% faster than the subroutine equivalent [getInvMatDet()](@ref getinvmatdet).
    !>
    !> \author
    !> Amir Shahmoradi, Apr 8, 2017, 1:54 PM, MTU
    function getInvMat(nd,Matrix) result(InverseMatrix)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getInvMat
#endif
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: Matrix(nd,nd)
        real(RK)                :: InverseMatrix(nd,nd),LU(nd,nd)
        integer(IK)             :: i,j,Permutation(nd)
        real(RK)                :: parity
        do i = 1,nd
            do j = 1,nd
                InverseMatrix(i,j) = 0._RK
            end do
            InverseMatrix(i,i) = 1._RK
        end do
        LU = Matrix
        call getLU(nd,LU,Permutation,parity)
        do j = 1,nd
            call solveLinearSystem(nd,LU,Permutation,InverseMatrix(1:nd,j))
        end do
    end function getInvMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Solve the set of `nd` linear equations `Matrix * X = InputSolution`.
    !>
    !> \param[in]       nd              :   The size of the square matrix and its `LU` decomposition - `nd` by `nd`.
    !> \param[in]       LU              :   The LU factorization of the matrix, determined by the routine [getLU()](@ref getlu).
    !> \param[in]       Permutation     :   The permutation vector returned by [getLU()](@ref getlu).
    !> \param[inout]    InputSolution   :   The right-hand-side vector of length `nd`.
    !>                                      On output, it is completely overwritten by the solution to the system.
    !>
    !> \remark
    !> The input `nd`, `LU`, and `Permutation` parameters are not modified by this routine and,
    !> can be left in place for successive calls with different right-hand sides `InputSolution`.
    !>
    !> \remark
    !> This routine takes into account the possibility that `InputSolution` will begin with many zero elements,
    !> and is therefore efficient for matrix inversion.
    !>
    !> \author
    !> Amir Shahmoradi, Apr 8, 2017, 1:54 PM, MTU
    subroutine solveLinearSystem(nd, LU, Permutation, InputSolution)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: solveLinearSystem
#endif
        use Constants_mod, only: RK, IK
        integer(IK), intent(in)    :: nd
        integer(IK), intent(in)    :: Permutation(nd)
        real(RK)   , intent(in)    :: LU(nd,nd)
        real(RK)   , intent(inout) :: InputSolution(nd)
        integer(IK)                :: i,ii
        real(RK)                   :: summ
        ii = 0
        do i=1,nd
            summ = InputSolution(Permutation(i))
            InputSolution(Permutation(i)) = InputSolution(i)
            if (ii /= 0) then
                summ = summ - dot_product( LU(i,ii:i-1) , InputSolution(ii:i-1) )
            else if (summ /= 0._RK) then
                ii = i
            end if
            InputSolution(i) = summ
        end do
        do i = nd,1,-1
            InputSolution(i) = ( InputSolution(i) - dot_product( LU(i,i+1:nd) , InputSolution(i+1:nd)) ) / LU(i,i)
        end do
    end subroutine solveLinearSystem

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the LU decomposition of the input matrix `MatrixLU(nd,nd)`.
    !>
    !> \param[in]       nd          :   The size of the square matrix and its `LU` decomposition - `nd` by `nd`.
    !> \param[inout]    MatrixLU    :   The input matrix. On output, it is completely overwritten by its LU decomposition.
    !> \param[out]      Permutation :   An output vector of length `nd` that records the row permutation effected by the partial pivoting.
    !> \param[out]      Parity      :   An output real as `+-1` depending on whether the number of row interchanges was even or odd, respectively.
    !>
    !> \remark
    !> This routine is used in combination with [solveLinearSystem](@ref solvelinearsystem) to solve linear equations or invert a matrix.
    !>
    !> \author
    !> Amir Shahmoradi, Apr 21, 2017, 1:43 PM, ICES, UT Austin
    subroutine getLU(nd,MatrixLU,Permutation,parity) ! ,errorOccurred)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLU
#endif
        use Constants_mod, only: RK, IK
        use, intrinsic :: iso_fortran_env, only: output_unit
        integer(IK) , intent(in)    :: nd
        integer(IK) , intent(out)   :: Permutation(nd)
        real(RK)    , intent(inout) :: MatrixLU(nd,nd)
        real(RK)    , intent(out)   :: parity
        !logical     , intent(out)   :: errorOccurred
        real(RK)    , parameter     :: TINY = 1.e-20_RK
        real(RK)                    :: aamax,dum,summ,vv(nd)
        integer(IK)                 :: i,imax,j,k
        !errorOccurred = .false.
        parity = 1._RK
        do i = 1,nd
            aamax = 0._RK
            do j=1,nd
                if ( abs(MatrixLU(i,j)) > aamax ) aamax = abs( MatrixLU(i,j) )
            end do
            if (aamax == 0._RK) then
                write(*,*) "Statistics@getLU() failed. Singular matrix detected."
                error stop
                !errorOccurred = .true.
                !return
            end if
            vv(i) = 1._RK/aamax
        end do
        do j=1,nd
            do i=1,j-1
                summ = MatrixLU(i,j)
                do k=1,i-1
                    summ = summ - MatrixLU(i,k) * MatrixLU(k,j)
                end do
                MatrixLU(i,j) = summ
            end do
            aamax = 0._RK
            do i = j, nd
                summ = MatrixLU(i,j)
                do k = 1, j-1
                    summ = summ - MatrixLU(i,k) * MatrixLU(k,j)
                end do
                MatrixLU(i,j) = summ
                dum = vv(i) * abs(summ)
                if (dum >= aamax) then
                    imax = i
                    aamax = dum
                end if
            end do
            if (j /= imax)then
                do k=1,nd
                    dum = MatrixLU(imax,k)
                    MatrixLU(imax,k) = MatrixLU(j,k)
                    MatrixLU(j,k) = dum
                end do
                parity = - parity
                vv(imax) = vv(j)
            end if
            Permutation(j) = imax
            if (MatrixLU(j,j) == 0._RK) MatrixLU(j,j) = TINY
            if (j /= nd) then
                dum = 1._RK / MatrixLU(j,j)
                do i = j+1,nd
                    MatrixLU(i,j) = MatrixLU(i,j) * dum
                end do
            endif
        end do
  end subroutine getLU

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the product of two matrices.
    !> The product of two matrices is defined by `c(i,j) = a(i,1)*b(1,j) + a(i,2)*b(2,j) + ... + a(i,n)*b(n,j)`.
    !>
    !> \remark
    !> There is not reason to use this routine as the Fortran intrinsic already provides an optimized implementation of matrix
    !> multiplication via `matmul()`. It is present here only for archival and legacy reasons.
    !>
    !> \author
    !> Amir Shahmoradi, Oct 20, 2009, 10:56 PM, MTU
    subroutine multiplyMatrix(A,rowsA,colsA,B,rowsB,colsB,C)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: multiplyMatrix
#endif
        use Constants_mod, only: RK, IK
        use, intrinsic :: iso_fortran_env, only: output_unit
        !use Misc  , only: abortProgram
        implicit none
        integer(IK), intent(in)  :: rowsA, colsA, rowsB, colsB !Matrix Dimensions
        real(RK)   , intent(in)  :: A(rowsA,colsA)
        real(RK)   , intent(in)  :: B(rowsB,colsB)
        real(RK)   , intent(out) :: C(rowsA,colsB)
        integer(IK)              :: i,j,k,rowsC,colsC !Counters
        if (colsA /= rowsB) then
            write(*,*) 'Matrix@multiplyMatrix() failed. dimensions of matrices do not match.'
            stop
            !call abortProgram( output_unit , 1 , 1 , 'Statistics@multiplyMatrix() failed. dimensions of matrices do not match.' )
        else
            rowsC = rowsA
            colsC = colsB
        end if
        ! Initialize product matrix to 0
        do i = 1, rowsC
            do j = 1, colsC
                C(i,j) = 0._RK
            end do
        end do
        ! Find product as per above formula
        do i = 1, rowsA
            do j = 1, colsB
                do k = 1, colsA
                    C(i,j) = C(i,j) + A(i,k)*B(k,j)
                end do
            end do
        end do
    end subroutine multiplyMatrix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the determinant of a given `nd * nd` matrix via LU factorization.
    !>
    !> \param[in]       nd      :   The size of the square matrix - `nd` by `nd`.
    !> \param[in]       Matrix  :   The input matrix.
    !>
    !> \return
    !> `determinant` : The determinant of the matrix.
    !>
    !> \author
    !> Amir Shahmoradi, Oct 18, 2009, 4:10 PM, MTU
    function getDeterminant(nd,Matrix) result(determinant)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getDeterminant
#endif
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: Matrix(nd,nd)
        real(RK)                :: DummyMat(nd,nd)
        real(RK)                :: determinant
        integer(IK)             :: Permutation(nd)
        integer(IK)             :: j
        DummyMat = Matrix
        call getLU(nd,DummyMat,Permutation,determinant) ! This returns determinant as +-1.
        do j=1,nd
            determinant = determinant * DummyMat(j,j)
        end do
    end function getDeterminant

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Return the square root of the determinant of a given positive-definite `nd * nd` matrix `PosDefMat` using Cholesky factorization.
    !>
    !> \param[in]       nd          :   The size of the square matrix - `nd` by `nd`.
    !> \param[in]       PosDefMat   :   The input matrix.
    !>
    !> \return
    !> `sqrtDetPosDefMat` : The square root of the determinant of the matrix.
    !>
    !> \warning
    !> If the input matrix is not positive definite, then `sqrtDetPosDefMat = -1._RK` upon return.
    !>
    !> \author
    !> Amir Shahmoradi, Apr 21, 2017, 4:10 PM, ICES, UT
    function getSqrtDetPosDefMat(nd,PosDefMat) result(sqrtDetPosDefMat)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getSqrtDetPosDefMat
#endif
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: PosDefMat(nd,nd)
        real(RK)                :: Diagonal(nd),DummyMat(nd,nd)
        real(RK)                :: sqrtDetPosDefMat
        integer(IK)             :: i,j
        do j=1,nd
            do i=1,j
                DummyMat(i,j) = PosDefMat(i,j)
            end do
        end do
        call getCholeskyFactor(nd,DummyMat,Diagonal)
        if (Diagonal(1)<0._RK) then
            sqrtDetPosDefMat = -1._RK
            return
        end if
        sqrtDetPosDefMat = product(Diagonal)
    end function getSqrtDetPosDefMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the natural logarithm of the square root of the determinant of a given positive-definite `nd * nd` matrix `PosDefMat`
    !> using Cholesky factorization.
    !>
    !> \param[in]       nd                  :   The size of the square matrix - `nd` by `nd`.
    !> \param[inout]    PosDefMat           :   The input matrix. On input, the upper triangle should be given, which remains intact on output.
    !> \param[out]      logSqrtDetPosDefMat :   The natural logarithm of the square root of the determinant of `PosDefMat`.
    !> \param[out]      failed              :   A logical value. If `.true.`, the determinant computation has failed.
    !>
    !> \author
    !> Amir Shahmoradi, Apr 21, 2017, 4:10 PM, ICES, UT
    pure subroutine getLogSqrtDetPosDefMat(nd,PosDefMat,logSqrtDetPosDefMat,failed)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSqrtDetPosDefMat
#endif
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK), intent(in)     :: nd
        real(RK)   , intent(inout)  :: PosDefMat(nd,nd) ! Upper triangle + diagonal is input matrix, lower is output.
        real(RK)   , intent(out)    :: logSqrtDetPosDefMat
        logical    , intent(out)    :: failed
        real(RK)                    :: Diagonal(nd)
        real(RK)                    :: summ
        integer(IK)                 :: i
        failed = .false.
        do i=1,nd
            summ = PosDefMat(i,i) - dot_product(PosDefMat(i,1:i-1),PosDefMat(i,1:i-1))
            if (summ <= 0._RK) then
                failed = .true.
                return
            end if
            Diagonal(i) = sqrt(summ)
            PosDefMat(i+1:nd,i) = ( PosDefMat(i,i+1:nd) - matmul(PosDefMat(i+1:nd,1:i-1),PosDefMat(i,1:i-1)) ) / Diagonal(i)
        end do
        logSqrtDetPosDefMat = sum(log(Diagonal))
    end subroutine getLogSqrtDetPosDefMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return `.false.` value for `isPosDef`, if the Cholesky decomposition of the input matrix fails
    !> (i.e. matrix is not positive definite), otherwise return `.true.`.
    !>
    !> \param[in]       nd      :   The size of the square matrix - `nd` by `nd`.
    !> \param[inout]    Matrix  :   The input matrix.
    !>
    !> \return
    !> `isPosDef` : A logical value indicating whether the input matrix is positive-definite.
    !>
    !> \warning
    !> Do not call this routine when `nd = 1`. In such a case, if `Matrix(1,1) <= 0`, then the matrix is not positive-definite.
    !>
    !> \author
    !> Amir Shahmoradi, Apr 21, 2017, 4:10 PM, ICES, UT
    pure logical function isPosDef(nd,Matrix)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: isPosDef
#endif
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: Matrix(nd,nd)
        real(RK)                :: Array(nd,nd),p(nd)
        real(RK)                :: dummySum
        integer(IK)             :: i,j,k
        isPosDef = .true.
        Array = Matrix
        do i = 1,nd
            do j = i,nd
                dummySum = Array(i,j)
                do k = i-1,1,-1
                    dummySum = dummySum - Array(i,k) * Array(j,k)
                end do
                if(i==j)then
                    if(dummySum<=0._RK) then
                        isPosDef = .false.
                        return
                    end if
                    p(i) = sqrt(dummySum)
                else
                    Array(j,i)=dummySum/p(i)
                endif
            end do
        end do
    end function isPosDef

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the outer product of the two input matrices.
    !>
    !> \param[in]       Vector1 : The first input vector.
    !> \param[in]       Vector2 : The second input vector.
    !>
    !> \return
    !> `OuterProd` : A matrix of size `( size(Vector1), size(Vector2) )`.
    !>
    !> \author
    !> Amir Shahmoradi, Apr 21, 2017, 4:10 PM, ICES, UT
    pure function getOuterProd(Vector1,Vector2) result(OuterProd)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getOuterProd
#endif
        use Constants_mod, only: RK, IK
        real(RK), intent(in) :: Vector1(:),Vector2(:)
        real(RK)             :: OuterProd(size(Vector1), size(Vector2))
        OuterProd = spread(Vector1,dim=2,ncopies=size(Vector2)) * spread(Vector2,dim=1,ncopies=size(Vector1))
    end function getOuterProd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return an ordered matrix of the input matrix by rearranging the columns corresponding to `ColIndx` to into the corresponding
    !> columns in `ColIndxMap`, while keeping the rest of the matrix structure intact,
    !> such that if it is positive-definite, it remains positive-definite.
    !> For example,
    !> ```
    !> 11.0    12.0    13.0    14.0    15.0
    !>  0.0    22.0    23.0    24.0    25.0
    !>  0.0     0.0    33.0    34.0    35.0
    !>  0.0     0.0     0.0    44.0    45.0
    !>  0.0 *******     0.0     0.0    55.0
    !> ```
    !> goes to the following by a `ColIndx=[4]`, `ColIndxMap=[5]`,
    !> ```
    !> 11.0    12.0    13.0    15.0    14.0
    !>  0.0    22.0    23.0    25.0    24.0
    !>  0.0     0.0    33.0    35.0    34.0
    !>  0.0     0.0     0.0    55.0    45.0
    !>  0.0     0.0     0.0     0.0    44.0
    !> ```
    !>
    !> \param[in]       rank            :   The number of columns (or rows) of the input square matrix `PosDefMatUpper`.
    !> \param[in]       PosDefMatUpper  :   The input matrix.
    !> \param[in]       lenColIndx      :   The length of the input `ColIndx` vector.
    !> \param[in]       ColIndx         :   An input array of indices indicating the order by which
    !>                                      the columns in the input matrix must be rearranged.
    !> \param[in]       ColIndxMap      :   An input array of length `lenColIndxs` that indicated the
    !>                                      new column index of the corresponding column indices in `ColIndx`.
    !>
    !> \return
    !> `SortedPosDefMatUpper` : An ordered matrix of size `( rank, rank )`.
    !>
    !> \author
    !> Amir Shahmoradi, Apr 21, 2017, 4:10 PM, ICES, UT
    pure function sortPosDefMat(rank,PosDefMatUpper,lenColIndx,ColIndx,ColIndxMap) result(SortedPosDefMatUpper)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: sortPosDefMat
#endif
        use Constants_mod, only: RK, IK
        integer(IK), intent(in) :: rank, lenColIndx, ColIndx(lenColIndx), ColIndxMap(lenColIndx)
        real(RK), intent(in)    :: PosDefMatUpper(rank,rank)
        real(RK)                :: SortedPosDefMatUpper(rank,rank)
        integer(IK)             :: iRow, iCol, iRowNew, iColNew, i
        loopIndx: do i = 1, lenColIndx
            do iCol = 1, rank
                iColNew = iCol
                if (iColNew==ColIndx(i)) then
                    iColNew = ColIndxMap(i)
                elseif (iColNew==ColIndxMap(i)) then
                    iColNew = ColIndx(i)
                end if
                do iRow = 1, iCol
                    iRowNew = iRow
                    if (iRowNew==ColIndx(i)) then
                        iRowNew = ColIndxMap(i)
                    elseif (iRowNew==ColIndxMap(i)) then
                        iRowNew = ColIndx(i)
                    end if
                    if (iRowNew>iColNew) then
                        SortedPosDefMatUpper(iRow,iCol) = PosDefMatUpper(iColNew,iRowNew)
                    else
                        SortedPosDefMatUpper(iRow,iCol) = PosDefMatUpper(iRowNew,iColNew)
                    end if
                end do
            end do
        end do loopIndx
    end function sortPosDefMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Symmetrize an input upper-triangular matrix by copying the upper to the lower triangle.
    !>
    !> \param[in]       nd                  :   The size of the input matrix `PosDefMatUpper`.
    !> \param[inout]    UpperSquareMatrix   :   The input upper, and output full matrix. On output, the matrix is symmetrized.
    !>
    !> \author
    !> Amir Shahmoradi, Apr 21, 2017, 4:10 PM, ICES, UT
    pure subroutine symmetrizeUpperSquareMatrix(nd,UpperSquareMatrix)
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK) , intent(in)    :: nd
        real(RK)    , intent(inout) :: UpperSquareMatrix(nd,nd)
        integer(IK)                 :: i, j
        do j = 1, nd
            do i = 1,j-1
                UpperSquareMatrix(j,i) = UpperSquareMatrix(i,j)
            end do
        end do
    end subroutine symmetrizeUpperSquareMatrix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the the Regression Coefficient Matrix, whose dimension is `rankS11` by `rankS22`, as well as optionally the
    !> Schur complement of the `S22` block of rank `rankS22` of the input matrix `PosDefMat` of rank `rankPDM`.
    !> For example, if,
    !> \f{equation}{
    !> \texttt{PosDefMat} = \begin{pmatrix}
    !>                          S11 & S12 \\
    !>                          S21 & S22
    !>                      \end{pmatrix}
    !> \f}
    !> then, the Regression Coefficient Matrix given `S22` is: \f$ S12 \times S22^{-1} \f$, whose rank is `rankS11` by `rankS22`.
    !> The Schur complement of `S22` is:
    !> \f{equation}
    !>     \texttt{SchurComplement} = S11 - S12 \times S22^{-1} \times S21 ~,
    !> \f}
    !> whose rank is `rankS11` by `rankS11`.
    !>
    !> \warning
    !> On input, the full matrix PosDefMat must be given.
    !>
    !> \remark
    !> For clarity, note that `rankS11 + rankS22 = rankPDM`.
    !>
    !> \author
    !> Amir Shahmoradi, Apr 21, 2017, 4:10 PM, ICES, UT
    pure subroutine getRegresCoef(rankPDM,rankS11,rankS22,PosDefMat,RegresCoefMat,SchurComplement)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRegresCoef
#endif
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK), intent(in)             :: rankPDM, rankS11, rankS22
        real(RK), intent(in)                :: PosDefMat(rankPDM,rankPDM)
        real(RK), intent(out)               :: RegresCoefMat(rankS11,rankS22)
        real(RK), intent(out), optional     :: SchurComplement(rankS11,rankS11)
        real(RK)                            :: InvS22(rankS22,rankS22), S22(rankS22,rankS22)
        integer(IK)                         :: startS22
        startS22 = rankS11 + 1
        S22 = PosDefMat(startS22:rankPDM,startS22:rankPDM)
        if (rankS22==1) then
            InvS22(1,1) = 1._RK/S22(1,1)
        else
            InvS22 = getInvPosDefMat(rankS22,S22)
        end if
        if (InvS22(1,1)<0._RK) then
            RegresCoefMat(1,1) = -1._RK
            return
        end if
        RegresCoefMat = matmul( PosDefMat(1:rankS11,startS22:rankPDM) , InvS22 )
        if (present(SchurComplement)) then
            SchurComplement = PosDefMat(1:rankS11,1:rankS11) - matmul( RegresCoefMat , PosDefMat(startS22:rankPDM,1:rankS11) )
        end if
    end subroutine getRegresCoef

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Matrix_mod