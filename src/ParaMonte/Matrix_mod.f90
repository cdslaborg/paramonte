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

module Matrix_mod

    implicit none

    character(*), parameter :: MODULE_NAME = "@Matrix_mod"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ATTN: if nd=1, PosDefMat will not be touched, only sqrt(PosDefMat) will be output to Diagonal.
    ! ATTN: if Cholesky factorization fails, Diagonal(1) will be set to -1 to indicate error on return.
    ! Returns in the lower triangle of PosDefMat, the Cholesky factorization L of PosDefMat=L.L^T.
    ! On input, the upper upper triangle of PosDefMat should be given, which remains intact on output.
    pure subroutine getCholeskyFactor(nd,PosDefMat,Diagonal)
#if defined DLL_ENABLED && !defined CFI_ENABLED
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

    ! Solves the linear equation system of the form: PosDefMat.InputSolution = Intercept
    ! PosDefMat and Diagonal are the output of subroutine getCholeskyFactor (i.g., only the lower triangle of PosDefMat is used).
    subroutine solveLinearPosDefSystem(nd,PosDefMat,Diagonal,Intercept,InputSolution)
#if defined DLL_ENABLED && !defined CFI_ENABLED
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

    ! ATTN: Do not call this routine when nd = 1. There is no need to call this routine when nd = 1.
    ! This code returns the inverse matrix of a symmetric-positive-definite input matrix.
    ! which is given in the upper triangle of MatInvMat.
    ! On output MatInvMat is completely overwritten by in the inverse of the matrix.
    ! Also returns: determinant of the inverse matrix.
    ! Amir Shahmoradi, Apr 21, 2017, 1:54 AM, ICES, UT Austin
    pure subroutine getInvPosDefMatSqrtDet(nd,MatInvMat,sqrtDetInvPosDefMat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getInvPosDefMatSqrtDet
#endif
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK), intent(in)    :: nd
        real(RK)   , intent(inout) :: MatInvMat(nd,nd)           ! input: upper half is covariance matrix, output: inverse matrix
        real(RK)   , intent(out)   :: sqrtDetInvPosDefMat        ! determinant of the inverse matrix
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

    ! ATTN: There is no need to call this routine when nd = 1. For nd==1: InvMatFromCholFac = 1._RK / Diagonal(1)**2
    ! Returns the inverse of a matrix whose Cholesky Lower triangle is given in the lower part of CholeskyLower,
    ! and its diagonal elements in Diagonal.
    pure function getInvMatFromCholFac(nd,CholeskyLower,Diagonal) result(InvMatFromCholFac)
#if defined DLL_ENABLED && !defined CFI_ENABLED
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

    ! ATTN: Do not call this routine when nd = 1. Results may be wrong. There is no need to call this routine when nd = 1.
    ! exact same thing as subroutine getInvPosDefMatSqrtDet, but returns the result as a function.
    ! NOTE: according to my timing tests, compare_InvMatRoutines_1(), 
    ! the function version seems to be 15-30% faster than the subroutine version above.
    ! In order to preserve the purity of the function, if an error occurs, for example,
    ! when the input matrix is not positive definite, then the first element of the output positive-definite inverse matrix
    ! will be set to -1._RK to indicate an error occurrence.
    ! On input, only the upper half of PosDefMat needs to be given.
    ! Amir Shahmoradi, Apr 8, 2017, 1:54 PM, ICES, UT
    pure function getInvPosDefMat(nd,PosDefMat) result(InvPosDefMat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
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
       !    InvPosDefMat(1) = sqrt(PosDefMat(1))
       !    return
       !end if
        do j=1,nd
            do i=1,j
                CholeskyLower(i,j) = PosDefMat(i,j)
            end do
        end do
        call getCholeskyFactor(nd,CholeskyLower,Diagonal)
        if (Diagonal(1)<0._RK) InvPosDefMat(1,1) = -1._RK   ! error occurred: getCholeskyFactor() failed in getInvPosDefMat()
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
            InvPosDefMat(i,i) = dot_product(CholeskyLower(j:nd,j),CholeskyLower(j:nd,i))
            do j = i+1,nd
                InvPosDefMat(j,i) = dot_product(CholeskyLower(j:nd,j),CholeskyLower(j:nd,i))
                InvPosDefMat(i,j) = InvPosDefMat(j,i)
            end do
        end do
  end function getInvPosDefMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ATTN: Do not call this routine when nd = 1. Results may be wrong. There is no need to call this routine when nd = 1.
    ! This code returns the inverse matrix InverseMatrix of a n*n input matrix LU, and its determinant, using LU factorization.
    ! NOTE: according to my timing tests, compare_InvMatRoutines_1(),
    ! the function version of this code below seems to be 15-30% faster than this subroutine version.
    ! Amir Shahmoradi, Oct 18, 2009, 1:54 AM, MTU
    subroutine getInvMatDet(n,LU,InverseMatrix,detInvMat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getInvMatDet
#endif
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK), intent(in)    :: n
        real(RK)   , intent(inout) :: LU(n,n)     ! on input it is the matrix, on output it is LU decomposition
        real(RK)   , intent(out)   :: InverseMatrix(n,n)
        real(RK)   , intent(out)   :: detInvMat   ! determinant of the inverse matrix
        integer(IK)                :: i,j,Permutation(n)
        do i = 1,n
            do j = 1,n
                InverseMatrix(i,j) = 0._RK
            end do
            InverseMatrix(i,i) = 1._RK
        end do
        call getLU(n,LU,Permutation,detInvMat)
        do j = 1,n
            detInvMat = detInvMat * LU(j,j)
            call solveLinearSystem(n,LU,Permutation,InverseMatrix(1:n,j))
        end do
        detInvMat = 1._RK/detInvMat
  end subroutine getInvMatDet

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ! ATTN: Do not call this routine when nd = 1. Results may be wrong. There is no need to call this routine when nd = 1.
    ! exact same thing as subroutine getInvMatDet, but returns the result as a function.
    ! NOTE: according to my timing tests, compare_InvMatRoutines_1(),
    ! the function version seems to be 15-30% faster than the subroutine version above.
    ! Amir Shahmoradi, Apr 8, 2017, 1:54 PM, MTU
    function getInvMat(n,S)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getInvMat
#endif
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK), intent(in) :: n
        real(RK)   , intent(in) :: S(n,n)
        real(RK)                :: getInvMat(n,n),DummyMat(n,n)
        integer(IK)             :: i,j,Permutation(n)
        real(RK)                :: parity
        do i = 1,n
            do j = 1,n
                getInvMat(i,j) = 0._RK
            end do
            getInvMat(i,i) = 1._RK
        end do
        DummyMat = S
        call getLU(n,DummyMat,Permutation,parity)
        do j = 1,n
            call solveLinearSystem(n,DummyMat,Permutation,getInvMat(1:n,j))
        end do
    end function getInvMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Solves the set of n linear equations Array.X = InputSolution. Here A is input, not as the matrix A but
    ! rather as its LU decomposition, determined by the routine getLU(). 
    ! Permutation is input as the permutation vector returned by getLU().
    ! InputSolution(1:n) is input as the right-hand side vector,
    ! and returns with the solution vector X.
    ! Array, nd, and Permutation are not modiffied by this routine and,
    ! can be left in place for successive calls with different right-hand sides InputSolution.
    ! This routine takes into account the possibility that InputSolution will begin with many zero elements,
    ! so it is efficient for use in matrix inversion.
    ! Amir Shahmoradi, Oct 18, 2009, 1:54 PM, MTU
    subroutine solveLinearSystem(nd,Array,Permutation,InputSolution)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: solveLinearSystem
#endif
        use Constants_mod, only: RK, IK
        integer(IK), intent(in)    :: nd
        integer(IK), intent(in)    :: Permutation(nd)
        real(RK)   , intent(in)    :: Array(nd,nd)
        real(RK)   , intent(inout) :: InputSolution(nd)
        integer(IK)                :: i,ii
        real(RK)                   :: summ
        ii = 0
        do i=1,nd
            summ = InputSolution(Permutation(i))
            InputSolution(Permutation(i)) = InputSolution(i)
            if (ii /= 0) then
                summ = summ - dot_product( Array(i,ii:i-1) , InputSolution(ii:i-1) )
            else if (summ /= 0._RK) then
                ii = i
            end if
            InputSolution(i) = summ
        end do
        do i = nd,1,-1
            InputSolution(i) = ( InputSolution(i) - dot_product( Array(i,i+1:nd) , InputSolution(i+1:nd)) ) / Array(i,i)
        end do
    end subroutine solveLinearSystem

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Returns the LU decomposition of the input matrix Array(nd,nd).
    ! Permutation(1:n) is an output vector that records the row permutation effected by the partial pivoting
    ! Parity is output as +-1 depending on whether the number of row interchanges was even or odd, respectively.
    ! This routine is used in combination with solveLinearSystem to solve linear equations or invert a matrix.
    ! Amir Shahmoradi, Apr 21, 2017, 1:43 PM, ICES, UT
    subroutine getLU(n,Array,Permutation,parity) !,errorMessage)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLU
#endif
        use Constants_mod, only: RK, IK
        use, intrinsic :: iso_fortran_env, only: output_unit
        !use Misc  , only: abortProgram
        integer(IK), intent(in)    :: n
        integer(IK), intent(out)   :: Permutation(n)
        real(RK)   , intent(inout) :: Array(n,n)
        real(RK)   , intent(out)   :: parity
        real(RK)   , parameter     :: TINY = 1.e-20_RK
        real(RK)                   :: aamax,dum,summ,vv(n)
        integer(IK)                :: i,imax,j,k
        parity = 1._RK
        do i = 1,n
            aamax = 0._RK
            do j=1,n
                if ( abs(Array(i,j)) > aamax ) aamax = abs( Array(i,j) )
            end do
            if (aamax == 0._RK) then
                write(*,*) 'Statistics@getLU() failed. Singular matrix detected.'
                stop
                !call abortProgram( output_unit , 1 , 1 , 'Statistics@getLU() failed. Singular matrix detected.' )
            end if
            vv(i) = 1._RK/aamax
        end do
        do j=1,n
            do i=1,j-1
                summ = Array(i,j)
                do k=1,i-1
                    summ = summ - Array(i,k)*Array(k,j)
                end do
                Array(i,j) = summ
            end do
            aamax = 0._RK
            do i = j, n
                summ = Array(i,j)
                do k = 1, j-1
                    summ = summ - Array(i,k)*Array(k,j)
                end do
                Array(i,j) = summ
                dum = vv(i) * abs(summ)
                if (dum >= aamax) then
                    imax = i
                    aamax = dum
                end if
            end do
            if (j /= imax)then
                do k=1,n
                    dum = Array(imax,k)
                    Array(imax,k) = Array(j,k)
                    Array(j,k) = dum
                end do
                parity = - parity
                vv(imax) = vv(j)
            end if
            Permutation(j) = imax
            if (Array(j,j) == 0._RK) Array(j,j) = TINY
            if (j /= n) then
                dum = 1._RK / Array(j,j)
                do i = j+1,n
                    Array(i,j) = Array(i,j) * dum
                end do
            endif
        end do
  end subroutine getLU

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! subroutine to find product of two matrices
    ! Product of two matrices is defined by c(i,j) = a(i,1)*b(1,j) + a(i,2)*b(2,j)+........+a(i,n)*b(n,j)
    ! Amir Shahmoradi, Oct 20, 2009, 10:56 PM, MTU
    subroutine multiplyMatrix(A,rowsA,colsA,B,rowsB,colsB,C)
#if defined DLL_ENABLED && !defined CFI_ENABLED
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

    ! Returns the determinant of a given nd*nd matrix InputMat, using LU factorization.
    ! Amir Shahmoradi, Oct 18, 2009, 4:10 PM, MTU
    function getDeterminant(nd,InputMat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getDeterminant
#endif
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: InputMat(nd,nd)
        real(RK)                :: DummyMat(nd,nd)
        real(RK)                :: getDeterminant
        integer(IK)             :: Permutation(nd)
        integer(IK)             :: j
        DummyMat = InputMat
        call getLU(nd,DummyMat,Permutation,getDeterminant) !    This returns getDeterminant as +-1.
        do j=1,nd
            getDeterminant = getDeterminant*DummyMat(j,j)
        end do
    end function getDeterminant

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Returns the determinant of a given positive-definite nd*nd matrix PosDefMat using Cholesky factorization.
    ! Amir Shahmoradi, Apr 21, 2017, 4:10 PM, ICES, UT
    function getSqrtDetPosDefMat(nd,PosDefMat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getSqrtDetPosDefMat
#endif
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: PosDefMat(nd,nd)
        real(RK)                :: Diagonal(nd),DummyMat(nd,nd)
        real(RK)                :: getSqrtDetPosDefMat
        integer(IK)             :: i,j
        do j=1,nd
            do i=1,j
                DummyMat(i,j) = PosDefMat(i,j)
            end do
        end do
        call getCholeskyFactor(nd,DummyMat,Diagonal)
        if (Diagonal(1)<0._RK) then
            getSqrtDetPosDefMat = -1._RK
            return
        end if
        getSqrtDetPosDefMat = product(Diagonal)
    end function getSqrtDetPosDefMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ATTN: if Cholesky factorization fails, logSqrtDetPosDefMat will be set to -1 to indicate error on return.
    ! Returns in the lower triangle of PosDefMat, the Cholesky factorization L of PosDefMat=L.L^T.
    ! On input, the upper upper triangle of PosDefMat should be given, which remains intact on output.
    pure subroutine getLogSqrtDetPosDefMat(nd,PosDefMat,logSqrtDetPosDefMat,failed)
#if defined DLL_ENABLED && !defined CFI_ENABLED
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

    ! ATTN: Do not call this routine when nd = 1. Results may be wrong. There is no need to call this routine when nd = 1.
    ! Returns False value for isPosDef, if the Cholesky decomposition fails (i.e. matrix is not positive definite),
    ! otherwise isPosDef=true
    pure logical function isPosDef(nd,ArrayIn)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: isPosDef
#endif
        use Constants_mod, only: RK, IK
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: ArrayIn(nd,nd)
        real(RK)                :: Array(nd,nd),p(nd)
        real(RK)                :: dummySum
        integer(IK)             :: i,j,k
        isPosDef = .true.
        Array = ArrayIn
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

    function getOuterProd(Array1,Array2)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getOuterProd
#endif
        use Constants_mod, only: RK, IK
        real(RK), intent(in) :: Array1(:),Array2(:)
        real(RK)             :: getOuterProd(size(Array1),size(Array2))
        getOuterProd = spread(Array1,dim=2,ncopies=size(Array2)) * spread(Array2,dim=1,ncopies=size(Array1))
    end function getOuterProd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns an ordered matrix of the input matrix, by rearranging the columns corresponding to ColIndx to into the corresponding
    ! columns in ColIndxMap, while keeping the rest of the matrix structure intact.
    ! such that if it is positive-definite, it remains positive-definite.
    ! Since input matrix is symmetric, only the upper triangle of the input matrix will be used,
    ! and only the upper triangle of the output matrix will be computed.
    ! example:
    ! 11.0    12.0    13.0    14.0    15.0
    !  0.0    22.0    23.0    24.0    25.0
    !  0.0     0.0    33.0    34.0    35.0
    !  0.0     0.0     0.0    44.0    45.0
    !  0.0 *******     0.0     0.0    55.0
    ! goes to the following by a ColIndx=[4], ColIndxMap=[5] 
    ! 11.0    12.0    13.0    15.0    14.0
    !  0.0    22.0    23.0    25.0    24.0
    !  0.0     0.0    33.0    35.0    34.0
    !  0.0     0.0     0.0    55.0    45.0
    !  0.0     0.0     0.0     0.0    44.0
    pure function sortPosDefMat(rank,PosDefMatUpper,lenColIndx,ColIndx,ColIndxMap) result(SortedPosDefMatUpper)
#if defined DLL_ENABLED && !defined CFI_ENABLED
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

    ! Note: on input, the full matrix PosDefMat must be given.
    ! returns the the Regression Coefficient Matrix, whose dimension is rankS11 by rankS22, as well as optionally the
    ! Schur complement of the S22 block of rank rankS22 of the input matrix PosDefMat of rank rankPDM.
    !  For example, if 
    ! PosDefMat = | S11 S12 |
    !             | S21 S22 |
    ! then, the Regression Coefficient Matrix given S22 is: S12 * S22^(-1), whose rank is rankS11 by rankS22.
    ! The Schur complement of S22 is: SchurComplement = S11 - S12 * S22^(-1) * S21, whose rank is rankS11 by rankS11.
    ! For clarity, note that the rank of rankS11 + rankS22 = rankPDM.
    pure subroutine getRegresCoef(rankPDM,rankS11,rankS22,PosDefMat,RegresCoefMat,SchurComplement)
    !subroutine getRegresCoef(rankPDM,rankS11,rankS22,PosDefMat,RegresCoefMat,SchurComplement)
#if defined DLL_ENABLED && !defined CFI_ENABLED
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