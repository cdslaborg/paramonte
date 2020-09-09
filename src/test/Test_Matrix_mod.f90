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

module Test_Matrix_mod

    use Matrix_mod
    use Err_mod, only: Err_type
    use Test_mod, only: Test_type
    implicit none

    private
    public :: test_Matrix

    type(Test_type) :: Test

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_Matrix()

        implicit none

        Test = Test_type(moduleName=MODULE_NAME)

        call test_sortPosDefMat()
        call test_getRegresCoef()
        call Test%finalize()
#ifdef CAF_ENABLED
        sync all
#endif


    end subroutine test_Matrix

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_sortPosDefMat()
        
        use Constants_mod, only: RK, IK

        implicit none
        integer                 :: rank
        real(RK), allocatable   :: PosDefMat(:,:), OutPosDefMat(:,:), OutPosDefMatRef(:,:)
        integer                 :: i,j

        if (Test%Image%isFirst) call Test%testing("sortPosDefMat")

        rank = 5
        allocate(PosDefMat(rank,rank))
        do j = 1,rank
            do i = 1,j
                PosDefMat(i,j) = i*10 + j
            end do
        end do

        ! switch the variable 4 with 5, such that the output remains a positive-definite matrix.
        OutPosDefMat = sortPosDefMat(rank,PosDefMat,1_IK,[4],[5])

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "Original Matrix:"
            do i = 1,rank
                write(Test%outputUnit,"(*(F7.1))") (PosDefMat(i,j),j=1,rank)
            end do
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "Output Positive-Definite Matrix with variables 4 and 5 swapped:"
            do i = 1,rank
                write(Test%outputUnit,"(*(F7.1))") (OutPosDefMat(i,j),j=1,rank)
            end do
            write(Test%outputUnit,"(*(g0))")
        end if

        OutPosDefMatRef = reshape(  [ 11._RK, 12._RK, 13._RK, 15._RK, 14._RK &
                                    ,  0._RK, 22._RK, 23._RK, 25._RK, 24._RK &
                                    ,  0._RK,  0._RK, 33._RK, 35._RK, 34._RK &
                                    ,  0._RK,  0._RK,  0._RK, 55._RK, 45._RK &
                                    ,  0._RK,  0._RK,  0._RK,  0._RK, 44._RK ] &
                                    , shape = [rank,rank] )
        OutPosDefMatRef = transpose(OutPosDefMatRef)
        do j = 1, rank
            do i = 1, j
                Test%assertion = OutPosDefMatRef(i,j) == OutPosDefMat(i,j)
               !write(*,*) OutPosDefMatRef(i,j) , OutPosDefMat(i,j)
                call Test%verify()
            end do
        end do
        !call Test%skipping()

    end subroutine test_sortPosDefMat

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine test_getRegresCoef()
        
        use Constants_mod, only: RK, IK

        implicit none
        real(RK)                ::  normalizedDifference
        integer , parameter     ::  rankPDM = 4, rankS11 = 3, rankS22 = 1
        real(RK), parameter     ::  PosDefMat(rankPDM,rankPDM) = reshape( &
                                    [ 4.414182620515998_RK, 1.173760167060120_RK, 0.757607629189287_RK, 5.075277296976230_RK &
                                    , 1.173760167060120_RK, 0.866956750091570_RK, 0.310654936099342_RK, 1.621274787164182_RK &
                                    , 0.757607629189287_RK, 0.310654936099342_RK, 0.955157221699132_RK, 1.254186231887444_RK &
                                    , 5.075277296976230_RK, 1.621274787164182_RK, 1.254186231887444_RK, 6.407791808961157_RK ] &
                                    , shape=shape(PosDefMat) )
        real(RK), parameter     ::  RegresCoefMatRef(rankS11,rankS22) = reshape( &
                                    [ 0.792047783119072 &
                                    , 0.253016145889269 &
                                    , 0.195728305363088 ] &
                                    , shape=shape(RegresCoefMatRef) )
        real(RK), parameter     ::  SchurComplementRef(rankS11,rankS11) = reshape( &
                                    [ 0.3943204887314210_RK, -0.110366933940115_RK, -0.235767795395625_RK &
                                    , -0.110366933940115_RK,  0.456748052015843_RK, -0.006674430520204_RK &
                                    , -0.235767795395625_RK, -0.006674430520204_RK,  0.709677475922085_RK ] &
                                    , shape=shape(SchurComplementRef) )
        real(RK)                :: SchurComplement(rankS11,rankS11)
        real(RK)                :: RegresCoefMat(rankS11,rankS22)
        integer                 :: i,j

        if (Test%Image%isFirst) call Test%testing("getRegresCoef")

        call getRegresCoef  ( rankPDM           = rankPDM           &
                            , rankS11           = rankS11           &
                            , rankS22           = rankS22           &
                            , PosDefMat         = PosDefMat         &
                            , RegresCoefMat     = RegresCoefMat     &
                            , SchurComplement   = SchurComplement   &
                            )

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0))")
            write(Test%outputUnit,"(*(g0))") "Original Covariance Matrix:"
            do i = 1,rankPDM
                write(Test%outputUnit,"(*(F22.15))") (PosDefMat(i,j),j=1,rankPDM)
            end do
        end if

        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0))"); write(Test%outputUnit,"(*(g0))") "RegresCoefMat, RegresCoefMatRef, difference, normalizedDifference:"
        end if
        do i = 1,rankS11
            !write(Test%outputUnit,"(*(F22.15))") (RegresCoefMat(i,j),j=1,rankS22), (RegresCoefMatRef(i,j),j=1,rankS22)
            do j = 1,rankS22
                normalizedDifference = abs(RegresCoefMatRef(i,j) - RegresCoefMat(i,j)) / (RegresCoefMatRef(i,j) + RegresCoefMat(i,j))
                if (Test%isDebugMode .and. Test%Image%isFirst) then
                    write(Test%outputUnit,"(*(F22.15))" ) RegresCoefMat(i,j) &
                                                        , RegresCoefMatRef(i,j) &
                                                        , abs(RegresCoefMat(i,j)-RegresCoefMatRef(i,j)) &
                                                        , normalizedDifference
                end if
                Test%assertion = normalizedDifference < 1.e-5_RK
                call Test%verify()
            end do
        end do
    
        if (Test%isDebugMode .and. Test%Image%isFirst) then
            write(Test%outputUnit,"(*(g0))"); write(Test%outputUnit,"(*(g0))") "SchurComplement, SchurComplementRef, difference, normalizedDifference:"
        end if
        do i = 1,rankS11
            do j = 1,rankS11
                normalizedDifference = abs(SchurComplementRef(i,j) - SchurComplement(i,j)) / (SchurComplementRef(i,j) + SchurComplement(i,j))
                if (Test%isDebugMode .and. Test%Image%isFirst) then
                    write(Test%outputUnit,"(*(F22.15))" ) SchurComplement(i,j) &
                                                        , SchurComplementRef(i,j) &
                                                        , abs(SchurComplement(i,j)-SchurComplementRef(i,j)) &
                                                        , normalizedDifference
                end if
                Test%assertion = normalizedDifference < 1.e-6_RK
                call Test%verify()
            end do
        end do

        write(Test%outputUnit,"(*(g0))")

        !call Test%skipping()

    end subroutine test_getRegresCoef

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Test_Matrix_mod