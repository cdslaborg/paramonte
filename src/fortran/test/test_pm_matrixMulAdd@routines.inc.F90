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
!>  This include file contains procedure implementations of the tests of [pm_matrixMulAdd](@ref pm_matrixMulAdd).
!>
!>  \fintest
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     IK_ENABLED
#define GET_CONJG(X) X
#define TYPE_KIND integer(IKC)
        TYPE_KIND   , parameter     :: lb = -10, ub = 10
        TYPE_KIND   , parameter     :: TOL = 0
#elif   CK_ENABLED
#define GET_CONJG(X) conjg(X)
#define TYPE_KIND complex(CKC)
        real(CKC)   , parameter     :: TOL = epsilon(0._CKC)**.66
        TYPE_KIND   , parameter     :: lb = (-1._CKC, -1._CKC), ub = (1._CKC, 1._CKC)
#elif   RK_ENABLED
#define GET_CONJG(X) X
#define TYPE_KIND real(RKC)
        real(RKC)   , parameter     :: TOL = epsilon(0._RKC)**.66
        TYPE_KIND   , parameter     :: lb = -1._RKC, ub = 1._RKC
#else
#error  "Unrecognized interface."
#endif
        logical(LK) :: renabled
        integer(IK) :: itry, begB, endB, begC, endC, irow, nrow, ncol, ndum, ndim, roffA, coffA, roffB, coffB, roffC, coffC, incB, incC
        integer(IK) , parameter     :: increments(*) = [(irow, irow = -3, -1), (irow, irow = 1, 3)]
        TYPE_KIND   , allocatable   :: matD(:,:), matT(:,:), matA(:,:), matB(:,:), matC(:,:), matO(:,:), matR(:,:)
        TYPE_KIND   , allocatable   :: vecA(:), vecB(:), vecC(:), vecO(:), vecR(:)
        TYPE_KIND   , allocatable   :: choices(:)
        TYPE_KIND   , parameter     :: ONE = 1, ZERO = 0
        TYPE_KIND                   :: alpha, beta

        assertion = .true._LK

        do itry = 1, 100

            ! Build the matrices.

            incB = getChoice(increments)
            incC = getChoice(increments)
            roffA = getUnifRand(0_IK, 3_IK)
            coffA = getUnifRand(0_IK, 3_IK)
            roffB = getUnifRand(0_IK, 3_IK)
            coffB = getUnifRand(0_IK, 3_IK)
            roffC = getUnifRand(0_IK, 3_IK)
            coffC = getUnifRand(0_IK, 3_IK)
            nrow = getUnifRand(0_IK, 10_IK)
            ncol = getUnifRand(0_IK, 10_IK)
            ndum = getUnifRand(0_IK, 10_IK)
            ndim = getUnifRand(0_IK, 10_IK)

            matA = getUnifRand(lb, ub, 2 * roffA + nrow, 2 * coffA + ndum)
            matB = getUnifRand(lb, ub, 2 * roffB + ndum, 2 * coffB + ncol)
            matC = getUnifRand(lb, ub, 2 * roffC + nrow, 2 * coffC + ncol)
            !print *, nrow, ndum, incB, incC
            vecB = getUnifRand(lb, ub, max(0, 1 + (ndum - 1) * abs(incB)))
            vecC = getUnifRand(lb, ub, max(0, 1 + (nrow - 1) * abs(incC)))
            call setBegEnd()

            choices = [ZERO, ONE, getUnifRand(lb, ub)]
            alpha = getChoice(choices)
            beta = getChoice(choices)

            ! BLAS - LEVEL 2: ?GEMV

            vecO = vecC
            vecR = vecC
            renabled = nrow == 0_IK .or. ndum == 0_IK .or. (alpha == ZERO .and. beta == ONE)
            ! \bug
            ! gfortran 13.2 yields incorrect results in the last element of the output vecR in the expression below. 
            ! This happens in release mode (with heap memory in serial shared library) but not in debug mode, indicating that aggressive optimization has something to do with this potential bug.
            ! For now, it seems like converting the `merge()` to an if block resolves the error.
            !vecR(begC:endC:incC) = merge(vecO(begC:endC:incC), matmul(alpha * matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum), vecB(begB:endB:incB)) + beta * vecO(begC:endC:incC), renabled)
            if (renabled) then
                vecR(begC:endC:incC) = vecO(begC:endC:incC)
            else
                vecR(begC:endC:incC) = matmul(alpha * matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum), vecB(begB:endB:incB)) + beta * vecO(begC:endC:incC)
            end if

            if (alpha == ONE .and. beta == ONE .and. getUnifRand()) call testGEMV()
            if (alpha == ONE .and. getUnifRand()) call testGEMV(beta = beta)
            if (beta == ONE .and. getUnifRand()) call testGEMV(alpha)
            call testGEMV(alpha, beta)

            ! BLAS - LEVEL 3: ?GEMM

            matO = matC
            matR = matO
            renabled = nrow == 0_IK .or. ncol == 0_IK .or. ((alpha == ZERO .or. ndum == 0_IK) .and. beta == ONE)
            if (.not. renabled) then
                matR(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) = & ! LCOV_EXCL_LINE
                matmul(alpha * matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum), matB(roffB + 1 : roffB + ndum, coffB + 1 : coffB + ncol)) + beta * matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol)
            end if

            if (alpha == ONE .and. beta == ONE .and. getUnifRand()) call testGEMM()
            if (alpha == ONE .and. getUnifRand()) call testGEMM(beta = beta)
            if (beta == ONE .and. getUnifRand()) call testGEMM(alpha)
            call testGEMM(alpha, beta)

            ! BLAS - LEVEL 2: SHMV - Symmetric/Hermitian matrix.

            ndim = ndum
            incC = incB
            renabled = ndim == 0_IK .or. (alpha == ZERO .and. beta == ONE)
            vecC = getUnifRand(lb, ub, max(0, 1 + (ndim - 1) * abs(incC)))
            matD = getUnifRand(lb, ub, ndim, ndim)
            call setBegEnd()
#if         CK_ENABLED
            do irow = 1, ndim
                matD(irow, irow) = matD(irow, irow)%re
            end do
#endif
            if (alpha == ONE .and. beta == ONE .and. getUnifRand()) call testSHMV()
            if (alpha == ONE .and. getUnifRand()) call testSHMV(beta = beta)
            if (beta == ONE .and. getUnifRand()) call testSHMV(alpha)
            call testSHMV(alpha, beta)

            ! BLAS - LEVEL 2: ?HPMV - Symmetric/Hermitian matrix.

            ndim = ndum
            incC = incB
            renabled = ndim == 0_IK .or. (alpha == ZERO .and. beta == ONE)
            vecC = getUnifRand(lb, ub, max(0, 1 + (ndim - 1) * abs(incC)))
            matD = getUnifRand(lb, ub, ndim, ndim)
            call setBegEnd()
#if         CK_ENABLED
            do irow = 1, ndim
                matD(irow, irow) = matD(irow, irow)%re
            end do
#endif
            if (alpha == ONE .and. beta == ONE .and. getUnifRand()) call testXPMV()
            if (alpha == ONE .and. getUnifRand()) call testXPMV(beta = beta)
            if (beta == ONE .and. getUnifRand()) call testXPMV(alpha)
            call testXPMV(alpha, beta)

            ! BLAS - LEVEL 3: ?SHMM - Symmetric/Hermitian matrix.

            matO = matC
            matR = matO
            renabled = nrow == 0_IK .or. ncol == 0_IK .or. (alpha == ZERO .and. beta == ONE)
            !if (.not. renabled) then
            !    matR(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) = & ! LCOV_EXCL_LINE
            !    matmul(alpha * matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum), matB(roffB + 1 : roffB + ndum, coffB + 1 : coffB + ncol)) + beta * matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol)
            !end if
            if (alpha == ONE .and. beta == ONE .and. getUnifRand()) call testSHMM()
            if (alpha == ONE .and. getUnifRand()) call testSHMM(beta = beta)
            if (beta == ONE .and. getUnifRand()) call testSHMM(alpha)
            call testSHMM(alpha, beta)

        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine setBegEnd()
            if (incB < 0) then
                begB = size(vecB)
                endB = 1
            else
                begB = 1
                endB = size(vecB)
            end if
            if (incC < 0) then
                begC = size(vecC)
                endC = 1
            else
                begC = 1
                endC = size(vecC)
            end if
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine testGEMV(alpha, beta)
            TYPE_KIND, intent(in), optional :: alpha, beta
            logical(LK) :: simplified
            simplified = .not. (present(alpha) .and. present(beta))
            ! nothing
            vecO = vecC
            call setMatMulAdd(matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum), vecB(begB:endB:incB), vecO(begC:endC:incC), alpha, beta)
            call reportGEMV(__LINE__, simplified, operationA = "nothing")
            ! transSymm
            vecO = vecC
            ! gfortran-11 bug: cannot pass transpose temporary array the second time for `transHerm`.
            matD = transpose(matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum))
            !print *, "transSymm", shape(matD), shape(vecB(begB:endB:incB)), shape(vecO(begC:endC:incC))
            call setMatMulAdd(matD, transSymm, vecB(begB:endB:incB), vecO(begC:endC:incC), alpha, beta)
            call reportGEMV(__LINE__, simplified, operationA = "transSymm")
            ! transHerm
            vecO = vecC
            matD = GET_CONJG(matD)
            call setMatMulAdd(matD, transHerm, vecB(begB:endB:incB), vecO(begC:endC:incC), alpha, beta)
            call reportGEMV(__LINE__, simplified, operationA = "transHerm")
            ! contiguous
            if (present(alpha) .and. present(beta)) then
                ! nothing
                vecO = vecC
                call setMatMulAdd(matA, vecB, vecO, alpha, beta, nrow, ndum, roffA, coffA, incB, incC)
                call reportGEMV(__LINE__, simplified, operationA = "nothing")
                ! transSymm
                vecO = vecC
                matD = transpose(matA)
                call setMatMulAdd(matD, transSymm, vecB, vecO, alpha, beta, ndum, nrow, coffA, roffA, incB, incC)
                call reportGEMV(__LINE__, simplified, operationA = "transSymm")
                ! transHerm
                vecO = vecC
                matD = GET_CONJG(matD)
                call setMatMulAdd(matD, transHerm, vecB, vecO, alpha, beta, ndum, nrow, coffA, roffA, incB, incC)
                call reportGEMV(__LINE__, simplified, operationA = "transHerm")
            end if
        end subroutine

        subroutine reportGEMV(line, simplified, operationA)
            integer, intent(in) :: line
            character(*), intent(in) :: operationA
            logical(LK) , intent(in) :: simplified
            logical(LK), allocatable :: tolerable(:)
            type(display_type) :: disp
#if         IK_ENABLED
            tolerable = vecR(begC:endC:incC) == vecO(begC:endC:incC)
#elif       CK_ENABLED || RK_ENABLED
            tolerable = isClose(vecR, vecO, abstol = TOL)
#else
#error      "Unrecognized interface."
#endif
            assertion = assertion .and. all(tolerable)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call disp%show("GEMV")
                call disp%show("renabled")
                call disp%show( renabled )
                call disp%show("simplified")
                call disp%show( simplified )
                call disp%show("operationA")
                call disp%show( operationA )
                call disp%show("shape(matA)")
                call disp%show( shape(matA) )
                call disp%show("shape(vecB)")
                call disp%show( shape(vecB) )
                call disp%show("shape(vecC)")
                call disp%show( shape(vecC) )
                call disp%show("[nrow, ndum]")
                call disp%show( [nrow, ndum] )
                call disp%show("[roffA, coffA, roffB, coffB, roffC, coffC, incB, incC]")
                call disp%show( [roffA, coffA, roffB, coffB, roffC, coffC, incB, incC] )
                call disp%show("[alpha, beta]")
                call disp%show( [alpha, beta] )
                call disp%show("matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum)")
                call disp%show( matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum) )
                call disp%show("vecB(begB:endB:incB)")
                call disp%show( vecB(begB:endB:incB) )
                call disp%show("vecC(begC:endC:incC)")
                call disp%show( vecC(begC:endC:incC) )
                call disp%show("vecO(begC:endC:incC)")
                call disp%show( vecO(begC:endC:incC) )
                call disp%show("vecR(begC:endC:incC)")
                call disp%show( vecR(begC:endC:incC) )
                call disp%show("vecB")
                call disp%show( vecB )
                call disp%show("vecC")
                call disp%show( vecC )
                call disp%show("vecO")
                call disp%show( vecO )
                call disp%show("vecR")
                call disp%show( vecR )
                call disp%show("tolerable")
                call disp%show( tolerable )
                call disp%show("pack(vecO, .not. tolerable)")
                call disp%show( pack(vecO, .not. tolerable) )
                call disp%show("pack(vecR, .not. tolerable)")
                call disp%show( pack(vecR, .not. tolerable) )
                call disp%show("pack(vecR - vecO, .not. tolerable)")
                call disp%show( pack(vecR - vecO, .not. tolerable) )
                call disp%show("TOL")
                call disp%show( TOL )
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"GEMV test with the above specifications must successfully pass.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine testGEMM(alpha, beta)
            TYPE_KIND, intent(in), optional :: alpha, beta
            logical(LK) :: simplified
            simplified = .not. (present(alpha) .and. present(beta))
            ! nothing, nothing
            matO = matC
            call setMatMulAdd(matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum), matB(roffB + 1 : roffB + ndum, coffB + 1 : coffB + ncol), matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol), alpha, beta)
            call reportGEMM(__LINE__, simplified, operationA = "nothing", operationB = "nothing")
            ! transSymm, nothing
            matO = matC
            matD = transpose(matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum))! gfortran-11 bug: cannot pass transpose temporary array the second time for `transHerm`.
            call setMatMulAdd(matD, transSymm, matB(roffB + 1 : roffB + ndum, coffB + 1 : coffB + ncol), matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol), alpha, beta)
            call reportGEMM(__LINE__, simplified, operationA = "transSymm", operationB = "nothing")
            ! transHerm, nothing
            matO = matC
            matD = GET_CONJG(transpose(matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum)))
            call setMatMulAdd(matD, transHerm, matB(roffB + 1 : roffB + ndum, coffB + 1 : coffB + ncol), matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol), alpha, beta)
            call reportGEMM(__LINE__, simplified, operationA = "transHerm", operationB = "nothing")
            ! nothing, transSymm
            matO = matC
            matT = transpose(matB(roffB + 1 : roffB + ndum, coffB + 1 : coffB + ncol))! gfortran-11 bug: cannot pass transpose temporary array the second time for `transHerm`.
            call setMatMulAdd(matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum), matT, transSymm, matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol), alpha, beta)
            call reportGEMM(__LINE__, simplified, operationA = "nothing", operationB = "transSymm")
            ! transSymm, transSymm
            matO = matC
            matD = transpose(matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum))! gfortran-11 bug: cannot pass transpose temporary array the second time for `transHerm`.
            matT = transpose(matB(roffB + 1 : roffB + ndum, coffB + 1 : coffB + ncol))! gfortran-11 bug: cannot pass transpose temporary array the second time for `transHerm`.
            call setMatMulAdd(matD, transSymm, matT, transSymm, matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol), alpha, beta)
            call reportGEMM(__LINE__, simplified, operationA = "transSymm", operationB = "transSymm")
            ! transHerm, transSymm
            matO = matC
            matD = GET_CONJG(transpose(matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum)))! gfortran-11 bug: cannot pass transpose temporary array the second time for `transHerm`.
            matT = transpose(matB(roffB + 1 : roffB + ndum, coffB + 1 : coffB + ncol))! gfortran-11 bug: cannot pass transpose temporary array the second time for `transHerm`.
            call setMatMulAdd(matD, transHerm, matT, transSymm, matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol), alpha, beta)
            call reportGEMM(__LINE__, simplified, operationA = "transHerm", operationB = "transSymm")
            ! nothing, transHerm
            matO = matC
            matD = matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum)
            matT = GET_CONJG(transpose(matB(roffB + 1 : roffB + ndum, coffB + 1 : coffB + ncol)))
            call setMatMulAdd(matD, matT, transHerm, matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol), alpha, beta)
            call reportGEMM(__LINE__, simplified, operationA = "transSymm", operationB = "transHerm")
            ! transSymm, transHerm
            matO = matC
            matD = transpose(matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum))
            matT = GET_CONJG(transpose(matB(roffB + 1 : roffB + ndum, coffB + 1 : coffB + ncol)))
            call setMatMulAdd(matD, transSymm, matT, transHerm, matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol), alpha, beta)
            call reportGEMM(__LINE__, simplified, operationA = "transSymm", operationB = "transHerm")
            ! transHerm, transHerm
            matO = matC
            matD = GET_CONJG(transpose(matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum)))
            matT = GET_CONJG(transpose(matB(roffB + 1 : roffB + ndum, coffB + 1 : coffB + ncol)))
            call setMatMulAdd(matD, transHerm, matT, transHerm, matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol), alpha, beta)
            call reportGEMM(__LINE__, simplified, operationA = "transHerm", operationB = "transHerm")
            ! contiguous
            if (present(alpha) .and. present(beta)) then
                ! nothing, nothing
                matO = matC
                call setMatMulAdd(matA, matB, matO, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
                call reportGEMM(__LINE__, simplified, operationA = "nothing", operationB = "nothing")
                ! transSymm, nothing
                matO = matC
                matD = transpose(matA)
                call setMatMulAdd(matD, transSymm, matB, matO, alpha, beta, nrow, ncol, ndum, coffA, roffA, roffB, coffB, roffC, coffC)
                call reportGEMM(__LINE__, simplified, operationA = "transSymm", operationB = "nothing")
                ! transHerm, nothing
                matO = matC
                matD = GET_CONJG(matD)
                call setMatMulAdd(matD, transHerm, matB, matO, alpha, beta, nrow, ncol, ndum, coffA, roffA, roffB, coffB, roffC, coffC)
                call reportGEMM(__LINE__, simplified, operationA = "transHerm", operationB = "nothing")
                ! nothing, transSymm
                matO = matC
                matT = transpose(matB)
                call setMatMulAdd(matA, matT, transSymm, matO, alpha, beta, nrow, ncol, ndum, roffA, coffA, coffB, roffB, roffC, coffC)
                call reportGEMM(__LINE__, simplified, operationA = "nothing", operationB = "transSymm")
                ! transSymm, transSymm
                matO = matC
                matD = transpose(matA)
                matT = transpose(matB)
                call setMatMulAdd(matD, transSymm, matT, transSymm, matO, alpha, beta, nrow, ncol, ndum, coffA, roffA, coffB, roffB, roffC, coffC)
                call reportGEMM(__LINE__, simplified, operationA = "transSymm", operationB = "transSymm")
                ! transHerm, transSymm
                matO = matC
                matD = GET_CONJG(transpose(matA))
                matT = transpose(matB)
                call setMatMulAdd(matD, transHerm, matT, transSymm, matO, alpha, beta, nrow, ncol, ndum, coffA, roffA, coffB, roffB, roffC, coffC)
                call reportGEMM(__LINE__, simplified, operationA = "transHerm", operationB = "transSymm")
                ! nothing, transHerm
                matO = matC
                matT = GET_CONJG(transpose(matB))
                call setMatMulAdd(matA, matT, transHerm, matO, alpha, beta, nrow, ncol, ndum, roffA, coffA, coffB, roffB, roffC, coffC)
                call reportGEMM(__LINE__, simplified, operationA = "nothing", operationB = "transHerm")
                ! transSymm, transHerm
                matO = matC
                matD = transpose(matA)
                matT = GET_CONJG(transpose(matB))
                call setMatMulAdd(matD, transSymm, matT, transHerm, matO, alpha, beta, nrow, ncol, ndum, coffA, roffA, coffB, roffB, roffC, coffC)
                call reportGEMM(__LINE__, simplified, operationA = "transSymm", operationB = "transHerm")
                ! transHerm, transHerm
                matO = matC
                matD = GET_CONJG(transpose(matA))
                matT = GET_CONJG(transpose(matB))
                call setMatMulAdd(matD, transHerm, matT, transHerm, matO, alpha, beta, nrow, ncol, ndum, coffA, roffA, coffB, roffB, roffC, coffC)
                call reportGEMM(__LINE__, simplified, operationA = "transHerm", operationB = "transHerm")
            end if
        end subroutine

        subroutine reportGEMM(line, simplified, operationA, operationB)
            integer, intent(in) :: line
            character(*), intent(in) :: operationA, operationB
            logical(LK) , intent(in) :: simplified
            logical(LK), allocatable :: tolerable(:,:)
            type(display_type) :: disp
#if         IK_ENABLED
            tolerable = matR == matO
#elif       CK_ENABLED || RK_ENABLED
            tolerable = isClose(matR, matO, abstol = TOL)
#else
#error      "Unrecognized interface."
#endif
            assertion = assertion .and. all(tolerable)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call disp%show("GEMM")
                call disp%show("renabled")
                call disp%show( renabled )
                call disp%show("simplified")
                call disp%show( simplified )
                call disp%show("operationA")
                call disp%show( operationA )
                call disp%show("operationB")
                call disp%show( operationB )
                call disp%show("shape(matA)")
                call disp%show( shape(matA) )
                call disp%show("shape(matB)")
                call disp%show( shape(matB) )
                call disp%show("shape(matC)")
                call disp%show( shape(matC) )
                call disp%show("[nrow, ndum, ncol]")
                call disp%show( [nrow, ndum, ncol] )
                call disp%show("[roffA, coffA, roffB, coffB, roffC, coffC]")
                call disp%show( [roffA, coffA, roffB, coffB, roffC, coffC] )
                call disp%show("[alpha, beta]")
                call disp%show( [alpha, beta] )
                call disp%show("matA") ! (roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum)
                call disp%show( matA ) ! (roffA + 1 : roffA + nrow, coffA + 1 : coffA + ndum)
                call disp%show("matB") ! (roffB + 1 : roffB + ndum, coffB + 1 : coffB + ncol)
                call disp%show( matB ) ! (roffB + 1 : roffB + ndum, coffB + 1 : coffB + ncol)
                call disp%show("matC(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol)")
                call disp%show( matC(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) )
                call disp%show("matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol)")
                call disp%show( matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) )
                call disp%show("matR(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol)")
                call disp%show( matR(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) )
                call disp%show("matR - matO")
                call disp%show( matR - matO )
                call disp%show("tolerable")
                call disp%show( tolerable )
                call disp%show("pack(matO, .not. tolerable)")
                call disp%show( pack(matO, .not. tolerable) )
                call disp%show("pack(matR, .not. tolerable)")
                call disp%show( pack(matR, .not. tolerable) )
                call disp%show("pack(matR - matO, .not. tolerable)")
                call disp%show( pack(matR - matO, .not. tolerable) )
                call disp%show("TOL")
                call disp%show( TOL )
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"GEMV test with the above specifications must successfully pass.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine testSHMV(alpha, beta)
            TYPE_KIND, intent(in), optional :: alpha, beta
            TYPE_KIND :: alpha_def, beta_def
            logical(LK) :: simplified
            alpha_def = 1; beta_def = 1
            if (present(beta)) beta_def = beta
            if (present(alpha)) alpha_def = alpha
            simplified = .not. (present(alpha) .and. present(beta))
            ! symmetric uppDia
            vecO = vecC
            vecR = vecC ! this initial assignment is important.
            call setMatCopy(matD, rdpack, matD, rdpack, uppDia, transSymm)
            ! \bug
            ! gfortran 13.2 yields incorrect results in the last element of the output vecR in the expression below. 
            ! This happens in release mode (with heap memory in serial shared library) but not in debug mode, indicating that aggressive optimization has something to do with this potential bug.
            ! For now, it seems like converting the `merge()` to an if block resolves the error.
            !vecR(begC:endC:incC) = merge(vecO(begC:endC:incC), matmul(alpha_def * matD, vecB(begB:endB:incB)) + beta_def * vecO(begC:endC:incC), renabled)
            if (renabled) then
                vecR(begC:endC:incC) = vecO(begC:endC:incC)
            else
                vecR(begC:endC:incC) = matmul(alpha_def * matD, vecB(begB:endB:incB)) + beta_def * vecO(begC:endC:incC)
            end if
            matT = getMatCopy(rdpack, matD, rdpack, uppDia)
            call setMatMulAdd(matT, symmetric, uppDia, vecB(begB:endB:incB), vecO(begC:endC:incC), alpha, beta)
            call reportSHMV(__LINE__, simplified, classA = "symmetric", subsetA = "uppDia")
            if (.not. simplified) then ! contiguous
                vecO = vecC
                call setRebilled(matT, fill = ZERO, lb = 1 - [roffA, coffA], ub = shape(matT))
                call setMatMulAdd(matT, symmetric, uppDia, vecB, vecO, alpha, beta, ndim, roffA, coffA, incB, incC)
                call reportSHMV(__LINE__, simplified, classA = "symmetric", subsetA = "uppDia")
            end if
            ! symmetric lowDia
            vecO = vecC
            vecR = vecC ! this initial assignment is important.
            call setMatCopy(matD, rdpack, matD, rdpack, lowDia, transSymm)
            ! \bug
            ! gfortran 13.2 yields incorrect results in the last element of the output vecR in the expression below. 
            ! This happens in release mode (with heap memory in serial shared library) but not in debug mode, indicating that aggressive optimization has something to do with this potential bug.
            ! For now, it seems like converting the `merge()` to an if block resolves the error.
            !vecR(begC:endC:incC) = merge(vecO(begC:endC:incC), matmul(alpha_def * matD, vecB(begB:endB:incB)) + beta_def * vecO(begC:endC:incC), renabled)
            if (renabled) then
                vecR(begC:endC:incC) = vecO(begC:endC:incC)
            else
                vecR(begC:endC:incC) = matmul(alpha_def * matD, vecB(begB:endB:incB)) + beta_def * vecO(begC:endC:incC)
            end if
            matT = getMatCopy(rdpack, matD, rdpack, lowDia)
            call setMatMulAdd(matT, symmetric, lowDia, vecB(begB:endB:incB), vecO(begC:endC:incC), alpha, beta)
            call reportSHMV(__LINE__, simplified, classA = "symmetric", subsetA = "lowDia")
            if (.not. simplified) then ! contiguous
                vecO = vecC
                call setRebilled(matT, fill = ZERO, lb = 1 - [roffA, coffA], ub = shape(matT))
                call setMatMulAdd(matT, symmetric, lowDia, vecB, vecO, alpha, beta, ndim, roffA, coffA, incB, incC)
                call reportSHMV(__LINE__, simplified, classA = "symmetric", subsetA = "lowDia")
            end if
            ! hermitian uppDia
            vecO = vecC
            vecR = vecC ! this initial assignment is important.
            call setMatCopy(matD, rdpack, matD, rdpack, uppDia, transHerm)
            ! \bug
            ! gfortran 13.2 yields incorrect results in the last element of the output vecR in the expression below. 
            ! This happens in release mode (with heap memory in serial shared library) but not in debug mode, indicating that aggressive optimization has something to do with this potential bug.
            ! For now, it seems like converting the `merge()` to an if block resolves the error.
            !vecR(begC:endC:incC) = merge(vecO(begC:endC:incC), matmul(alpha_def * matD, vecB(begB:endB:incB)) + beta_def * vecO(begC:endC:incC), renabled)
            if (renabled) then
                vecR(begC:endC:incC) = vecO(begC:endC:incC)
            else
                vecR(begC:endC:incC) = matmul(alpha_def * matD, vecB(begB:endB:incB)) + beta_def * vecO(begC:endC:incC)
            end if
            matT = getMatCopy(rdpack, matD, rdpack, uppDia)
            call setMatMulAdd(matT, hermitian, uppDia, vecB(begB:endB:incB), vecO(begC:endC:incC), alpha, beta)
            call reportSHMV(__LINE__, simplified, classA = "hermitian", subsetA = "uppDia")
            if (.not. simplified) then ! contiguous
                vecO = vecC
                call setRebilled(matT, fill = ZERO, lb = 1 - [roffA, coffA], ub = shape(matT))
                call setMatMulAdd(matT, hermitian, uppDia, vecB, vecO, alpha, beta, ndim, roffA, coffA, incB, incC)
                call reportSHMV(__LINE__, simplified, classA = "hermitian", subsetA = "uppDia")
            end if
            ! hermitian lowDia
            vecO = vecC
            vecR = vecC ! this initial assignment is important.
            call setMatCopy(matD, rdpack, matD, rdpack, lowDia, transHerm)
            ! \bug
            ! gfortran 13.2 yields incorrect results in the last element of the output vecR in the expression below. 
            ! This happens in release mode (with heap memory in serial shared library) but not in debug mode, indicating that aggressive optimization has something to do with this potential bug.
            ! For now, it seems like converting the `merge()` to an if block resolves the error.
            !vecR(begC:endC:incC) = merge(vecO(begC:endC:incC), matmul(alpha_def * matD, vecB(begB:endB:incB)) + beta_def * vecO(begC:endC:incC), renabled)
            if (renabled) then
                vecR(begC:endC:incC) = vecO(begC:endC:incC)
            else
                vecR(begC:endC:incC) = matmul(alpha_def * matD, vecB(begB:endB:incB)) + beta_def * vecO(begC:endC:incC)
            end if
            matT = getMatCopy(rdpack, matD, rdpack, lowDia)
            call setMatMulAdd(matT, hermitian, lowDia, vecB(begB:endB:incB), vecO(begC:endC:incC), alpha, beta)
            call reportSHMV(__LINE__, simplified, classA = "hermitian", subsetA = "lowDia")
            if (.not. simplified) then ! contiguous
                vecO = vecC
                call setRebilled(matT, fill = ZERO, lb = 1 - [roffA, coffA], ub = shape(matT))
                call setMatMulAdd(matT, hermitian, lowDia, vecB, vecO, alpha, beta, ndim, roffA, coffA, incB, incC)
                call reportSHMV(__LINE__, simplified, classA = "hermitian", subsetA = "lowDia")
            end if
        end subroutine

        subroutine reportSHMV(line, simplified, classA, subsetA)
            integer, intent(in) :: line
            character(*), intent(in) :: classA, subsetA
            logical(LK) , intent(in) :: simplified
            logical(LK), allocatable :: tolerable(:)
            type(display_type) :: disp
#if         IK_ENABLED
            tolerable = vecR == vecO
#elif       CK_ENABLED || RK_ENABLED
            tolerable = isClose(vecR, vecO, abstol = TOL)
#else
#error      "Unrecognized interface."
#endif
            assertion = assertion .and. all(tolerable)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call disp%show("SHMV")
                call disp%show("renabled")
                call disp%show( renabled )
                call disp%show("simplified")
                call disp%show( simplified )
                call disp%show("classA")
                call disp%show( classA )
                call disp%show("subsetA")
                call disp%show( subsetA )
                call disp%show("shape(matT)")
                call disp%show( shape(matT) )
                call disp%show("shape(vecB)")
                call disp%show( shape(vecB) )
                call disp%show("shape(vecC)")
                call disp%show( shape(vecC) )
                call disp%show("[roffA, coffA]")
                call disp%show( [roffA, coffA] )
                call disp%show("[ndim, incB, incC]")
                call disp%show( [ndim, incB, incC] )
                call disp%show("[alpha, beta]")
                call disp%show( [alpha, beta] )
                call disp%show("matD")
                call disp%show( matD )
                call disp%show("matT")
                call disp%show( matT )
                call disp%show("vecB(begB:endB:incB)")
                call disp%show( vecB(begB:endB:incB) )
                call disp%show("vecC(begC:endC:incC)")
                call disp%show( vecC(begC:endC:incC) )
                call disp%show("vecO(begC:endC:incC)")
                call disp%show( vecO(begC:endC:incC) )
                call disp%show("vecR(begC:endC:incC)")
                call disp%show( vecR(begC:endC:incC) )
                call disp%show("vecR(begC:endC:incC) - vecO(begC:endC:incC)")
                call disp%show( vecR(begC:endC:incC) - vecO(begC:endC:incC) )
                call disp%show("tolerable")
                call disp%show( tolerable )
                call disp%show("pack(vecO, .not. tolerable)")
                call disp%show( pack(vecO, .not. tolerable) )
                call disp%show("pack(vecR, .not. tolerable)")
                call disp%show( pack(vecR, .not. tolerable) )
                call disp%show("pack(vecR - vecO, .not. tolerable)")
                call disp%show( pack(vecR - vecO, .not. tolerable) )
                call disp%show("TOL")
                call disp%show( TOL )
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"SHMV test with the above specifications must successfully pass.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine testSHMM(alpha, beta)
            TYPE_KIND, intent(in), optional :: alpha, beta
            TYPE_KIND :: alpha_def, beta_def
            logical(LK) :: simplified
            alpha_def = 1; beta_def = 1
            if (present(beta)) beta_def = beta
            if (present(alpha)) alpha_def = alpha
            simplified = .not. (present(alpha) .and. present(beta))
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! matA symmetric uppDia
            matO = matC
            matR = matC ! this initial assignment is important.
            matA = getUnifRand(lb, ub, 2 * roffA + nrow, 2 * coffA + nrow)
            matB = getUnifRand(lb, ub, 2 * roffB + nrow, 2 * coffB + ncol)
            associate(sliceA => matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + nrow), sliceB => matB(roffB + 1 : roffB + nrow, coffB + 1 : coffB + ncol), sliceO => matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol))
                matD = sliceA
                call setMatCopy(matD, rdpack, matD, rdpack, uppDia, transSymm)
                if (.not. renabled) matR(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) = matmul(alpha_def * matD, sliceB) + beta_def * sliceO
                call setMatMulAdd(sliceA, symmetric, uppDia, sliceB, sliceO, alpha, beta)
            end associate
            call reportSHMM(__LINE__, simplified, classA = "symmetric", subsetA = "uppDia", classB = "matrix", subsetB = "uppLowDia")
            if (.not. simplified) then ! contiguous
                matO = matC
                call setMatMulAdd(matA, symmetric, uppDia, matB, matO, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
                call reportSHMM(__LINE__, simplified, classA = "symmetric", subsetA = "uppDia", classB = "matrix", subsetB = "uppLowDia")
            end if
            ! matA symmetric lowDia
            matO = matC
            matR = matC ! this initial assignment is important.
            matA = getUnifRand(lb, ub, 2 * roffA + nrow, 2 * coffA + nrow)
            matB = getUnifRand(lb, ub, 2 * roffB + nrow, 2 * coffB + ncol)
            associate(sliceA => matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + nrow), sliceB => matB(roffB + 1 : roffB + nrow, coffB + 1 : coffB + ncol), sliceO => matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol))
                matD = sliceA
                call setMatCopy(matD, rdpack, matD, rdpack, lowDia, transSymm)
                if (.not. renabled) matR(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) = matmul(alpha_def * matD, sliceB) + beta_def * sliceO
                call setMatMulAdd(sliceA, symmetric, lowDia, sliceB, sliceO, alpha, beta)
                call reportSHMM(__LINE__, simplified, classA = "symmetric", subsetA = "lowDia", classB = "matrix", subsetB = "uppLowDia")
            end associate
            if (.not. simplified) then ! contiguous
                matO = matC
                call setMatMulAdd(matA, symmetric, lowDia, matB, matO, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
                call reportSHMM(__LINE__, simplified, classA = "symmetric", subsetA = "lowDia", classB = "matrix", subsetB = "uppLowDia")
            end if
            ! matA hermitian uppDia
            matO = matC
            matR = matC ! this initial assignment is important.
            matA = getUnifRand(lb, ub, 2 * roffA + nrow, 2 * coffA + nrow)
#if         CK_ENABLED
            do irow = 1, nrow; matA(roffA + irow, coffA + irow) = matA(roffA + irow, coffA + irow)%re; end do
#endif
            matB = getUnifRand(lb, ub, 2 * roffB + nrow, 2 * coffB + ncol)
            associate(sliceA => matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + nrow), sliceB => matB(roffB + 1 : roffB + nrow, coffB + 1 : coffB + ncol), sliceO => matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol))
                matD = sliceA
                call setMatCopy(matD, rdpack, matD, rdpack, uppDia, transHerm)
                if (.not. renabled) matR(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) = matmul(alpha_def * matD, sliceB) + beta_def * sliceO
                call setMatMulAdd(sliceA, hermitian, uppDia, sliceB, sliceO, alpha, beta)
            end associate
            call reportSHMM(__LINE__, simplified, classA = "hermitian", subsetA = "uppDia", classB = "matrix", subsetB = "uppLowDia")
            if (.not. simplified) then ! contiguous
                matO = matC
                call setMatMulAdd(matA, hermitian, uppDia, matB, matO, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
                call reportSHMM(__LINE__, simplified, classA = "hermitian", subsetA = "uppDia", classB = "matrix", subsetB = "uppLowDia")
            end if
            ! matA hermitian lowDia
            matO = matC
            matR = matC ! this initial assignment is important.
            matA = getUnifRand(lb, ub, 2 * roffA + nrow, 2 * coffA + nrow)
#if         CK_ENABLED
            do irow = 1, nrow; matA(roffA + irow, coffA + irow) = matA(roffA + irow, coffA + irow)%re; end do
#endif
            matB = getUnifRand(lb, ub, 2 * roffB + nrow, 2 * coffB + ncol)
            associate(sliceA => matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + nrow), sliceB => matB(roffB + 1 : roffB + nrow, coffB + 1 : coffB + ncol), sliceO => matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol))
                matD = sliceA
                call setMatCopy(matD, rdpack, matD, rdpack, lowDia, transHerm)
                if (.not. renabled) matR(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) = matmul(alpha_def * matD, sliceB) + beta_def * sliceO
                call setMatMulAdd(sliceA, hermitian, lowDia, sliceB, sliceO, alpha, beta)
            end associate
            call reportSHMM(__LINE__, simplified, classA = "hermitian", subsetA = "lowDia", classB = "matrix", subsetB = "uppLowDia")
            if (.not. simplified) then ! contiguous
                matO = matC
                call setMatMulAdd(matA, hermitian, lowDia, matB, matO, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
                call reportSHMM(__LINE__, simplified, classA = "hermitian", subsetA = "lowDia", classB = "matrix", subsetB = "uppLowDia")
            end if
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! matB symmetric uppDia
            matO = matC
            matR = matC ! this initial assignment is important.
            matA = getUnifRand(lb, ub, 2 * roffA + nrow, 2 * coffA + ncol)
            matB = getUnifRand(lb, ub, 2 * roffB + ncol, 2 * coffB + ncol)
            associate(sliceA => matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ncol), sliceB => matB(roffB + 1 : roffB + ncol, coffB + 1 : coffB + ncol), sliceO => matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol))
                matD = sliceB
                call setMatCopy(matD, rdpack, matD, rdpack, uppDia, transSymm)
                if (.not. renabled) matR(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) = matmul(alpha_def * sliceA, matD) + beta_def * sliceO
                call setMatMulAdd(sliceA, sliceB, symmetric, uppDia, sliceO, alpha, beta)
            end associate
            call reportSHMM(__LINE__, simplified, classB = "symmetric", subsetB = "uppDia", classA = "matrix", subsetA = "uppLowDia")
            if (.not. simplified) then ! contiguous
                matO = matC
                call setMatMulAdd(matA, matB, symmetric, uppDia, matO, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
                call reportSHMM(__LINE__, simplified, classB = "symmetric", subsetB = "uppDia", classA = "matrix", subsetA = "uppLowDia")
            end if
            ! matB symmetric lowDia
            matO = matC
            matR = matC ! this initial assignment is important.
            matA = getUnifRand(lb, ub, 2 * roffA + nrow, 2 * coffA + ncol)
            matB = getUnifRand(lb, ub, 2 * roffB + ncol, 2 * coffB + ncol)
            associate(sliceA => matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ncol), sliceB => matB(roffB + 1 : roffB + ncol, coffB + 1 : coffB + ncol), sliceO => matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol))
                matD = sliceB
                call setMatCopy(matD, rdpack, matD, rdpack, lowDia, transSymm)
                if (.not. renabled) matR(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) = matmul(alpha_def * sliceA, matD) + beta_def * sliceO
                call setMatMulAdd(sliceA, sliceB, symmetric, lowDia, sliceO, alpha, beta)
                call reportSHMM(__LINE__, simplified, classB = "symmetric", subsetB = "lowDia", classA = "matrix", subsetA = "uppLowDia")
            end associate
            if (.not. simplified) then ! contiguous
                matO = matC
                call setMatMulAdd(matA, matB, symmetric, lowDia, matO, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
                call reportSHMM(__LINE__, simplified, classB = "symmetric", subsetB = "lowDia", classA = "matrix", subsetA = "uppLowDia")
            end if
            ! matB hermitian uppDia
            matO = matC
            matR = matC ! this initial assignment is important.
            matA = getUnifRand(lb, ub, 2 * roffA + nrow, 2 * coffA + ncol)
            matB = getUnifRand(lb, ub, 2 * roffB + ncol, 2 * coffB + ncol)
#if         CK_ENABLED
            do irow = 1, ncol; matB(roffB + irow, coffB + irow) = matB(roffB + irow, coffB + irow)%re; end do
#endif
            associate(sliceA => matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ncol), sliceB => matB(roffB + 1 : roffB + ncol, coffB + 1 : coffB + ncol), sliceO => matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol))
                matD = sliceB
                call setMatCopy(matD, rdpack, matD, rdpack, uppDia, transHerm)
                if (.not. renabled) matR(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) = matmul(alpha_def * sliceA, matD) + beta_def * sliceO
                call setMatMulAdd(sliceA, sliceB, hermitian, uppDia, sliceO, alpha, beta)
            end associate
            call reportSHMM(__LINE__, simplified, classB = "hermitian", subsetB = "uppDia", classA = "matrix", subsetA = "uppLowDia")
            if (.not. simplified) then ! contiguous
                matO = matC
                call setMatMulAdd(matA, matB, hermitian, uppDia, matO, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
                call reportSHMM(__LINE__, simplified, classB = "hermitian", subsetB = "uppDia", classA = "matrix", subsetA = "uppLowDia")
            end if
            ! matB hermitian lowDia
            matO = matC
            matR = matC ! this initial assignment is important.
            matA = getUnifRand(lb, ub, 2 * roffA + nrow, 2 * coffA + ncol)
            matB = getUnifRand(lb, ub, 2 * roffB + ncol, 2 * coffB + ncol)
#if         CK_ENABLED
            do irow = 1, ncol; matB(roffB + irow, coffB + irow) = matB(roffB + irow, coffB + irow)%re; end do
#endif
            associate(sliceA => matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + ncol), sliceB => matB(roffB + 1 : roffB + ncol, coffB + 1 : coffB + ncol), sliceO => matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol))
                matD = sliceB
                call setMatCopy(matD, rdpack, matD, rdpack, lowDia, transHerm)
                if (.not. renabled) matR(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) = matmul(alpha_def * sliceA, matD) + beta_def * sliceO
                call setMatMulAdd(sliceA, sliceB, hermitian, lowDia, sliceO, alpha, beta)
            end associate
            call reportSHMM(__LINE__, simplified, classB = "hermitian", subsetB = "lowDia", classA = "matrix", subsetA = "uppLowDia")
            if (.not. simplified) then ! contiguous
                matO = matC
                call setMatMulAdd(matA, matB, hermitian, lowDia, matO, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
                call reportSHMM(__LINE__, simplified, classB = "hermitian", subsetB = "lowDia", classA = "matrix", subsetA = "uppLowDia")
            end if
        end subroutine

        subroutine reportSHMM(line, simplified, classA, subsetA, classB, subsetB)
            integer, intent(in) :: line
            character(*), intent(in) :: classA, subsetA, classB, subsetB
            logical(LK) , intent(in) :: simplified
            logical(LK), allocatable :: tolerable(:,:)
            type(display_type) :: disp
#if         IK_ENABLED
            tolerable = matR == matO
#elif       CK_ENABLED || RK_ENABLED
            tolerable = isClose(matR, matO, abstol = TOL)
#else
#error      "Unrecognized interface."
#endif
            assertion = assertion .and. all(tolerable)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call disp%show("SHMM")
                call disp%show("renabled")
                call disp%show( renabled )
                call disp%show("simplified")
                call disp%show( simplified )
                call disp%show("classA")
                call disp%show( classA )
                call disp%show("subsetA")
                call disp%show( subsetA )
                call disp%show("classB")
                call disp%show( classB )
                call disp%show("subsetB")
                call disp%show( subsetB )
                call disp%show("shape(matA)")
                call disp%show( shape(matA) )
                call disp%show("shape(matB)")
                call disp%show( shape(matB) )
                call disp%show("shape(matC)")
                call disp%show( shape(matC) )
                call disp%show("[nrow, ncol]")
                call disp%show( [nrow, ncol] )
                call disp%show("[roffA, coffA, roffB, coffB, roffC, coffC]")
                call disp%show( [roffA, coffA, roffB, coffB, roffC, coffC] )
                call disp%show("[alpha, beta]")
                call disp%show( [alpha, beta] )
                call disp%show("matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + nrow)")
                call disp%show( matA(roffA + 1 : roffA + nrow, coffA + 1 : coffA + nrow) )
                call disp%show("matB(roffB + 1 : roffB + nrow, coffB + 1 : coffB + ncol)")
                call disp%show( matB(roffB + 1 : roffB + nrow, coffB + 1 : coffB + ncol) )
                call disp%show("matC(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol)")
                call disp%show( matC(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) )
                call disp%show("matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol)")
                call disp%show( matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) )
                call disp%show("matR(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol)")
                call disp%show( matR(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) )
                call disp%show("matR(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) - matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol)")
                call disp%show( matR(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) - matO(roffC + 1 : roffC + nrow, coffC + 1 : coffC + ncol) )
                call disp%show("tolerable")
                call disp%show( tolerable )
                call disp%show("pack(matO, .not. tolerable)")
                call disp%show( pack(matO, .not. tolerable) )
                call disp%show("pack(matR, .not. tolerable)")
                call disp%show( pack(matR, .not. tolerable) )
                call disp%show("pack(matR - matO, .not. tolerable)")
                call disp%show( pack(matR - matO, .not. tolerable) )
                call disp%show("TOL")
                call disp%show( TOL )
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"SHMM test with the above specifications must successfully pass.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine testXPMV(alpha, beta)
            TYPE_KIND, intent(in), optional :: alpha, beta
            TYPE_KIND :: alpha_def, beta_def
            logical(LK) :: simplified
            alpha_def = 1; beta_def = 1
            if (present(beta)) beta_def = beta
            if (present(alpha)) alpha_def = alpha
            simplified = .not. (present(alpha) .and. present(beta))
            ! symmetric uppDia
            vecO = vecC
            vecR = vecC ! this initial assignment is important.
            call setMatCopy(matD, rdpack, matD, rdpack, uppDia, transSymm)
            ! \bug
            ! gfortran 13.2 yields incorrect results in the last element of the output vecR in the expression below. 
            ! This happens in release mode (with heap memory in serial shared library) but not in debug mode, indicating that aggressive optimization has something to do with this potential bug.
            ! For now, it seems like converting the `merge()` to an if block resolves the error.
            !vecR(begC:endC:incC) = merge(vecO(begC:endC:incC), matmul(alpha_def * matD, vecB(begB:endB:incB)) + beta_def * vecO(begC:endC:incC), renabled)
            if (renabled) then
                vecR(begC:endC:incC) = vecO(begC:endC:incC)
            else
                vecR(begC:endC:incC) = matmul(alpha_def * matD, vecB(begB:endB:incB)) + beta_def * vecO(begC:endC:incC)
            end if
            vecA = getMatCopy(lfpack, matD, rdpack, uppDia)
            call setMatMulAdd(vecA, symmetric, uppDia, lfpack, vecB(begB:endB:incB), vecO(begC:endC:incC), alpha, beta)
            call reportXPMV(__LINE__, simplified, classA = "symmetric", subsetA = "uppDia")
            if (.not. simplified) then ! contiguous
                vecO = vecC
                call setMatMulAdd(vecA, symmetric, uppDia, lfpack, vecB, vecO, alpha, beta, ndim, incB, incC)
                call reportXPMV(__LINE__, simplified, classA = "symmetric", subsetA = "uppDia")
            end if
            ! symmetric lowDia
            vecO = vecC
            vecR = vecC ! this initial assignment is important.
            call setMatCopy(matD, rdpack, matD, rdpack, lowDia, transSymm)
            ! \bug
            ! gfortran 13.2 yields incorrect results in the last element of the output vecR in the expression below. 
            ! This happens in release mode (with heap memory in serial shared library) but not in debug mode, indicating that aggressive optimization has something to do with this potential bug.
            ! For now, it seems like converting the `merge()` to an if block resolves the error.
            !vecR(begC:endC:incC) = merge(vecO(begC:endC:incC), matmul(alpha_def * matD, vecB(begB:endB:incB)) + beta_def * vecO(begC:endC:incC), renabled)
            if (renabled) then
                vecR(begC:endC:incC) = vecO(begC:endC:incC)
            else
                vecR(begC:endC:incC) = matmul(alpha_def * matD, vecB(begB:endB:incB)) + beta_def * vecO(begC:endC:incC)
            end if
            vecA = getMatCopy(lfpack, matD, rdpack, lowDia)
            call setMatMulAdd(vecA, symmetric, lowDia, lfpack, vecB(begB:endB:incB), vecO(begC:endC:incC), alpha, beta)
            call reportXPMV(__LINE__, simplified, classA = "symmetric", subsetA = "lowDia")
            if (.not. simplified) then ! contiguous
                vecO = vecC
                call setMatMulAdd(vecA, symmetric, lowDia, lfpack, vecB, vecO, alpha, beta, ndim, incB, incC)
                call reportXPMV(__LINE__, simplified, classA = "symmetric", subsetA = "lowDia")
            end if
            ! hermitian uppDia
            vecO = vecC
            vecR = vecC ! this initial assignment is important.
            call setMatCopy(matD, rdpack, matD, rdpack, uppDia, transHerm)
            ! \bug
            ! gfortran 13.2 yields incorrect results in the last element of the output vecR in the expression below. 
            ! This happens in release mode (with heap memory in serial shared library) but not in debug mode, indicating that aggressive optimization has something to do with this potential bug.
            ! For now, it seems like converting the `merge()` to an if block resolves the error.
            !vecR(begC:endC:incC) = merge(vecO(begC:endC:incC), matmul(alpha_def * matD, vecB(begB:endB:incB)) + beta_def * vecO(begC:endC:incC), renabled)
            if (renabled) then
                vecR(begC:endC:incC) = vecO(begC:endC:incC)
            else
                vecR(begC:endC:incC) = matmul(alpha_def * matD, vecB(begB:endB:incB)) + beta_def * vecO(begC:endC:incC)
            end if
            vecA = getMatCopy(lfpack, matD, rdpack, uppDia)
            call setMatMulAdd(vecA, hermitian, uppDia, lfpack, vecB(begB:endB:incB), vecO(begC:endC:incC), alpha, beta)
            call reportXPMV(__LINE__, simplified, classA = "hermitian", subsetA = "uppDia")
            if (.not. simplified) then ! contiguous
                vecO = vecC
                call setMatMulAdd(vecA, hermitian, uppDia, lfpack, vecB, vecO, alpha, beta, ndim, incB, incC)
                call reportXPMV(__LINE__, simplified, classA = "hermitian", subsetA = "uppDia")
            end if
            ! hermitian lowDia
            vecO = vecC
            vecR = vecC ! this initial assignment is important.
            call setMatCopy(matD, rdpack, matD, rdpack, lowDia, transHerm)
            ! \bug
            ! gfortran 13.2 yields incorrect results in the last element of the output vecR in the expression below. 
            ! This happens in release mode (with heap memory in serial shared library) but not in debug mode, indicating that aggressive optimization has something to do with this potential bug.
            ! For now, it seems like converting the `merge()` to an if block resolves the error.
            !vecR(begC:endC:incC) = merge(vecO(begC:endC:incC), matmul(alpha_def * matD, vecB(begB:endB:incB)) + beta_def * vecO(begC:endC:incC), renabled)
            if (renabled) then
                vecR(begC:endC:incC) = vecO(begC:endC:incC)
            else
                vecR(begC:endC:incC) = matmul(alpha_def * matD, vecB(begB:endB:incB)) + beta_def * vecO(begC:endC:incC)
            end if
            vecA = getMatCopy(lfpack, matD, rdpack, lowDia)
            call setMatMulAdd(vecA, hermitian, lowDia, lfpack, vecB(begB:endB:incB), vecO(begC:endC:incC), alpha, beta)
            call reportXPMV(__LINE__, simplified, classA = "hermitian", subsetA = "lowDia")
            if (.not. simplified) then ! contiguous
                vecO = vecC
                call setMatMulAdd(vecA, hermitian, lowDia, lfpack, vecB, vecO, alpha, beta, ndim, incB, incC)
                call reportXPMV(__LINE__, simplified, classA = "hermitian", subsetA = "lowDia")
            end if
        end subroutine

        subroutine reportXPMV(line, simplified, classA, subsetA)
            character(*), intent(in) :: classA, subsetA
            logical(LK) , intent(in) :: simplified
            integer, intent(in) :: line
            logical(LK), allocatable :: tolerable(:)
            type(display_type) :: disp
#if         IK_ENABLED
            tolerable = vecR == vecO
#elif       CK_ENABLED || RK_ENABLED
            tolerable = isClose(vecR, vecO, abstol = TOL)
#else
#error      "Unrecognized interface."
#endif
            assertion = assertion .and. all(tolerable)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call disp%show("XPMV")
                call disp%show("renabled")
                call disp%show( renabled )
                call disp%show("simplified")
                call disp%show( simplified )
                call disp%show("classA")
                call disp%show( classA )
                call disp%show("subsetA")
                call disp%show( subsetA )
                call disp%show("shape(vecA)")
                call disp%show( shape(vecA) )
                call disp%show("shape(vecB)")
                call disp%show( shape(vecB) )
                call disp%show("shape(vecC)")
                call disp%show( shape(vecC) )
                call disp%show("[ndim, incB, incC]")
                call disp%show( [ndim, incB, incC] )
                call disp%show("[alpha, beta]")
                call disp%show( [alpha, beta] )
                call disp%show("matD")
                call disp%show( matD )
                call disp%show("vecA")
                call disp%show( vecA )
                call disp%show("vecB(begB:endB:incB)")
                call disp%show( vecB(begB:endB:incB) )
                call disp%show("vecC(begC:endC:incC)")
                call disp%show( vecC(begC:endC:incC) )
                call disp%show("vecO(begC:endC:incC)")
                call disp%show( vecO(begC:endC:incC) )
                call disp%show("vecR(begC:endC:incC)")
                call disp%show( vecR(begC:endC:incC) )
                call disp%show("vecR(begC:endC:incC) - vecO(begC:endC:incC)")
                call disp%show( vecR(begC:endC:incC) - vecO(begC:endC:incC) )
                call disp%show("tolerable")
                call disp%show( tolerable )
                call disp%show("pack(vecO, .not. tolerable)")
                call disp%show( pack(vecO, .not. tolerable) )
                call disp%show("pack(vecR, .not. tolerable)")
                call disp%show( pack(vecR, .not. tolerable) )
                call disp%show("pack(vecR - vecO, .not. tolerable)")
                call disp%show( pack(vecR - vecO, .not. tolerable) )
                call disp%show("TOL")
                call disp%show( TOL )
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"XPMV test with the above specifications must successfully pass.", int(line, IK))
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  TYPE_KIND
#undef  GET_CONJG