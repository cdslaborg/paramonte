! Test the performance of `setMatMulAdd()` vs. optimized BLAS.
program benchmark

    use pm_matrixInit, only: getMatInit, uppLowDia
    use pm_matrixMulAdd, only: setMatMulAdd
    use pm_distUnif, only: setUnifRand
    use iso_fortran_env, only: error_unit
    use pm_kind, only: SK, IK, RKC => RK
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                                 :: itry, ntry
    integer(IK)                                 :: ibench                       !<  The procedure benchmark counter.
    integer(IK)                                 :: irank                        !<  The matrix rank counter.
    integer(IK)                                 :: rank                         !<  The matrix rank.
    integer(IK)                                 :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK) , parameter                     :: NUM_RANK = 10_IK             !<  The number of benchmark ranks.
    integer(IK) , parameter                     :: MAX_RANK = 2**NUM_RANK       !<  The maximum possible ranks of matrices.
    integer(IK) , parameter                     :: MAX_ITER = 10000             !<  The maximum number of iterations.
    real(RKC)   , parameter                     :: ALPHA = 1._RKC               !<  The alpha parameter in the multiplication.
    real(RKC)   , parameter                     :: BETA = 1._RKC                !<  The beta  parameter in the multiplication.
    real(RKC)                                   :: dummySum = 0._RKC            !<  The dummy variable to prevent aggressive compiler optimizations.
    real(RKC)   , dimension(:,:), allocatable   :: matA, matB, matC             !<  The test matrices.
    type(bench_type)            , allocatable   :: bench(:)                     !<  The Benchmark array.

    bench = [ bench_type(name = SK_"setMatMulAddExplicit", exec = setMatMulAddExplicit, overhead = setOverhead) &
            , bench_type(name = SK_"setMatMulAddAssumed", exec = setMatMulAddAssumed, overhead = setOverhead) &
            , bench_type(name = SK_"matmul", exec = setMatMulTri, overhead = setOverhead) &
#if         BLAS_ENABLED
            , bench_type(name = SK_"GEMM", exec = setBlasGEMM, overhead = setOverhead) &
#endif
            ]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setMatMulAdd() vs. matmul() vs. GEMM()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")
        write(fileUnit, "(*(g0,:,', '))") "Matrix Rank", (bench(ibench)%name, ibench = 1, size(bench)) !, (bench(ibench)%name, ibench = 1, size(bench))
        loopOverMatDiaRank: do irank = 1, NUM_RANK

            rank = 2**irank
            ntry = MAX_ITER / rank
            write(*,"(*(g0,:,' '))") "Benchmarking with rank:", rank
            matA = getMatInit([rank, rank], uppLowDia, 0._RKC, 0._RKC, 0._RKC); !call setUnifRand(matA)
            matB = getMatInit([rank, rank], uppLowDia, 0._RKC, 0._RKC, 0._RKC); !call setUnifRand(matB)
            matC = getMatInit([rank, rank], uppLowDia, 0._RKC, 0._RKC, 0._RKC); !call setUnifRand(matC)

            ! warmup
            call setMatMulAddExplicit()
            call setMatMulAddAssumed()
            call setMatMulTri()
#if         BLAS_ENABLED
            call setBlasGEMM()
#endif

            do ibench = 1, size(bench)
                bench(ibench)%timing = bench(ibench)%getTiming(minsec = 0.07_RKC)
            end do
            write(fileUnit, "(*(g0,:,', '))") rank &
                                            , (bench(ibench)%timing%mean / ntry, ibench = 1, size(bench)) !&
                                           !, (bench(1)%timing%mean/bench(ibench)%timing%mean, ibench = 1, NBENCH)
        end do loopOverMatDiaRank
        write(*,"(*(g0,:,' '))") sum(matC)
        write(*,"(*(g0,:,' '))")
    close(fileUnit)

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        do itry = 1, ntry
            call getDummy()
        end do
    end subroutine

    subroutine getDummy()
        dummySum = dummySum + matC(1,1)
    end subroutine

    subroutine setMatMulAddExplicit()
        do itry = 1, ntry
            call setMatMulAdd(matA, matB, matC, alpha, beta, rank, rank, rank, 0_IK, 0_IK, 0_IK, 0_IK, 0_IK, 0_IK)
            call getDummy()
        end do
    end subroutine

    subroutine setMatMulAddAssumed()
        do itry = 1, ntry
            call setMatMulAdd(matA, matB, matC)
            call getDummy()
        end do
    end subroutine

    subroutine setMatMulTri()
        do itry = 1, ntry
            matC = alpha * matmul(matA(1:rank, 1:rank), matB(1:rank, 1:rank)) + beta * matC(1:rank, 1:rank)
            call getDummy()
        end do
    end subroutine

#if BLAS_ENABLED
    subroutine setBlasGEMM()
        do itry = 1, ntry
            call DGEMM  ( "N" & ! transa
                        , "N" & ! transb
                        , rank & ! m
                        , rank & ! n
                        , rank & ! k
                        , alpha & ! alpha
                        , matA & ! a
                        , rank & ! lda
                        , matB & ! b
                        , rank & ! ldb
                        , beta & ! beta
                        , matC & ! c
                        , rank & ! ldc
                        )
            call getDummy()
        end do
    end subroutine
#endif

end program benchmark