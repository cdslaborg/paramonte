#define MatB_ENABLED 0
! Test the performance of `transpose()` vs. `setMatTrans()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, RKC => RK, RK, SK
    use pm_distUnif, only: setUnifRand
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK)                         :: rank, irank                  !<  The matrix rank and its counter.
    integer(IK) , parameter             :: NRANK = 20_IK                !<  The number of benchmark ranks.
    integer(IK) , parameter             :: NBENCH = 2_IK                !<  The number of benchmark procedures.
    real(RKC)                           :: dummySum = 0._RKC            !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    real(RKC)   , allocatable           :: matA(:,:)                    !<  The matrix.
#if MatB_ENABLED
    real(RKC)   , allocatable           :: matB(:,:)                    !<  The matrix transposed.
#endif
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"setMatTrans", exec = setMatTrans , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"transpose", exec = transpose , overhead = setOverhead)


    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setMatTrans() vs. transpose()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,', '))") "MatrixRank", (bench(i)%name, i = 1, NBENCH)

        loopOverMatrixRank: do irank = 1, NRANK

            rank = 1.5**irank
            allocate(matA(rank, rank)); call setUnifRand(matA)
#if         MatB_ENABLED
            allocate(matB(rank, rank)); call setUnifRand(matB)
#endif
            write(*,"(*(g0,:,' '))") "Benchmarking with rank", rank

            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.07_RK)
            end do

            write(fileUnit,"(*(g0,:,', '))") rank, (bench(i)%timing%mean, i = 1, NBENCH)
            deallocate(matA)
#if         MatB_ENABLED
            deallocate(matB)
#endif
        end do loopOverMatrixRank
        write(*,"(*(g0,:,' '))") dummySum
        write(*,"(*(g0,:,' '))")

    close(fileUnit)

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        call getDummy()
    end subroutine

    subroutine getDummy()
#if     MatB_ENABLED
        dummySum = dummySum + matB(1,1)
#else
        dummySum = dummySum + matA(1,1)
#endif
    end subroutine

    subroutine setMatTrans()
        block
            use pm_matrixTrans, only: setMatTrans
#if         MatB_ENABLED
            call setMatTrans(matA, matB)
#else
            call setMatTrans(matA)
#endif
            call getDummy()
        end block
    end subroutine

    subroutine transpose()
        block
            intrinsic :: transpose
#if         MatB_ENABLED
            matB = transpose(matA)
#else
            matA = transpose(matA)
#endif
            call getDummy()
        end block
    end subroutine

end program benchmark