! Test the performance of `row-major()` vs. `column-major()` matrix multiplication.
program benchmark

    use pm_bench, only: bench_type
    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, RKG => RK, RK, SK
    use pm_distUnif, only: getUnifRand

    implicit none

    integer(IK)                         :: i                                    !<  The procedure benchmark counter.
    integer(IK)                         :: fileUnit                             !<  The output file unit for benchmark results.
    integer(IK)                         :: rank, irank                          !<  The matrix rank and its counter.
    integer(IK)     , parameter         :: NRANK = 11_IK                        !<  The number of benchmark ranks.
    real(RKG)                           :: dummySum = 0._RKG                    !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    real(RKG)       , allocatable       :: matA(:,:), matB(:), matC(:), matD(:) !<  The matrix.
    type(bench_type), allocatable       :: bench(:)                             !<  The Benchmark array.

    bench = [ bench_type(name = SK_"matmulCol", exec = matmulCol, overhead = setOverhead) &
            , bench_type(name = SK_"matmulRow", exec = matmulRow, overhead = setOverhead) &
            ]

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,', '))") "MatrixRank", (bench(i)%name, i = 1, size(bench))

        loopOverMatrixRank: do irank = 1, NRANK

            rank = 2_IK**irank
            matD = getUnifRand(0._RKG, 1._RKG, rank)
            matC = getUnifRand(0._RKG, 1._RKG, rank)
            matB = getUnifRand(0._RKG, 1._RKG, rank)
            matA = getUnifRand(0._RKG, 1._RKG, rank, rank)

            write(*,"(*(g0,:,' '))") "Benchmarking with rank", rank

            do i = 1, size(bench)
                bench(i)%timing = bench(i)%getTiming(minsec = 0.07_RK)
            end do

            write(fileUnit,"(*(g0,:,', '))") rank, (bench(i)%timing%mean, i = 1, size(bench))

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
        if (all(matD == matC)) dummySum = dummySum + matC(1) + matD(1)
    end subroutine

    subroutine matmulRow()
        matC = matmul(matA, matB)
        call getDummy()
    end subroutine

    subroutine matmulCol()
        matD = matmul(matB, matA)
        call getDummy()
    end subroutine

end program benchmark