program benchmark

    use pm_bench, only: bench_type
    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, RKC => RK, RK, SK
    use pm_distUnif, only: getUnifRand

    implicit none

    integer(IK)                         :: fileUnit                             !<  The output file unit for benchmark results.
    integer(IK)                         :: i, isim, nsim                        !<  The procedure benchmark counter.
    integer(IK)                         :: rank, irank                          !<  The matrix rank and its counter.
    integer(IK)     , parameter         :: NRANK = 10_IK                        !<  The number of benchmark ranks.
    real(RKC)                           :: dummySum = 0._RKC                    !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    real(RKC)                           :: dummyOne = 0._RKC                    !<  The dummy result.
    real(RKC)                           :: dummyTwo = 0._RKC                    !<  The dummy result.
    real(RKC)       , allocatable       :: matA(:,:), matB(:)                   !<  The matrix.
    type(bench_type), allocatable       :: bench(:)                             !<  The Benchmark array.

    bench = [ bench_type(name = SK_"loop_and_dotp", exec = loop_and_dotp, overhead = setOverhead) &
            , bench_type(name = SK_"dotp_matmul", exec = dotp_matmul, overhead = setOverhead) &
            ]

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,', '))") "MatrixRank", (bench(i)%name, i = 1, size(bench))

        loopOverMatrixRank: do irank = 1, NRANK

            rank = 2_IK**irank
            nsim = nint(2.**NRANK / rank)
            matB = getUnifRand(0._RKC, 1._RKC, rank)
            matA = getUnifRand(0._RKC, 1._RKC, rank, rank)

            write(*,"(*(g0,:,' '))") "Benchmarking with rank", rank

            do i = 1, size(bench)
                bench(i)%timing = bench(i)%getTiming(minsec = 0.07_RK)
            end do

            write(fileUnit,"(*(g0,:,', '))") rank, (bench(i)%timing%mean / nsim, i = 1, size(bench))

        end do loopOverMatrixRank
        write(*,"(*(g0,:,' '))") dummySum
        write(*,"(*(g0,:,' '))")

    close(fileUnit)

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        do isim = 1, nsim
            call getDummy()
        end do
    end subroutine

    subroutine getDummy()
        dummySum = dummySum + dummyOne + dummyTwo
    end subroutine

    subroutine loop_and_dotp()
        integer(IK) :: i, sizeB
        sizeB = size(matB, 1, IK)
        dummyOne = 0._RKC
        do isim = 1, nsim
            do i = 1, sizeB
                dummyOne = dummyOne + matB(i) * dot_product(matB, matA(1:sizeB, i))
            end do
            call getDummy()
        end do
    end subroutine

    subroutine dotp_matmul()
        do isim = 1, nsim
            dummyTwo = dot_product(matB, matmul(matB, matA))
            call getDummy()
        end do
    end subroutine

end program benchmark