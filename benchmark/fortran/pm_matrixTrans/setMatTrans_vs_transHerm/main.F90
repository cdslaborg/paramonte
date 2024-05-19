! Test the performance of `transpose()` vs. `setMatTrans()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, RKG => RKS, RK, SK
    use pm_distUnif, only: setUnifRand
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK)                         :: rank, irank                  !<  The matrix rank and its counter.
    integer(IK) , parameter             :: NRANK = 20_IK                !<  The number of benchmark ranks.
    integer(IK) , parameter             :: NBENCH = 2_IK                !<  The number of benchmark procedures.
    complex(RKG)                        :: dummySum = 0._RKG            !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    complex(RKG), allocatable           :: matA(:,:)                    !<  The matrix.
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"setMatTrans(transHerm)", exec = setMatTrans , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"transpose(conjg())", exec = transpose , overhead = setOverhead)


    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setMatTrans() vs. transpose(conjg())"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,', '))") "MatrixRank", (bench(i)%name, i = 1, NBENCH)

        loopOverMatrixRank: do irank = 1, NRANK

            rank = 1.5**irank
            allocate(matA(rank, rank))
            write(*,"(*(g0,:,' '))") "Benchmarking with rank", rank
            call setUnifRand(matA)

            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.07_RK)
            end do

            write(fileUnit,"(*(g0,:,', '))") rank, (bench(i)%timing%mean, i = 1, NBENCH)
            deallocate(matA)

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
        dummySum = dummySum + matA(1,1)
    end subroutine

    subroutine setMatTrans()
        block
            use pm_matrixTrans, only: setMatTrans, transHerm
            call setMatTrans(matA, operation = transHerm)
            call getDummy()
        end block
    end subroutine

    subroutine transpose()
        block
            intrinsic :: transpose
            matA = transpose(conjg(matA))
            call getDummy()
        end block
    end subroutine

end program benchmark