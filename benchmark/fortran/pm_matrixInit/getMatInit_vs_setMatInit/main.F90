! Test the performance of `getMatInit()` vs. `setMatInit()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_bench, only: bench_type
    use pm_kind, only: IK, RKG => RK, RK, SK

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK)                         :: rank, irank                  !<  The matrix rank and its counter.
    integer(IK) , parameter             :: NRANK = 11_IK                !<  The number of benchmark ranks.
    integer(IK) , parameter             :: NBENCH = 3_IK                !<  The number of benchmark procedures.
    real(RKG)                           :: dummySum = 0._RKG            !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    real(RKG)   , allocatable           :: MatInit(:,:)                 !<  The matrix.
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"setMatInit", exec = setMatInit , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"getMatInit", exec = getMatInit , overhead = setOverhead)
    bench(3) = bench_type(name = SK_"reshape_looping", exec = getMatInit_D2_reshape , overhead = setOverhead)


    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setMatInit() vs. getMatInit() vs. getMatInit_D2_reshape()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,', '))") "MatrixRank", (bench(i)%name, i = 1, NBENCH)

        loopOverMatrixRank: do irank = 1, NRANK

            rank = 2_IK**irank
            allocate(MatInit(rank, rank))
            write(*,"(*(g0,:,' '))") "Benchmarking with rank", rank

            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.07_RK)
            end do

            write(fileUnit,"(*(g0,:,', '))") rank, (bench(i)%timing%mean, i = 1, NBENCH)
            deallocate(MatInit)

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
        dummySum = dummySum + MatInit(1,1)
    end subroutine

    subroutine setMatInit()
        block
            use pm_matrixInit, only: setMatInit, dia_type
            call setMatInit(MatInit, dia_type(), 1._RKG, rank, 0_IK, 0_IK)
            call getDummy()
        end block
    end subroutine

    subroutine getMatInit()
        block
            use pm_matrixInit, only: getMatInit, dia
            MatInit = getMatInit([rank, rank], dia, 1._RKG)
            call getDummy()
        end block
    end subroutine

    subroutine getMatInit_D2_reshape()
        MatInit = reshape_looping(rank)
        call getDummy()
    end subroutine

    pure function reshape_looping(n) result(MatInit)
        integer(IK), intent(in) :: n
        real(RKG)               :: MatInit(n,n)
        integer(IK)             :: k, j
        MatInit = reshape([1._RKG, ([(0._RKG, k = 1, n)], 1._RKG, j = 1, n - 1)], shape(MatInit))
    end function

end program benchmark