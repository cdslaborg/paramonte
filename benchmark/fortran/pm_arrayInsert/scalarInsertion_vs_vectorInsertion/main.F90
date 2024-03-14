! Test the performance of setInserted() with a vector `insertion` vs. scalar `insertion`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, LK, RK, SK
    use pm_bench, only: bench_type
    use pm_arrayRange, only: getRange

    implicit none

    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: isize                        !<  The array size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK) , parameter             :: NSIZE = 18_IK                !<  The number of benchmark ranks.
    integer(IK) , parameter             :: NBENCH = 2_IK                !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(NSIZE)             !<  The sizes of the benchmark array.
    logical(LK)                         :: dummy = .true._LK            !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    integer(IK) , allocatable           :: index(:)                     !<  The positions of insertions.
    real(RK)    , allocatable           :: array(:)                     !<  The array to insert insertion to.
    real(RK)    , allocatable           :: arrayNew(:)                  !<  The array to insert insertion to.
    real(RK)                            :: insertion(1)                 !<  The insertion to insert.
    type(bench_type)                    :: bench(NBENCH)                !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"scalarInsertion", exec = scalarInsertion , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"vectorInsertion", exec = vectorInsertion , overhead = setOverhead)

    arraySize = [( 2_IK**isize, isize = 1_IK, NSIZE )]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "scalarInsertion() vs. vectorInsertion()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverArraySize: do isize = 1, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)
            allocate(array(arraySize(isize)), arrayNew(arraySize(isize)*2))
            index = getRange(1_IK, arraySize(isize))

            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.05_RK)
            end do

            deallocate(array, arrayNew)
            write(fileUnit,"(*(g0,:,','))") arraySize(isize), (bench(i)%timing%mean, i = 1, NBENCH)

        end do loopOverArraySize
        write(*,"(*(g0,:,' '))") dummy
        write(*,"(*(g0,:,' '))")

    close(fileUnit)

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        call initialize()
        call finalize()
    end subroutine

    subroutine initialize()
        call random_number(insertion)
    end subroutine

    subroutine finalize()
        dummy = dummy .and. insertion(1) < 0.5_RK
    end subroutine

    subroutine scalarInsertion()
        use pm_arrayInsert, only: setInserted
        call initialize()
        call setInserted(arrayNew, array, insertion(1), index)
        call finalize()
    end subroutine

    subroutine vectorInsertion()
        block
            use pm_arrayInsert, only: setInserted
            call initialize()
            call setInserted(arrayNew, array, insertion, index)
            call finalize()
        end block
    end subroutine

end program benchmark