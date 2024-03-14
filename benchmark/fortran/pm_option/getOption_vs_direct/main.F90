program benchmark

    use iso_fortran_env, only: error_unit
    use pm_bench, only: bench_type
    use pm_kind, only: IK, RK, SK

    implicit none

    integer(IK)                         :: i                    !<  The procedure benchmark counter.
    integer(IK)                         :: isize                !<  The array size counter.
    integer(IK)                         :: fileUnit             !<  The output file unit for benchmark results.
    integer(IK) , parameter             :: NSIZE = 14_IK        !<  The number of benchmark sizes.
    integer(IK) , parameter             :: NBENCH = 2_IK        !<  The number of benchmark procedures.
    integer(IK)                         :: arraySize(0:NSIZE)   !<  The sizes of the benchmark array.
    real(RK)                            :: dummy = 0._RK        !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    real(RK)                            :: unifrnd              !<  The uniform random number.
    real(RK)    , allocatable           :: Array(:)             !<  The array to be passed to getOption.
    real(RK)    , allocatable           :: Default(:)           !<  The array of default values.
    type(bench_type)                    :: bench(NBENCH)        !<  The Benchmark array.

    bench(1) = bench_type(name = SK_"getOption", exec = getOption , overhead = setOverhead)
    bench(2) = bench_type(name = SK_"direct", exec = direct , overhead = setOverhead)

    arraySize = [0_IK, ( 2_IK**isize, isize = 0_IK, NSIZE - 1_IK)]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "getOption() vs. direct()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "arraySize", (bench(i)%name, i = 1, NBENCH)

        loopOverArraySize: do isize = 0_IK, NSIZE

            write(*,"(*(g0,:,' '))") "Benchmarking with size", arraySize(isize)

            allocate(Array(arraySize(isize)), Default(arraySize(isize)))
            call random_number(Default)
            call random_number(Array)
            do i = 1, NBENCH
                bench(i)%timing = bench(i)%getTiming(minsec = 0.07_RK)
            end do
            deallocate(Array, Default)

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
        call random_number(unifrnd)
        if (arraySize(isize) > 0_IK) then
            !if (unifrnd < 0.5_RK) then
            !    call addOverhead_D1()
            !else
            !    call addOverhead_D1()
            !end if
            dummy = dummy + sum(Default)
        else
            dummy = dummy + unifrnd
            !if (unifrnd < 0.5_RK) then
            !    call addOverhead_D0()
            !else
            !    call addOverhead_D0()
            !end if
        end if
    end subroutine

    !subroutine addOverhead_D0()
    !    dummy = dummy + unifrnd
    !end subroutine
    !
    !subroutine addOverhead_D1()
    !    dummy = dummy + sum(Default)
    !end subroutine

    subroutine getOption()
        call random_number(unifrnd)
        if (arraySize(isize) > 0_IK) then
            if (unifrnd < 0.5_RK) then
                call add_option_D1(Array)
            else
                call add_option_D1()
            end if
        else
            if (unifrnd < 0.5_RK) then
                call add_option_D0(unifrnd)
            else
                call add_option_D0()
            end if
        end if
    end subroutine

    subroutine direct()
        call random_number(unifrnd)
        if (arraySize(isize) > 0_IK) then
            if (unifrnd < 0.5_RK) then
                call add_direct_D1(Array)
            else
                call add_direct_D1()
            end if
        else
            if (unifrnd < 0.5_RK) then
                call add_direct_D0(unifrnd)
            else
                call add_direct_D0()
            end if
        end if
    end subroutine

    subroutine add_option_D0(value)
        use pm_option, only: getOption
        real(RK), intent(in), optional :: value
        dummy = dummy + getOption(unifrnd, value)
    end subroutine

    subroutine add_direct_D0(value)
        real(RK), intent(in), optional :: value
        if (present(value)) then
            dummy = dummy + value
        else
            dummy = dummy + unifrnd
        end if
    end subroutine

    subroutine add_option_D1(Value)
        use pm_option, only: getOption
        real(RK), intent(in), optional :: Value(:)
        dummy = dummy + sum(getOption(Default, Value))
    end subroutine

    subroutine add_direct_D1(Value)
        real(RK), intent(in), optional :: Value(:)
        if (present(Value)) then
            dummy = dummy + sum(Value)
        else
            dummy = dummy + sum(Default)
        end if
    end subroutine

end program benchmark