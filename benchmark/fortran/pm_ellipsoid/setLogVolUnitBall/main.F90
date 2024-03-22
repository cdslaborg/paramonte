! Test the performance of `setLogVolUnitBall()` with and without the selection `control` argument.
program benchmark

    use pm_bench, only: bench_type
    use pm_arraySpace, only: getLinSpace
    use pm_kind, only: IK, LK, RK1 => RKD, RK2 => RKS, RK, SK
    use iso_fortran_env, only: error_unit

    implicit none

    integer(IK)                         :: i                        !<  The procedure benchmark counter.
    integer(IK)                         :: idim                     !<  The array size counter.
    integer(IK)                         :: fileUnit                 !<  The output file unit for benchmark results.
    integer(IK)                         :: ndim_IK                  !<  The number of dimensions as `integer`.

    real(RK1)                           :: logVolUnitBall_RK1       !<  The procedure output.
    real(RK1)                           :: dum_RK1 = 0._RK1         !<  The dummy computation to prevent the compiler from aggressive optimizations.
    real(RK1)                           :: ndim_RK1                 !<  The number of dimensions as `real`.

    real(RK2)                           :: logVolUnitBall_RK2       !<  The procedure output.
    real(RK2)                           :: dum_RK2 = 0._RK2         !<  The dummy computation to prevent the compiler from aggressive optimizations.
    real(RK2)                           :: ndim_RK2                 !<  The number of dimensions as `real`.
    type(bench_type)    , allocatable   :: bench(:)                 !<  The Benchmark array.

    bench = [ bench_type(name = SK_"ndimInt2RK1", exec = ndimInt2RK1, overhead = setOverhead_RK1) &
            , bench_type(name = SK_"ndimInt2RK2", exec = ndimInt2RK2, overhead = setOverhead_RK1) &
            , bench_type(name = SK_"ndimRealRK1", exec = ndimRealRK1, overhead = setOverhead_RK1) &
            , bench_type(name = SK_"ndimRealRK2", exec = ndimRealRK2, overhead = setOverhead_RK2) &
            ]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setLogVolUnitBallInt() vs. setLogVolUnitBallReal()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "ndim", (bench(i)%name, i = 1, size(bench))
        loopOverArraySize: do idim = 0, 20_IK

            ndim_IK   = idim
            ndim_RK1 = real(idim, RK1)
            ndim_RK2 = real(idim, RK2)
            write(*,"(*(g0,:,' '))") "Benchmarking with ndim = ", idim
            do i = 1, size(bench)
                bench(i)%timing = bench(i)%getTiming(minsec = 0.04_RK)
            end do
            write(fileUnit,"(*(g0,:,','))") idim, (bench(i)%timing%mean, i = 1, size(bench))

        end do loopOverArraySize
        write(*,"(*(g0,:,' '))") dum_RK1
        write(*,"(*(g0,:,' '))")

    close(fileUnit)

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead_RK1()
        call getDummy_RK1()
    end subroutine

    subroutine setOverhead_RK2()
        call getDummy_RK2()
    end subroutine

    subroutine getDummy_RK1()
        dum_RK1 = dum_RK1 + logVolUnitBall_RK1
    end subroutine

    subroutine getDummy_RK2()
        dum_RK2 = dum_RK2 + logVolUnitBall_RK2
    end subroutine

    subroutine ndimInt2RK1()
        use pm_ellipsoid, only: setLogVolUnitBall
        call setLogVolUnitBall(logVolUnitBall_RK1, ndim_IK)
        call getDummy_RK1()
    end subroutine

    subroutine ndimInt2RK2()
        use pm_ellipsoid, only: setLogVolUnitBall
        call setLogVolUnitBall(logVolUnitBall_RK2, ndim_IK)
        call getDummy_RK2()
    end subroutine

    subroutine ndimRealRK1()
        use pm_ellipsoid, only: setLogVolUnitBall
        call setLogVolUnitBall(logVolUnitBall_RK1, ndim_RK1)
        call getDummy_RK1()
    end subroutine

    subroutine ndimRealRK2()
        use pm_ellipsoid, only: setLogVolUnitBall
        call setLogVolUnitBall(logVolUnitBall_RK2, ndim_RK2)
        call getDummy_RK2()
    end subroutine

end program benchmark