! Test the performance of `setPolyRoot()` with and without the selection `control` argument.
program benchmark

    use pm_bench, only: bench_type
    use pm_distUnif, only: setUnifRand
    use pm_arrayResize, only: setResized
    use pm_kind, only: SK, IK, LK, RKH, RK, RKG => RKD
    use iso_fortran_env, only: error_unit

    implicit none

    integer(IK)         , parameter     :: degmax = 9_IK        !<  The maximum polynomial degree.
    integer(IK)         , allocatable   :: rootCount(:)         !<  The root count.
    integer(IK)                         :: ibench               !<  The procedure benchmark counter.
    integer(IK)                         :: ideg                 !<  The array size counter.
    integer(IK)                         :: degree               !<  The array size.
    integer(IK)                         :: fileUnit             !<  The output file unit for benchmark results.
    integer(IK)                         :: rootUnit             !<  The output file unit for benchmark results.
    real(RKG)                           :: dumsum = 0._RKG      !<  The dummy computation to prevent the compiler from aggressive optimizations.
    complex(RKG)                        :: coef(0:2**degmax)    !<  The benchmark values.
    complex(RKG)                        :: root(2**degmax)      !<  The benchmark roots.
    type(bench_type)    , allocatable   :: bench(:)             !<  The Benchmark array.

    bench = [ bench_type(name = SK_"Eigen", exec = setPolyRootEigen, overhead = setOverhead) &
            , bench_type(name = SK_"Jenkins", exec = setPolyRootJenkins, overhead = setOverhead) &
            , bench_type(name = SK_"Laguerre", exec = setPolyRootLaguerre, overhead = setOverhead) &
            , bench_type(name = SK_"SGL", exec = setPolyRootSGL, overhead = setOverhead) &
            ]
    allocate(rootCount(size(bench)))

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "Benchmarking setPolyRoot()"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")
    open(newunit = rootUnit, file = "root.count", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "Polynomial Degree", (bench(ibench)%name, ibench = 1, size(bench))
        write(rootUnit, "(*(g0,:,','))") "Polynomial Degree", (bench(ibench)%name, ibench = 1, size(bench))
        loopOverArraySize: do ideg = 0, degmax

            degree = 2**ideg
            call setUnifRand(coef(0 : degree), (1._RKG, 1._RKG), (2._RKG, 2._RKG))
            write(*,"(*(g0,:,' '))") "Benchmarking with coef size", degree
            do ibench = 1, size(bench)
                bench(ibench)%timing = bench(ibench)%getTiming(minsec = 0.07_RK)
            end do
            write(rootUnit,"(*(g0,:,','))") degree, rootCount
            write(fileUnit,"(*(g0,:,','))") degree, (bench(ibench)%timing%mean, ibench = 1, size(bench))

        end do loopOverArraySize
        write(*,"(*(g0,:,' '))") dumsum
        write(*,"(*(g0,:,' '))")

    close(fileUnit)

contains

    !%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        call getDummy()
    end subroutine

    subroutine getDummy()
        dumsum = dumsum + rootCount(ibench)
    end subroutine

    subroutine setPolyRootEigen()
        use pm_polynomial, only: setPolyRoot, method => eigen
        call setPolyRoot(root(1 : degree), rootCount(ibench), coef(0 : degree), method)
        call getDummy()
    end subroutine

    subroutine setPolyRootJenkins()
        use pm_polynomial, only: setPolyRoot, method => jenkins
        call setPolyRoot(root(1 : degree), rootCount(ibench), coef(0 : degree), method)
        call getDummy()
    end subroutine

    subroutine setPolyRootLaguerre()
        use pm_polynomial, only: setPolyRoot, method => laguerre
        call setPolyRoot(root(1 : degree), rootCount(ibench), coef(0 : degree), method)
        call getDummy()
    end subroutine

    subroutine setPolyRootSGL()
        use pm_polynomial, only: setPolyRoot, method => sgl
        call setPolyRoot(root(1 : degree), rootCount(ibench), coef(0 : degree), method)
        rootCount(ibench) = root(1)%re
        call getDummy()
    end subroutine

end program benchmark