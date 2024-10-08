! Test the performance of Cholesky factorization computation using an assumed-shape interface vs. explicit-shape interface.
program benchmark

    use pm_kind, only: IK, LK, RKG => RKD, SK
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: itry, ntry
    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: iarr                         !<  The array size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: NARR = 18_IK                 !<  The number of benchmark array sizes.
    real(RKG)       , allocatable       :: samdim1(:,:)                 !<  The positive-definite matrix.
    real(RKG)       , allocatable       :: samdim2(:,:)                 !<  The positive-definite matrix.
    type(bench_type), allocatable       :: bench(:)                     !<  The Benchmark array.
    integer(IK)     , parameter         :: nsammax = 2**NARR            !<  The maximum number of calculation repeats.
    integer(IK)     , parameter         :: ndim = 5_IK                  !<  The number of data attributes.
    real(RKG)                           :: mean(ndim)                   !<  The sample mean.
    integer(IK)                         :: isam, nsam                   !<  The benchmarking array size.
    real(RKG)                           :: dumm

    bench = [ bench_type(name = SK_"intrinsicDIM1", exec = intrinsicDIM1, overhead = setOverhead) &
            , bench_type(name = SK_"intrinsicDIM2", exec = intrinsicDIM2, overhead = setOverhead) &
            , bench_type(name = SK_"setMeanDIM1", exec = setMeanDIM1, overhead = setOverhead) &
            , bench_type(name = SK_"setMeanDIM2", exec = setMeanDIM2, overhead = setOverhead) &
            ]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "sample mean benchmarking..."
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "nsam", (bench(i)%name, i = 1, size(bench))

        dumm = 0._RKG
        loopOverMatrixSize: do iarr = 1, NARR - 1

            nsam = 2**iarr
            ntry = nsammax / nsam
            allocate(samdim1(nsam, ndim), samdim2(ndim, nsam))
            write(*,"(*(g0,:,' '))") "Benchmarking setMeanDIM1() vs. setMeanDIM2()", nsam, ntry

            do i = 1, size(bench)
                bench(i)%timing = bench(i)%getTiming()
            end do

            write(fileUnit,"(*(g0,:,','))") nsam, (bench(i)%timing%mean / ntry, i = 1, size(bench))
            deallocate(samdim1, samdim2)

        end do loopOverMatrixSize
        write(*,"(*(g0,:,' '))") dumm

    close(fileUnit)

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        do itry = 1, ntry
            call random_number(mean)
            dumm = dumm - mean(1)
        end do
    end subroutine

    subroutine setSamDIM1()
        do isam = 1, nsam
            call random_number(samdim1(isam, 1 : ndim))
        end do
    end subroutine

    subroutine setSamDIM2()
        do isam = 1, nsam
            call random_number(samdim2(1 : ndim, isam))
        end do
    end subroutine

    subroutine setMeanDIM1()
        use pm_sampleMean, only: setMean
        do itry = 1, ntry
            call setSamDIM1()
            call setMean(mean, samdim1, dim = 1_IK)
            dumm = dumm + mean(1)
        end do
    end subroutine

    subroutine setMeanDIM2()
        use pm_sampleMean, only: setMean
        do itry = 1, ntry
            call setSamDIM2()
            call setMean(mean, samdim2, dim = 2_IK)
            dumm = dumm + mean(1)
        end do
    end subroutine

    subroutine intrinsicDIM1()
        do itry = 1, ntry
            call setSamDIM1()
            mean = sum(samdim1, dim = 1) / size(samdim1, dim = 1)
            dumm = dumm + mean(1)
        end do
    end subroutine

    subroutine intrinsicDIM2()
        do itry = 1, ntry
            call setSamDIM2()
            mean = sum(samdim2, dim = 2) / size(samdim1, dim = 2)
            dumm = dumm + mean(1)
        end do
    end subroutine

end program benchmark