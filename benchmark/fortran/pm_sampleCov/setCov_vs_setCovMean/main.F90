! Test the performance of Cholesky factorization computation using an assumed-shape interface vs. explicit-shape interface.
program benchmark

    use pm_kind, only: IK, LK, RKG => RKD, SK
    use pm_sampleCov, only: uppDia
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: itry, ntry
    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: iarr                         !<  The array size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: NARR = 18_IK                 !<  The number of benchmark array sizes.
    real(RKG)       , allocatable       :: sample(:,:)                  !<  The positive-definite matrix.
    type(bench_type), allocatable       :: bench(:)                     !<  The Benchmark array.
    integer(IK)     , parameter         :: nsammax = 2**NARR            !<  The maximum number of calculation repeats.
    integer(IK)     , parameter         :: ndim = 5_IK, dim = 2         !<  The number of data attributes.
    real(RKG)                           :: mean(ndim), cov(ndim,ndim)   !<  The positive-definite matrix.
    integer(IK)                         :: nsam                         !<  The benchmarking array size.
    real(RKG)                           :: dumm

    bench = [ bench_type(name = SK_"setCov", exec = setCov, overhead = setOverhead) &
            , bench_type(name = SK_"setCovMean", exec = setCovMean, overhead = setOverhead) &
            ]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "sample covariance benchmarking..."
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "nsam", (bench(i)%name, i = 1, size(bench))

        dumm = 0._RKG
        loopOverMatrixSize: do iarr = 1, NARR - 1

            nsam = 2**iarr
            ntry = nsammax / nsam
            allocate(sample(ndim, nsam))
            write(*,"(*(g0,:,' '))") "Benchmarking setCov() vs. setCovMean()", nsam, ntry

            do i = 1, size(bench)
                bench(i)%timing = bench(i)%getTiming()
            end do

            write(fileUnit,"(*(g0,:,','))") nsam, (bench(i)%timing%mean / ntry, i = 1, size(bench))
            deallocate(sample)

        end do loopOverMatrixSize
        write(*,"(*(g0,:,' '))") dumm

    close(fileUnit)

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        do itry = 1, ntry
            call setSample()
        end do
    end subroutine

    subroutine setSample()
        integer(IK) :: i
        call random_number(sample)
    end subroutine

    subroutine setCov()
        block
            use pm_sampleCov, only: setCov
            use pm_sampleMean, only: setMean
            do itry = 1, ntry
                call setSample()
                call setMean(mean, sample, dim)
                call setCov(cov, uppDia, mean, sample, dim)
                dumm = dumm + cov(1,1) - mean(ndim)
            end do
        end block
    end subroutine

    subroutine setCovMean()
        block
            use pm_sampleCov, only: setCovMean
            do itry = 1, ntry
                call setSample()
                call setCovMean(cov, uppDia, mean, sample, dim, sample(1:ndim, 1))
                dumm = dumm + cov(1,1) - mean(ndim)
            end do
        end block
    end subroutine

end program benchmark