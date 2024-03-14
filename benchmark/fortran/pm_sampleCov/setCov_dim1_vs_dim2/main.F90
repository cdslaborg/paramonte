! Test the performance of Cholesky factorization computation using an assumed-shape interface vs. explicit-shape interface.
program benchmark

    use pm_kind, only: IK, LK, RKC => RKD, SK
    use pm_sampleCov, only: uppDia
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: itry, ntry
    integer(IK)                         :: i                            !<  The procedure benchmark counter.
    integer(IK)                         :: iarr                         !<  The array size counter.
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: NARR = 18_IK                 !<  The number of benchmark array sizes.
    integer(IK)     , allocatable       :: rperm(:)                     !<  The permutation vector for LUP factorization.
    real(RKC)       , allocatable       :: samdim1(:,:)                 !<  The positive-definite matrix.
    real(RKC)       , allocatable       :: samdim2(:,:)                 !<  The positive-definite matrix.
    type(bench_type), allocatable       :: bench(:)                     !<  The Benchmark array.
    integer(IK)     , parameter         :: nsammax = 2**NARR            !<  The maximum number of calculation repeats.
    integer(IK)     , parameter         :: ndim = 5_IK                  !<  The number of data attributes.
    real(RKC)                           :: mean(ndim), cov(ndim, ndim)  !<  The positive-definite matrix.
    integer(IK)                         :: idim, jdim, isam, nsam       !<  The benchmarking array size.
    real(RKC)                           :: dumm

    bench = [ bench_type(name = SK_"intrinsicDIM1", exec = intrinsicDIM1, overhead = setOverhead) &
            , bench_type(name = SK_"intrinsicDIM2", exec = intrinsicDIM2, overhead = setOverhead) &
            , bench_type(name = SK_"setCovDIM1", exec = setCovDIM1, overhead = setOverhead) &
            , bench_type(name = SK_"setCovDIM2", exec = setCovDIM2, overhead = setOverhead) &
            ]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "sample covariance benchmarking..."
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "nsam", (bench(i)%name, i = 1, size(bench))

        dumm = 0._RKC
        loopOverMatrixSize: do iarr = 1, NARR - 1

            nsam = 2**iarr
            ntry = nsammax / nsam
            allocate(samdim1(nsam, ndim), samdim2(ndim, nsam))
            write(*,"(*(g0,:,' '))") "Benchmarking setCovDIM1() vs. setCovDIM2()", nsam, ntry

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

    subroutine setCovDIM1()
        block
            use pm_sampleCov, only: setCov
            do itry = 1, ntry
                call setSamDIM1()
                call setCov(cov, uppDia, samdim1, dim = 1_IK)
                dumm = dumm + cov(1,1)
            end do
        end block
    end subroutine

    subroutine setCovDIM2()
        block
            use pm_sampleCov, only: setCov
            do itry = 1, ntry
                call setSamDIM2()
                call setCov(cov, uppDia, samdim2, dim = 2_IK)
                dumm = dumm + cov(1,1)
            end do
        end block
    end subroutine

    subroutine intrinsicDIM1()
        do itry = 1, ntry
            call setSamDIM1()
            do jdim = 1, ndim
                do idim = 1, jdim
                    cov(idim, jdim) = dot_product(samdim1(1 : nsam, idim), samdim1(1 : nsam, jdim)) / size(samdim1, dim = 1)
                end do
            end do
            dumm = dumm + cov(1,1)
        end do
    end subroutine

    subroutine intrinsicDIM2()
        do itry = 1, ntry
            call setSamDIM2()
            do jdim = 1, ndim
                do idim = 1, jdim
                    cov(idim, jdim) = dot_product(samdim2(idim, 1 : nsam), samdim2(jdim, 1 : nsam)) / size(samdim2, dim = 2)
                end do
            end do
            dumm = dumm + cov(1,1)
        end do
    end subroutine

end program benchmark