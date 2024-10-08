! Test the performance of Cholesky factorization computation using an assumed-shape interface vs. explicit-shape interface.
program benchmark

    use pm_kind, only: IK, LK, RKG => RKD, SK
    use pm_distUnif, only: setUnifRand
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: icls
    integer(IK)                         :: ibench
    integer(IK)                         :: fileUnit                     !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: nsam = 1000                  !<  The number of kmeans clusters.
    integer(IK)     , parameter         :: ndim = 3_IK                  !<  The number of data attributes.
    integer(IK)     , parameter         :: ncls = 20_IK                 !<  The number of data attributes.
    integer(IK)                         :: membership(nsam)             !<  cluster memberships.
    real(RKG)                           :: sample(ndim, nsam)           !<  The sample.
    real(RKG)                           :: center(ndim, ncls)           !<  The cluster centers.
    real(RKG)                           :: disq(nsam)                   !<  The distances squared.
    type(bench_type), allocatable       :: bench(:)                     !<  The Benchmark array.
    real(RKG)                           :: mean(ndim)                   !<  The sample mean.
    integer(IK)                         :: idum = 0
    logical(LK)                         :: mchanged

    bench = [ bench_type(name = SK_"default", exec = default, overhead = setOverhead) &
            , bench_type(name = SK_"changed", exec = changed, overhead = setOverhead) &
            ]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "membership benchmarking..."
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "ClusterCount", (bench(ibench)%name, ibench = 1, size(bench))

        do icls = 1, 20

            write(*,"(*(g0,:,' '))") "Benchmarking default() vs. changed()", nsam

            call random_number(disq)
            call random_number(center)
            call random_number(sample)
            call setUnifRand(membership, 1, icls)
            do ibench = 1, size(bench)
                bench(ibench)%timing = bench(ibench)%getTiming()
            end do

            write(fileUnit,"(*(g0,:,','))") icls, (bench(ibench)%timing%mean, ibench = 1, size(bench))

        end do
        write(*,"(*(g0,:,' '))") idum

    close(fileUnit)

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        if (mchanged) idum = idum + 1
    end subroutine

    subroutine default()
        use pm_clusKmeans, only: setMember
        integer(IK) :: membersnew(size(membership))
        call setMember(membersnew, disq, sample, center(:, 1 : icls))
        mchanged = all(membersnew == membership)
    end subroutine

    subroutine changed()
        use pm_clusKmeans, only: setMember
        call setMember(membership, disq, sample, center(:, 1 : icls), mchanged)
    end subroutine

end program benchmark