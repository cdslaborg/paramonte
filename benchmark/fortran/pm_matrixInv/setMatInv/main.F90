! Test the performance of Cholesky factorization computation using an assumed-shape interface vs. explicit-shape interface.
program benchmark

    use pm_kind, only: IK, LK, RKC => RKD, SK
    use pm_matrixCopy, only: setMatCopy, rdpack, uppDia, lowDia, transHerm
    use pm_distUnif, only: rngx_type => xoshiro256ssw_type
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: itry, ntry
    integer(IK)                         :: i                    !<  The procedure benchmark counter.
    integer(IK)                         :: iarr                 !<  The array size counter.
    integer(IK)                         :: fileUnit             !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: NARR = 10_IK         !<  The number of benchmark array sizes.
    integer(IK)     , allocatable       :: rperm(:)             !<  The permutation vector for LUP factorization.
    real(RKC)       , allocatable       :: mat(:,:), inv(:,:)   !<  The positive-definite matrix.
    type(bench_type), allocatable       :: bench(:)             !<  The Benchmark array.
    integer(IK)     , parameter         :: nsim = 2**NARR       !<  The maximum number of calculation repeats.
    integer(IK)                         :: rank                 !<  The benchmarking array size.
    real(RKC)                           :: dumm
    type(rngx_type)                     :: rngx

    rngx = rngx_type()
    bench = [ bench_type(name = SK_"setMatInvLow", exec = setMatInvLow, overhead = setOverhead) &
            , bench_type(name = SK_"setMatInvUpp", exec = setMatInvUpp, overhead = setOverhead) &
            , bench_type(name = SK_"setMatInvLUP", exec = setMatInvLUP, overhead = setOverhead) &
#if         LAPACK_ENABLED
            , bench_type(name = SK_"lapack_dpotri", exec = lapack_dpotri, overhead = setOverhead) &
#endif
            ]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "inverse matrix benchmarking..."
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "rank", (bench(i)%name, i = 1, size(bench))

        dumm = 0._RKC
        loopOverMatrixSize: do iarr = 1, NARR

            rank = 2**iarr
            ntry = nsim / rank
            allocate(mat(rank, rank + 1), inv(rank, rank), rperm(rank))
            write(*,"(*(g0,:,' '))") "Benchmarking setMatInv() algorithms with array size", rank, ntry

            do i = 1, size(bench)
                bench(i)%timing = bench(i)%getTiming()
            end do

            write(fileUnit,"(*(g0,:,','))") rank, (bench(i)%timing%mean / ntry, i = 1, size(bench))
            deallocate(mat, inv, rperm)

        end do loopOverMatrixSize
        write(*,"(*(g0,:,' '))") dumm

    close(fileUnit)

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        do itry = 1, ntry
            call setMatrix()
        end do
    end subroutine

    subroutine setMatrix()
        !integer(IK) :: i
        !call random_number(mat)
        !mat = mat * 1.e-4_RKC
        !do i = 1, size(mat, dim = 1, kind = IK)
        !    mat(i,i+1) = 1._RKC
        !end do
        !use pm_distCov, only: setCovRand
        !call setCovRand(rngx, mat(:,2:rank+1)) ! causes numeric overflow for large matrix ranks.
        use pm_matrixInit, only: setMatInit, uppLowDia
        call setMatInit(mat(:,2:rank+1), uppLowDia, 0._RKC, 0._RKC, 1._RKC)
        dumm = dumm - mat(rank, rank) - mat(rank-1, rank-1)
    end subroutine

#if LAPACK_ENABLED
    subroutine lapack_dpotri()
        integer(IK) :: info
        do itry = 1, ntry
            call setMatrix()
            call dpotrf("U", rank, mat(:,2:rank+1), rank, info)
            if (info /= 0_IK) error stop
            call dpotri("U", rank, mat(:,2:rank+1), rank, info)
            if (info /= 0_IK) error stop
            call setMatCopy(mat(:,2:rank+1), rdpack, mat(:,1:rank), rdpack, uppDia, transHerm) ! symmetrize inverse matrix.
            !dumm = dumm + mat(rank,rank) + mat(rank-1,rank-1)
        end do
    end subroutine
#endif

    subroutine setMatInvLow()
        use pm_matrixInv, only: setMatInv, choUpp
        use pm_matrixChol, only: setMatChol, nothing
        integer(IK) :: info
        do itry = 1, ntry
            call setMatrix()
            call setMatChol(mat(:,2:rank+1), uppDia, info, mat(:,2:rank+1), nothing)
            if (info /= 0_IK) error stop
            call setMatInv(mat(:,1:rank), mat(:,2:rank+1), choUpp)
            !dumm = dumm + mat(rank,rank) + mat(rank-1,rank-1)
        end do
    end subroutine

    subroutine setMatInvUpp()
        use pm_matrixInv, only: setMatInv, choLow, choUpp
        use pm_matrixChol, only: setMatChol
        integer(IK) :: info
        do itry = 1, ntry
            call setMatrix()
            call setMatChol(mat(:,2:rank+1), uppDia, info, mat(:,1:rank), transHerm)
            if (info /= 0_IK) error stop
            call setMatInv(mat(:,2:rank+1), mat(:,1:rank), choLow)
            !dumm = dumm + mat(rank,rank) + mat(rank-1,rank-1)
        end do
    end subroutine

    subroutine setMatInvLUP()
        use pm_matrixLUP, only: setMatLUP
        use pm_matrixInv, only: setMatInv
        integer(IK), parameter :: offset = 1
        integer(IK) :: info
        do itry = 1, ntry
            call setMatrix()
            call setMatLUP(mat(:,2:rank+1), rperm, info)
            if (info /= 0_IK) error stop
            call setMatInv(inv, mat(:,2:rank+1), rperm)
            !dumm = dumm + mat(rank,rank) + mat(rank-1,rank-1)
        end do
    end subroutine

end program benchmark