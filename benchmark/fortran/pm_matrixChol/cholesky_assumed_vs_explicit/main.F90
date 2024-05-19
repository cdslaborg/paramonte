! Test the performance of Cholesky factorization computation using an assumed-shape interface vs. explicit-shape interface.
program benchmark

    use pm_kind, only: IK, LK, RKG => RKD, SK
    use pm_matrixChol, only: lowDia_type, uppDia_type
    use pm_matrixChol, only: subset_type => lowDia_type!uppDia_type
    use pm_bench, only: bench_type

    implicit none

    integer(IK)                         :: itry, ntry
    integer(IK)                         :: i                    !<  The procedure benchmark counter.
    integer(IK)                         :: iarr                 !<  The array size counter.
    integer(IK)                         :: fileUnit             !<  The output file unit for benchmark results.
    integer(IK)     , parameter         :: NARR = 11_IK         !<  The number of benchmark array sizes.
    real(RKG)       , allocatable       :: mat(:,:), choDia(:)  !<  The positive-definite mat.
    type(bench_type), allocatable       :: bench(:)             !<  The Benchmark array.
    integer(IK)     , parameter         :: nsim = 2**NARR       !<  The maximum number of calculation repeats.
    integer(IK)                         :: rank                 !<  The benchmarking array size.
    type(subset_type), parameter        :: subset = subset_type()
    integer(IK)                         :: offset
   !real(RKG)                           :: dumm
    offset = merge(1, 0, same_type_as(subset, lowDia_type()))

    bench = [ bench_type(name = SK_"setMatCholComplement", exec = setMatCholComplement, overhead = setOverhead) &
            , bench_type(name = SK_"setMatCholOverwrite", exec = setMatCholOverwrite, overhead = setOverhead) &
            , bench_type(name = SK_"unsafeExplicitShape", exec = unsafeExplicitShape, overhead = setOverhead) &
            , bench_type(name = SK_"setMatCholRecursive", exec = setMatCholRecursive, overhead = setOverhead) &
            , bench_type(name = SK_"setMatCholLooping", exec = setMatCholLooping, overhead = setOverhead) &
#if         LAPACK_ENABLED
            , bench_type(name = SK_"lapack_dpotrf", exec = lapack_dpotrf, overhead = setOverhead) &
#endif
            ]

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "assumed-shape vs. explicit-shape setChoLow()."
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,','))") "rank", (bench(i)%name, i = 1, size(bench))

        !dumm = 0._RKG
        loopOverMatrixSize: do iarr = 1, NARR

            rank = 2**iarr
            ntry = nsim / rank
            allocate(mat(rank, 0 : rank), choDia(rank))
            write(*,"(*(g0,:,' '))") "Benchmarking setChoLow() algorithms with array size", rank, ntry

            do i = 1, size(bench)
                bench(i)%timing = bench(i)%getTiming()
            end do

            write(fileUnit,"(*(g0,:,','))") rank, (bench(i)%timing%mean / ntry, i = 1, size(bench))
            deallocate(mat, choDia)

        end do loopOverMatrixSize
        write(*,"(*(g0,:,' '))")

    close(fileUnit)

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        do itry = 1, ntry
            call getMatrix()
        end do
    end subroutine

    subroutine getMatrix()
        integer(IK) :: i
        call random_number(mat)
        mat = mat * 1.e-5_RKG
        do i = 1, size(mat, dim = 1, kind = IK)
            mat(i, i - offset) = 1._RKG ! lowDia
           !mat(i, i) = 1._RKG ! uppDia
        end do
    end subroutine

#if LAPACK_ENABLED
    subroutine lapack_dpotrf()
        integer(IK) :: info
        do itry = 1, ntry
            call getMatrix()
            call dpotrf("U", rank, mat(1,1), rank, info)
            if (info /= 0_IK) error stop
        end do
    end subroutine
#endif

    subroutine unsafeExplicitShape()
        use pm_matrixChol, only: setChoLow
        logical(LK) :: failed
        do itry = 1, ntry
            call getMatrix()
            call setChoLow(mat(:,1-offset:rank-offset), choDia, rank)
            if (choDia(1) < 0._RKG) error stop
        end do
    end subroutine

    subroutine setMatCholOverwrite()
        use pm_matrixChol, only: setMatChol, nothing
        integer(IK) :: info
        do itry = 1, ntry
            call getMatrix()
            call setMatChol(mat(:,1-offset:rank-offset), subset, info, mat(:,1-offset:rank-offset), nothing)
            if (info /= 0_IK) error stop
        end do
    end subroutine

    subroutine setMatCholComplement()
        use pm_matrixChol, only: setMatChol, transHerm
        integer(IK) :: info
        do itry = 1, ntry
            call getMatrix()
           !call setMatChol(mat(:,0:rank-1), subset, info, mat(:,1:rank), transHerm)
            call setMatChol(mat(:,1-offset:rank-offset), subset, info, mat(:,offset:rank), transHerm)
            if (info /= 0_IK) error stop
        end do
    end subroutine

    subroutine setMatCholLooping()
        use pm_matrixChol, only: setMatChol, iteration
        integer(IK) :: info
        do itry = 1, ntry
            call getMatrix()
            call setMatChol(mat(:,1-offset:rank-offset), subset, info, iteration)
            if (info /= 0_IK) error stop
        end do
    end subroutine

    subroutine setMatCholRecursive()
        use pm_matrixChol, only: setMatChol, recursion
        integer(IK) :: info
        do itry = 1, ntry
            call getMatrix()
            call setMatChol(mat(:,1-offset:rank-offset), subset, info, recursion)
            if (info /= 0_IK) error stop
        end do
    end subroutine

end program benchmark