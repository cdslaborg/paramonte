#define MatB_ENABLED 0
! Test the performance of `transpose()` vs. `setMatTrans()`.
program benchmark

    use iso_fortran_env, only: error_unit
    use pm_kind, only: IK, RKG => RK, RK, SK
    use pm_distUnif, only: setUnifRand
    use pm_arraySpace, only: getLogSpace
    use pm_arrayReplace, only: getReplaced
    use pm_arrayUnique, only: getUnique
    use pm_bench, only: bench_type
    use pm_val2str, only: getStr

    implicit none

    integer(IK)                     :: i                    !<  The procedure benchmark counter.
    integer(IK)                     :: fileUnit             !<  The output file unit for benchmark results.
    integer(IK)                     :: iblock               !<  The matrix rank and its counter.
    integer(IK)                     :: bsize                !<  The matrix rank and its counter.
    integer(IK)     , parameter     :: RANK = 1000_IK       !<  The matrix rank.
    real(RKG)                       :: dummySum = 0._RKG    !<  The dummy computation to prevent the compiler from doing aggressive optimizations.
    integer(IK)     , allocatable   :: BlockSize(:)         !<  The vector of block sizes.
    type(bench_type), allocatable   :: bench(:)             !<  The Benchmark array.
    real(RKG)       , allocatable   :: matA(:,:)            !<  The matrix.
#if MatB_ENABLED
    real(RKG)       , allocatable   :: matB(:,:)            !<  The matrix transposed.
    allocate(matB(RANK, RANK))
    call setUnifRand(matB)
#endif
    allocate(matA(RANK, RANK))
    call setUnifRand(matA)

    bench = [ bench_type(name = getReplaced(SK_"setMatTrans(matA(RANK,RANK))", SK_"RANK", getStr(RANK)), exec = setMatTrans, overhead = setOverhead) &
            , bench_type(name = getReplaced(  SK_"transpose(matA(RANK,RANK))", SK_"RANK", getStr(RANK)), exec = transpose, overhead = setOverhead) &
            ]

    BlockSize = getUnique(int(getLogSpace(log(1._RKG), log(real(2*RANK, RKG)), count = 50_IK), IK))

    write(*,"(*(g0,:,' '))")
    write(*,"(*(g0,:,' '))") "setMatTransBlock"
    write(*,"(*(g0,:,' '))")

    open(newunit = fileUnit, file = "main.out", status = "replace")

        write(fileUnit, "(*(g0,:,', '))") "BlockSize", (bench(i)%name, i = 1, size(bench))

        loopOverMatrixRank: do iblock = 1, size(BlockSize)

            bsize = BlockSize(iblock)
            write(*,"(*(g0,:,' '))") "Benchmarking with block size", bsize

            do i = 1, size(bench)
                bench(i)%timing = bench(i)%getTiming(minsec = 0.07_RK)
            end do

            write(fileUnit,"(*(g0,:,', '))") bsize, (bench(i)%timing%mean, i = 1, size(bench))

        end do loopOverMatrixRank

    close(fileUnit)

    write(*,"(*(g0,:,' '))") dummySum
    write(*,"(*(g0,:,' '))")

contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! procedure wrappers.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOverhead()
        call getDummy()
    end subroutine

    subroutine getDummy()
#if     MatB_ENABLED
        dummySum = dummySum + matB(1,1)
#else
        dummySum = dummySum + matA(1,1)
#endif
    end subroutine

    subroutine setMatTrans()
        block
            use pm_matrixTrans, only: setMatTrans
#if         MatB_ENABLED
            call setMatTrans(matA, matB, bsize)
#else
            call setMatTrans(matA, bsize)
#endif
            call getDummy()
        end block
    end subroutine

    subroutine transpose()
        block
            intrinsic :: transpose
#if         MatB_ENABLED
            matB = transpose(matA)
#else
            matA = transpose(matA)
#endif
            call getDummy()
        end block
    end subroutine

end program benchmark