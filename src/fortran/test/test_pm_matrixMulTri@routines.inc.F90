!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief
!>  This include file contains procedure implementations of the tests of [pm_matrixMulTri](@ref pm_matrixMulTri).
!>
!>  \final
!>
!>  \author
!>  \FatemehBagheri, 12:27 AM Tuesday, February 22, 2022, Dallas, TX

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if     CK_ENABLED
#define TYPE_KIND complex(TKG)
        complex(TKG), parameter :: lb = (1._TKG, -2._TKG), ub = (2._TKG, -1._TKG), ONE = (1._TKG, 0._TKG), ZERO = (0._TKG, 0._TKG)
#elif   RK_ENABLED
#define TYPE_KIND real(TKG)
        real(TKG), parameter :: lb = 1._TKG, ub = 2._TKG, ONE = 1._TKG, ZERO = 0._TKG
#else
#error  "Unrecognized interface."
#endif
        logical(LK) :: simplified
        integer(IK) :: itry, begG, incG, endG, irow, nrow, ncol, ndim, roffT, coffT, roffG, coffG
        integer(IK), parameter :: increments(*) = [(irow, irow = -3, -1), (irow, irow = 1, 3)]
        TYPE_KIND, allocatable :: solmat(:,:), trimat(:,:), genmat(:,:), genvec(:), dumvec(:), choices(:)
        real(TKG), parameter :: TOL = epsilon(0._TKG)**.5
        TYPE_KIND :: alpha, alinv

        assertion = .true.

        do itry = 1, 100

            ! Build the matrices.

            incG = getChoice(increments)
            nrow = getUnifRand(0_IK, 7_IK)
            ncol = getUnifRand(0_IK, 7_IK)
            roffT = getUnifRand(0_IK, 3_IK)
            coffT = getUnifRand(0_IK, 3_IK)
            roffG = getUnifRand(0_IK, 3_IK)
            coffG = getUnifRand(0_IK, 3_IK)
            choices = [ZERO, ONE, getUnifRand(lb, ub)]
            alpha = getChoice(choices)
            if (alpha == ZERO) then
                alinv = alpha
            else
                alinv = 1 / alpha
            end if

            ndim = nrow
            genvec = getChoice([-1, 1]) * getUnifRand(lb, ub, max(0, 1 + (nrow - 1) * abs(incG)))
            genmat = getChoice([-1, 1]) * getUnifRand(lb, ub, 2 * roffG + nrow, 2 * coffG + ncol)
            trimat = getChoice([-1, 1]) * getUnifRand(lb, ub, 2 * roffT + ndim, 2 * coffT + ndim)
            if (incG < 0) then
                begG = size(genvec)
                endG = 1
            else
                begG = 1
                endG = size(genvec)
            end if

            ! BLAS - LEVEL 2: ?TRMV / ?TRSV - simplified interface

            simplified = .true.

            ! upperDiag, nothing/inversion
            dumvec = genvec
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), upperDiag, nothing, dumvec(begG:endG:incG))
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), upperDiag, inversion, dumvec(begG:endG:incG))
            call reportTRXV(simplified, "nothing/inversion")

            ! upperDiag, transSymm/transOrth
            dumvec = genvec
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), upperDiag, transSymm, dumvec(begG:endG:incG))
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), upperDiag, transOrth, dumvec(begG:endG:incG))
            call reportTRXV(simplified, "transSymm/transOrth")

            ! upperDiag, transHerm/transUnit
            dumvec = genvec
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), upperDiag, transHerm, dumvec(begG:endG:incG))
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), upperDiag, transUnit, dumvec(begG:endG:incG))
            call reportTRXV(simplified, "transHerm/transUnit")

            ! lowerDiag, nothing/inversion
            dumvec = genvec
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), lowerDiag, nothing, dumvec(begG:endG:incG))
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), lowerDiag, inversion, dumvec(begG:endG:incG))
            call reportTRXV(simplified, "nothing/inversion")

            ! lowerDiag, transSymm/transOrth
            dumvec = genvec
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), lowerDiag, transSymm, dumvec(begG:endG:incG))
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), lowerDiag, transOrth, dumvec(begG:endG:incG))
            call reportTRXV(simplified, "transSymm/transOrth")

            ! lowerDiag, transHerm/transUnit
            dumvec = genvec
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), lowerDiag, transHerm, dumvec(begG:endG:incG))
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), lowerDiag, transUnit, dumvec(begG:endG:incG))
            call reportTRXV(simplified, "transHerm/transUnit")

            ! upperUnit, nothing/inversion
            dumvec = genvec
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), upperUnit, nothing, dumvec(begG:endG:incG))
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), upperUnit, inversion, dumvec(begG:endG:incG))
            call reportTRXV(simplified, "nothing/inversion")

            ! upperUnit, transSymm/transOrth
            dumvec = genvec
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), upperUnit, transSymm, dumvec(begG:endG:incG))
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), upperUnit, transOrth, dumvec(begG:endG:incG))
            call reportTRXV(simplified, "transSymm/transOrth")

            ! upperUnit, transHerm/transUnit
            dumvec = genvec
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), upperUnit, transHerm, dumvec(begG:endG:incG))
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), upperUnit, transUnit, dumvec(begG:endG:incG))
            call reportTRXV(simplified, "transHerm/transUnit")

            ! lowerUnit, nothing/inversion
            dumvec = genvec
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), lowerUnit, nothing, dumvec(begG:endG:incG))
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), lowerUnit, inversion, dumvec(begG:endG:incG))
            call reportTRXV(simplified, "nothing/inversion")

            ! lowerUnit, transSymm/transOrth
            dumvec = genvec
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), lowerUnit, transSymm, dumvec(begG:endG:incG))
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), lowerUnit, transOrth, dumvec(begG:endG:incG))
            call reportTRXV(simplified, "transSymm/transOrth")

            ! lowerUnit, transHerm/transUnit
            dumvec = genvec
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), lowerUnit, transHerm, dumvec(begG:endG:incG))
            call setMatMulTri(trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim), lowerUnit, transUnit, dumvec(begG:endG:incG))
            call reportTRXV(simplified, "transHerm/transUnit")

            ! BLAS - LEVEL 2: ?TRMV / ?TRSV - contiguous interface

            simplified = .false.

            ! upperDiag, nothing/inversion
            dumvec = genvec
            call setMatMulTri(trimat, upperDiag, nothing, dumvec, ndim, roffT, coffT, incG)
            call setMatMulTri(trimat, upperDiag, inversion, dumvec, ndim, roffT, coffT, incG)
            call reportTRXV(simplified, "nothing/inversion")

            ! upperDiag, transSymm/transOrth
            dumvec = genvec
            call setMatMulTri(trimat, upperDiag, transSymm, dumvec, ndim, roffT, coffT, incG)
            call setMatMulTri(trimat, upperDiag, transOrth, dumvec, ndim, roffT, coffT, incG)
            call reportTRXV(simplified, "transSymm/transOrth")

            ! upperDiag, transHerm/transUnit
            dumvec = genvec
            call setMatMulTri(trimat, upperDiag, transHerm, dumvec, ndim, roffT, coffT, incG)
            call setMatMulTri(trimat, upperDiag, transUnit, dumvec, ndim, roffT, coffT, incG)
            call reportTRXV(simplified, "transHerm/transUnit")

            ! lowerDiag, nothing/inversion
            dumvec = genvec
            call setMatMulTri(trimat, lowerDiag, nothing, dumvec, ndim, roffT, coffT, incG)
            call setMatMulTri(trimat, lowerDiag, inversion, dumvec, ndim, roffT, coffT, incG)
            call reportTRXV(simplified, "nothing/inversion")

            ! lowerDiag, transSymm/transOrth
            dumvec = genvec
            call setMatMulTri(trimat, lowerDiag, transSymm, dumvec, ndim, roffT, coffT, incG)
            call setMatMulTri(trimat, lowerDiag, transOrth, dumvec, ndim, roffT, coffT, incG)
            call reportTRXV(simplified, "transSymm/transOrth")

            ! lowerDiag, transHerm/transUnit
            dumvec = genvec
            call setMatMulTri(trimat, lowerDiag, transHerm, dumvec, ndim, roffT, coffT, incG)
            call setMatMulTri(trimat, lowerDiag, transUnit, dumvec, ndim, roffT, coffT, incG)
            call reportTRXV(simplified, "transHerm/transUnit")

            ! upperUnit, nothing/inversion
            dumvec = genvec
            call setMatMulTri(trimat, upperUnit, nothing, dumvec, ndim, roffT, coffT, incG)
            call setMatMulTri(trimat, upperUnit, inversion, dumvec, ndim, roffT, coffT, incG)
            call reportTRXV(simplified, "nothing/inversion")

            ! upperUnit, transSymm/transOrth
            dumvec = genvec
            call setMatMulTri(trimat, upperUnit, transSymm, dumvec, ndim, roffT, coffT, incG)
            call setMatMulTri(trimat, upperUnit, transOrth, dumvec, ndim, roffT, coffT, incG)
            call reportTRXV(simplified, "transSymm/transOrth")

            ! upperUnit, transHerm/transUnit
            dumvec = genvec
            call setMatMulTri(trimat, upperUnit, transHerm, dumvec, ndim, roffT, coffT, incG)
            call setMatMulTri(trimat, upperUnit, transUnit, dumvec, ndim, roffT, coffT, incG)
            call reportTRXV(simplified, "transHerm/transUnit")

            ! lowerUnit, nothing/inversion
            dumvec = genvec
            call setMatMulTri(trimat, lowerUnit, nothing, dumvec, ndim, roffT, coffT, incG)
            call setMatMulTri(trimat, lowerUnit, inversion, dumvec, ndim, roffT, coffT, incG)
            call reportTRXV(simplified, "nothing/inversion")

            ! lowerUnit, transSymm/transOrth
            dumvec = genvec
            call setMatMulTri(trimat, lowerUnit, transSymm, dumvec, ndim, roffT, coffT, incG)
            call setMatMulTri(trimat, lowerUnit, transOrth, dumvec, ndim, roffT, coffT, incG)
            call reportTRXV(simplified, "transSymm/transOrth")

            ! lowerUnit, transHerm/transUnit
            dumvec = genvec
            call setMatMulTri(trimat, lowerUnit, transHerm, dumvec, ndim, roffT, coffT, incG)
            call setMatMulTri(trimat, lowerUnit, transUnit, dumvec, ndim, roffT, coffT, incG)
            call reportTRXV(simplified, "transHerm/transUnit")

            ! BLAS - LEVEL 3: ?TRMM / ?TRSM

            call testTRXM()
            call testTRXM(alpha, alinv)

        end do

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reportTRXV(simplified, operationA)
            character(*), intent(in) :: operationA
            logical(LK) , intent(in) :: simplified
            logical(LK), allocatable :: tolerable(:)
            type(display_type) :: disp
            tolerable = isClose(genvec, dumvec, abstol = TOL)
            assertion = assertion .and. all(tolerable)
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call disp%show("TRMV/TRSV")
                call disp%show("simplified")
                call disp%show( simplified )
                call disp%show("operationA")
                call disp%show( operationA )
                call disp%show("shape(trimat)")
                call disp%show( shape(trimat) )
                call disp%show("shape(genvec)")
                call disp%show( shape(genvec) )
                call disp%show("shape(dumvec)")
                call disp%show( shape(dumvec) )
                call disp%show("[nrow, ndim]")
                call disp%show( [nrow, ndim] )
                call disp%show("[roffT, coffT, begG, incG, endG]")
                call disp%show( [roffT, coffT, begG, incG, endG] )
                call disp%show("trimat(roffT + 1 : roffT + nrow, coffT + 1 : coffT + nrow)")
                call disp%show( trimat(roffT + 1 : roffT + nrow, coffT + 1 : coffT + nrow) )
                call disp%show("genvec(begG:endG:incG)")
                call disp%show( genvec(begG:endG:incG) )
                call disp%show("dumvec(begG:endG:incG)")
                call disp%show( dumvec(begG:endG:incG) )
                call disp%show("tolerable")
                call disp%show( tolerable )
                call disp%show("pack(dumvec, .not. tolerable)")
                call disp%show( pack(dumvec, .not. tolerable) )
                call disp%show("pack(genvec, .not. tolerable)")
                call disp%show( pack(genvec, .not. tolerable) )
                call disp%show("pack(abs(genvec - dumvec), .not. tolerable)")
                call disp%show( pack(abs(genvec - dumvec), .not. tolerable) )
                call disp%show("TOL")
                call disp%show( TOL )
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"TRMV/TRSV test with the above specifications must successfully pass.")
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine testTRXM(alpha, alinv)
            TYPE_KIND, intent(in), optional :: alpha, alinv
            character(:, SK), allocatable :: subset
            if (present(alpha) .neqv. present(alinv)) error stop "`alpha` and `alinv` must be both present or both missing." ! LCOV_EXCL_LINE
#define     sliceT trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim)
#define     sliceS solmat(roffG + 1 : roffG + nrow, coffG + 1 : coffG + ncol)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ndim = nrow
            trimat = getChoice([-1, 1]) * getUnifRand(lb, ub, 2 * roffT + ndim, 2 * coffT + ndim)
            genmat = getChoice([-1, 1]) * getUnifRand(lb, ub, 2 * roffG + nrow, 2 * coffG + ncol)

            subset = "upperDiag"

            ! upperDiag, nothing/inversion
            solmat = genmat
            call setMatMulTri(sliceT, upperDiag, nothing, sliceS, alpha)
            call setMatMulTri(sliceT, upperDiag, inversion, sliceS, alinv)
            call reportTRXM(simplified, "nothing/inversion", subset, alpha)

            ! upperDiag, transSymm/transOrth
            solmat = genmat
            call setMatMulTri(sliceT, upperDiag, transSymm, sliceS, alpha)
            call setMatMulTri(sliceT, upperDiag, transOrth, sliceS, alinv)
            call reportTRXM(simplified, "transSymm/transOrth", subset, alpha)

            ! upperDiag, transHerm/transUnit
            solmat = genmat
            call setMatMulTri(sliceT, upperDiag, transHerm, sliceS, alpha)
            call setMatMulTri(sliceT, upperDiag, transUnit, sliceS, alinv)
            call reportTRXM(simplified, "transHerm/transUnit", subset, alpha)

            subset = "upperDiag"

            ! lowerDiag, nothing/inversion
            solmat = genmat
            call setMatMulTri(sliceT, lowerDiag, nothing, sliceS, alpha)
            call setMatMulTri(sliceT, lowerDiag, inversion, sliceS, alinv)
            call reportTRXM(simplified, "nothing/inversion", subset, alpha)

            ! lowerDiag, transSymm/transOrth
            solmat = genmat
            call setMatMulTri(sliceT, lowerDiag, transSymm, sliceS, alpha)
            call setMatMulTri(sliceT, lowerDiag, transOrth, sliceS, alinv)
            call reportTRXM(simplified, "transSymm/transOrth", subset, alpha)

            ! lowerDiag, transHerm/transUnit
            solmat = genmat
            call setMatMulTri(sliceT, lowerDiag, transHerm, sliceS, alpha)
            call setMatMulTri(sliceT, lowerDiag, transUnit, sliceS, alinv)
            call reportTRXM(simplified, "transHerm/transUnit", subset, alpha)

            subset = "upperUnit"

            ! upperUnit, nothing/inversion
            solmat = genmat
            call setMatMulTri(sliceT, upperUnit, nothing, sliceS, alpha)
            call setMatMulTri(sliceT, upperUnit, inversion, sliceS, alinv)
            call reportTRXM(simplified, "nothing/inversion", subset, alpha)

            ! upperUnit, transSymm/transOrth
            solmat = genmat
            call setMatMulTri(sliceT, upperUnit, transSymm, sliceS, alpha)
            call setMatMulTri(sliceT, upperUnit, transOrth, sliceS, alinv)
            call reportTRXM(simplified, "transSymm/transOrth", subset, alpha)

            ! upperUnit, transHerm/transUnit
            solmat = genmat
            call setMatMulTri(sliceT, upperUnit, transHerm, sliceS, alpha)
            call setMatMulTri(sliceT, upperUnit, transUnit, sliceS, alinv)
            call reportTRXM(simplified, "transHerm/transUnit", subset, alpha)

            subset = "upperUnit"

            ! lowerUnit, nothing/inversion
            solmat = genmat
            call setMatMulTri(sliceT, lowerUnit, nothing, sliceS, alpha)
            call setMatMulTri(sliceT, lowerUnit, inversion, sliceS, alinv)
            call reportTRXM(simplified, "nothing/inversion", subset, alpha)

            ! lowerUnit, transSymm/transOrth
            solmat = genmat
            call setMatMulTri(sliceT, lowerUnit, transSymm, sliceS, alpha)
            call setMatMulTri(sliceT, lowerUnit, transOrth, sliceS, alinv)
            call reportTRXM(simplified, "transSymm/transOrth", subset, alpha)

            ! lowerUnit, transHerm/transUnit
            solmat = genmat
            call setMatMulTri(sliceT, lowerUnit, transHerm, sliceS, alpha)
            call setMatMulTri(sliceT, lowerUnit, transUnit, sliceS, alinv)
            call reportTRXM(simplified, "transHerm/transUnit", subset, alpha)

            if (present(alpha)) then

            subset = "upperDiag"

            ! upperDiag, nothing/inversion
            solmat = genmat
            call setMatMulTri(trimat, upperDiag, nothing, solmat, alpha, nrow, ncol, roffT, coffT, roffG, coffG)
            call setMatMulTri(trimat, upperDiag, inversion, solmat, alinv, nrow, ncol, roffT, coffT, roffG, coffG)
            call reportTRXM(simplified, "nothing/inversion", subset, alpha)

            ! upperDiag, transSymm/transOrth
            solmat = genmat
            call setMatMulTri(trimat, upperDiag, transSymm, solmat, alpha, nrow, ncol, roffT, coffT, roffG, coffG)
            call setMatMulTri(trimat, upperDiag, transOrth, solmat, alinv, nrow, ncol, roffT, coffT, roffG, coffG)
            call reportTRXM(simplified, "transSymm/transOrth", subset, alpha)

            ! upperDiag, transHerm/transUnit
            solmat = genmat
            call setMatMulTri(trimat, upperDiag, transHerm, solmat, alpha, nrow, ncol, roffT, coffT, roffG, coffG)
            call setMatMulTri(trimat, upperDiag, transUnit, solmat, alinv, nrow, ncol, roffT, coffT, roffG, coffG)
            call reportTRXM(simplified, "transHerm/transUnit", subset, alpha)

            subset = "upperDiag"

            ! lowerDiag, nothing/inversion
            solmat = genmat
            call setMatMulTri(trimat, lowerDiag, nothing, solmat, alpha, nrow, ncol, roffT, coffT, roffG, coffG)
            call setMatMulTri(trimat, lowerDiag, inversion, solmat, alinv, nrow, ncol, roffT, coffT, roffG, coffG)
            call reportTRXM(simplified, "nothing/inversion", subset, alpha)

            ! lowerDiag, transSymm/transOrth
            solmat = genmat
            call setMatMulTri(trimat, lowerDiag, transSymm, solmat, alpha, nrow, ncol, roffT, coffT, roffG, coffG)
            call setMatMulTri(trimat, lowerDiag, transOrth, solmat, alinv, nrow, ncol, roffT, coffT, roffG, coffG)
            call reportTRXM(simplified, "transSymm/transOrth", subset, alpha)

            ! lowerDiag, transHerm/transUnit
            solmat = genmat
            call setMatMulTri(trimat, lowerDiag, transHerm, solmat, alpha, nrow, ncol, roffT, coffT, roffG, coffG)
            call setMatMulTri(trimat, lowerDiag, transUnit, solmat, alinv, nrow, ncol, roffT, coffT, roffG, coffG)
            call reportTRXM(simplified, "transHerm/transUnit", subset, alpha)

            subset = "upperUnit"

            ! upperUnit, nothing/inversion
            solmat = genmat
            call setMatMulTri(trimat, upperUnit, nothing, solmat, alpha, nrow, ncol, roffT, coffT, roffG, coffG)
            call setMatMulTri(trimat, upperUnit, inversion, solmat, alinv, nrow, ncol, roffT, coffT, roffG, coffG)
            call reportTRXM(simplified, "nothing/inversion", subset, alpha)

            ! upperUnit, transSymm/transOrth
            solmat = genmat
            call setMatMulTri(trimat, upperUnit, transSymm, solmat, alpha, nrow, ncol, roffT, coffT, roffG, coffG)
            call setMatMulTri(trimat, upperUnit, transOrth, solmat, alinv, nrow, ncol, roffT, coffT, roffG, coffG)
            call reportTRXM(simplified, "transSymm/transOrth", subset, alpha)

            ! upperUnit, transHerm/transUnit
            solmat = genmat
            call setMatMulTri(trimat, upperUnit, transHerm, solmat, alpha, nrow, ncol, roffT, coffT, roffG, coffG)
            call setMatMulTri(trimat, upperUnit, transUnit, solmat, alinv, nrow, ncol, roffT, coffT, roffG, coffG)
            call reportTRXM(simplified, "transHerm/transUnit", subset, alpha)

            subset = "upperUnit"

            ! lowerUnit, nothing/inversion
            solmat = genmat
            call setMatMulTri(trimat, lowerUnit, nothing, solmat, alpha, nrow, ncol, roffT, coffT, roffG, coffG)
            call setMatMulTri(trimat, lowerUnit, inversion, solmat, alinv, nrow, ncol, roffT, coffT, roffG, coffG)
            call reportTRXM(simplified, "nothing/inversion", subset, alpha)

            ! lowerUnit, transSymm/transOrth
            solmat = genmat
            call setMatMulTri(trimat, lowerUnit, transSymm, solmat, alpha, nrow, ncol, roffT, coffT, roffG, coffG)
            call setMatMulTri(trimat, lowerUnit, transOrth, solmat, alinv, nrow, ncol, roffT, coffT, roffG, coffG)
            call reportTRXM(simplified, "transSymm/transOrth", subset, alpha)

            ! lowerUnit, transHerm/transUnit
            solmat = genmat
            call setMatMulTri(trimat, lowerUnit, transHerm, solmat, alpha, nrow, ncol, roffT, coffT, roffG, coffG)
            call setMatMulTri(trimat, lowerUnit, transUnit, solmat, alinv, nrow, ncol, roffT, coffT, roffG, coffG)
            call reportTRXM(simplified, "transHerm/transUnit", subset, alpha)

            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ndim = ncol
            genmat = getChoice([-1, 1]) * getUnifRand(lb, ub, 2 * roffG + nrow, 2 * coffG + ncol)
            trimat = getChoice([-1, 1]) * getUnifRand(lb, ub, 2 * roffT + ndim, 2 * coffT + ndim)
            solmat = genmat

            subset = "upperDiag"

            ! upperDiag, nothing/inversion
            solmat = genmat
            call setMatMulTri(sliceS, sliceT, upperDiag, nothing, alpha)
            call setMatMulTri(sliceS, sliceT, upperDiag, inversion, alinv)
            call reportTRXM(simplified, "nothing/inversion", subset, alpha)

            ! upperDiag, transSymm/transOrth
            solmat = genmat
            call setMatMulTri(sliceS, sliceT, upperDiag, transSymm, alpha)
            call setMatMulTri(sliceS, sliceT, upperDiag, transOrth, alinv)
            call reportTRXM(simplified, "transSymm/transOrth", subset, alpha)

            ! upperDiag, transHerm/transUnit
            solmat = genmat
            call setMatMulTri(sliceS, sliceT, upperDiag, transHerm, alpha)
            call setMatMulTri(sliceS, sliceT, upperDiag, transUnit, alinv)
            call reportTRXM(simplified, "transHerm/transUnit", subset, alpha)

            subset = "lowerDiag"

            ! lowerDiag, nothing/inversion
            solmat = genmat
            call setMatMulTri(sliceS, sliceT, lowerDiag, nothing, alpha)
            call setMatMulTri(sliceS, sliceT, lowerDiag, inversion, alinv)
            call reportTRXM(simplified, "nothing/inversion", subset, alpha)

            ! lowerDiag, transSymm/transOrth
            solmat = genmat
            call setMatMulTri(sliceS, sliceT, lowerDiag, transSymm, alpha)
            call setMatMulTri(sliceS, sliceT, lowerDiag, transOrth, alinv)
            call reportTRXM(simplified, "transSymm/transOrth", subset, alpha)

            ! lowerDiag, transHerm/transUnit
            solmat = genmat
            call setMatMulTri(sliceS, sliceT, lowerDiag, transHerm, alpha)
            call setMatMulTri(sliceS, sliceT, lowerDiag, transUnit, alinv)
            call reportTRXM(simplified, "transHerm/transUnit", subset, alpha)

            subset = "upperUnit"

            ! upperUnit, nothing/inversion
            solmat = genmat
            call setMatMulTri(sliceS, sliceT, upperUnit, nothing, alpha)
            call setMatMulTri(sliceS, sliceT, upperUnit, inversion, alinv)
            call reportTRXM(simplified, "nothing/inversion", subset, alpha)

            ! upperUnit, transSymm/transOrth
            solmat = genmat
            call setMatMulTri(sliceS, sliceT, upperUnit, transSymm, alpha)
            call setMatMulTri(sliceS, sliceT, upperUnit, transOrth, alinv)
            call reportTRXM(simplified, "transSymm/transOrth", subset, alpha)

            ! upperUnit, transHerm/transUnit
            solmat = genmat
            call setMatMulTri(sliceS, sliceT, upperUnit, transHerm, alpha)
            call setMatMulTri(sliceS, sliceT, upperUnit, transUnit, alinv)
            call reportTRXM(simplified, "transHerm/transUnit", subset, alpha)

            subset = "lowerUnit"

            ! lowerUnit, nothing/inversion
            solmat = genmat
            call setMatMulTri(sliceS, sliceT, lowerUnit, nothing, alpha)
            call setMatMulTri(sliceS, sliceT, lowerUnit, inversion, alinv)
            call reportTRXM(simplified, "nothing/inversion", subset, alpha)

            ! lowerUnit, transSymm/transOrth
            solmat = genmat
            call setMatMulTri(sliceS, sliceT, lowerUnit, transSymm, alpha)
            call setMatMulTri(sliceS, sliceT, lowerUnit, transOrth, alinv)
            call reportTRXM(simplified, "transSymm/transOrth", subset, alpha)

            ! lowerUnit, transHerm/transUnit
            solmat = genmat
            call setMatMulTri(sliceS, sliceT, lowerUnit, transHerm, alpha)
            call setMatMulTri(sliceS, sliceT, lowerUnit, transUnit, alinv)
            call reportTRXM(simplified, "transHerm/transUnit", subset, alpha)

            if (present(alpha)) then

            subset = "upperDiag"

            ! upperDiag, nothing/inversion
            solmat = genmat
            call setMatMulTri(solmat, trimat, upperDiag, nothing, alpha, nrow, ncol, roffG, coffG, roffT, coffT)
            call setMatMulTri(solmat, trimat, upperDiag, inversion, alinv, nrow, ncol, roffG, coffG, roffT, coffT)
            call reportTRXM(simplified, "nothing/inversion", subset, alpha)

            ! upperDiag, transSymm/transOrth
            solmat = genmat
            call setMatMulTri(solmat, trimat, upperDiag, transSymm, alpha, nrow, ncol, roffG, coffG, roffT, coffT)
            call setMatMulTri(solmat, trimat, upperDiag, transOrth, alinv, nrow, ncol, roffG, coffG, roffT, coffT)
            call reportTRXM(simplified, "transSymm/transOrth", subset, alpha)

            ! upperDiag, transHerm/transUnit
            solmat = genmat
            call setMatMulTri(solmat, trimat, upperDiag, transHerm, alpha, nrow, ncol, roffG, coffG, roffT, coffT)
            call setMatMulTri(solmat, trimat, upperDiag, transUnit, alinv, nrow, ncol, roffG, coffG, roffT, coffT)
            call reportTRXM(simplified, "transHerm/transUnit", subset, alpha)

            subset = "lowerDiag"

            ! lowerDiag, nothing/inversion
            solmat = genmat
            call setMatMulTri(solmat, trimat, lowerDiag, nothing, alpha, nrow, ncol, roffG, coffG, roffT, coffT)
            call setMatMulTri(solmat, trimat, lowerDiag, inversion, alinv, nrow, ncol, roffG, coffG, roffT, coffT)
            call reportTRXM(simplified, "nothing/inversion", subset, alpha)

            ! lowerDiag, transSymm/transOrth
            solmat = genmat
            call setMatMulTri(solmat, trimat, lowerDiag, transSymm, alpha, nrow, ncol, roffG, coffG, roffT, coffT)
            call setMatMulTri(solmat, trimat, lowerDiag, transOrth, alinv, nrow, ncol, roffG, coffG, roffT, coffT)
            call reportTRXM(simplified, "transSymm/transOrth", subset, alpha)

            ! lowerDiag, transHerm/transUnit
            solmat = genmat
            call setMatMulTri(solmat, trimat, lowerDiag, transHerm, alpha, nrow, ncol, roffG, coffG, roffT, coffT)
            call setMatMulTri(solmat, trimat, lowerDiag, transUnit, alinv, nrow, ncol, roffG, coffG, roffT, coffT)
            call reportTRXM(simplified, "transHerm/transUnit", subset, alpha)

            subset = "upperUnit"

            ! upperUnit, nothing/inversion
            solmat = genmat
            call setMatMulTri(solmat, trimat, upperUnit, nothing, alpha, nrow, ncol, roffG, coffG, roffT, coffT)
            call setMatMulTri(solmat, trimat, upperUnit, inversion, alinv, nrow, ncol, roffG, coffG, roffT, coffT)
            call reportTRXM(simplified, "nothing/inversion", subset, alpha)

            ! upperUnit, transSymm/transOrth
            solmat = genmat
            call setMatMulTri(solmat, trimat, upperUnit, transSymm, alpha, nrow, ncol, roffG, coffG, roffT, coffT)
            call setMatMulTri(solmat, trimat, upperUnit, transOrth, alinv, nrow, ncol, roffG, coffG, roffT, coffT)
            call reportTRXM(simplified, "transSymm/transOrth", subset, alpha)

            ! upperUnit, transHerm/transUnit
            solmat = genmat
            call setMatMulTri(solmat, trimat, upperUnit, transHerm, alpha, nrow, ncol, roffG, coffG, roffT, coffT)
            call setMatMulTri(solmat, trimat, upperUnit, transUnit, alinv, nrow, ncol, roffG, coffG, roffT, coffT)
            call reportTRXM(simplified, "transHerm/transUnit", subset, alpha)

            subset = "lowerUnit"

            ! lowerUnit, nothing/inversion
            solmat = genmat
            call setMatMulTri(solmat, trimat, lowerUnit, nothing, alpha, nrow, ncol, roffG, coffG, roffT, coffT)
            call setMatMulTri(solmat, trimat, lowerUnit, inversion, alinv, nrow, ncol, roffG, coffG, roffT, coffT)
            call reportTRXM(simplified, "nothing/inversion", subset, alpha)

            ! lowerUnit, transSymm/transOrth
            solmat = genmat
            call setMatMulTri(solmat, trimat, lowerUnit, transSymm, alpha, nrow, ncol, roffG, coffG, roffT, coffT)
            call setMatMulTri(solmat, trimat, lowerUnit, transOrth, alinv, nrow, ncol, roffG, coffG, roffT, coffT)
            call reportTRXM(simplified, "transSymm/transOrth", subset, alpha)

            ! lowerUnit, transHerm/transUnit
            solmat = genmat
            call setMatMulTri(solmat, trimat, lowerUnit, transHerm, alpha, nrow, ncol, roffG, coffG, roffT, coffT)
            call setMatMulTri(solmat, trimat, lowerUnit, transUnit, alinv, nrow, ncol, roffG, coffG, roffT, coffT)
            call reportTRXM(simplified, "transHerm/transUnit", subset, alpha)

            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef      sliceT
#undef      sliceS
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subroutine reportTRXM(simplified, operation, subset, alpha)
            logical(LK) , intent(in) :: simplified
            character(*), intent(in) :: operation, subset
            TYPE_KIND, intent(in), optional :: alpha
            logical(LK), allocatable :: tolerable(:,:)
            type(display_type) :: disp
            tolerable = isClose(genmat, solmat, abstol = TOL)
            assertion = assertion .and. (all(tolerable) .or. (all(solmat(roffG + 1 : roffG + nrow, coffG + 1 : coffG + ncol) == ZERO) .and. getOption(ONE, alpha) == ZERO))
            if (test%traceable .and. .not. assertion) then
                ! LCOV_EXCL_START
                call disp%show("TRMM/TRSM")
                call disp%show("simplified")
                call disp%show( simplified )
                call disp%show("subset")
                call disp%show( subset )
                call disp%show("operation")
                call disp%show( operation )
                call disp%show("shape(trimat)")
                call disp%show( shape(trimat) )
                call disp%show("shape(genmat)")
                call disp%show( shape(genmat) )
                call disp%show("shape(solmat)")
                call disp%show( shape(solmat) )
                call disp%show("present(alpha)")
                call disp%show( present(alpha) )
                call disp%show("getOption(ONE, alpha)")
                call disp%show( getOption(ONE, alpha) )
                call disp%show("[nrow, ncol, ndim]")
                call disp%show( [nrow, ncol, ndim] )
                call disp%show("[roffT, coffT, roffG, coffG]")
                call disp%show( [roffT, coffT, roffG, coffG] )
                call disp%show("trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim)")
                call disp%show( trimat(roffT + 1 : roffT + ndim, coffT + 1 : coffT + ndim) )
                call disp%show("genmat(roffG + 1 : roffG + nrow, coffG + 1 : coffG + ncol)")
                call disp%show( genmat(roffG + 1 : roffG + nrow, coffG + 1 : coffG + ncol) )
                call disp%show("solmat(roffG + 1 : roffG + nrow, coffG + 1 : coffG + ncol)")
                call disp%show( solmat(roffG + 1 : roffG + nrow, coffG + 1 : coffG + ncol) )
                call disp%show("tolerable")
                call disp%show( tolerable )
                call disp%show("pack(solmat, .not. tolerable)")
                call disp%show( pack(solmat, .not. tolerable) )
                call disp%show("pack(genmat, .not. tolerable)")
                call disp%show( pack(genmat, .not. tolerable) )
                call disp%show("pack(abs(genmat - solmat), .not. tolerable)")
                call disp%show( pack(abs(genmat - solmat), .not. tolerable) )
                call disp%show("TOL")
                call disp%show( TOL )
                ! LCOV_EXCL_STOP
            end if
            call test%assert(assertion, SK_"TRMM/TRSM test with the above specifications must successfully pass.")
        end subroutine

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef  TYPE_KIND