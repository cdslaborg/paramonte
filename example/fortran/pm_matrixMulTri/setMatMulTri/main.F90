program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: CKC => CK32 ! all processor types and kinds are supported.
    use pm_kind, only: RKC => RK32 ! all processor types and kinds are supported.
    use pm_matrixMulTri, only: upperDiag, lowerDiag
    use pm_matrixMulTri, only: upperUnit, lowerUnit
    use pm_matrixMulTri, only: transSymm, transHerm
    use pm_matrixMulTri, only: transOrth, transUnit
    use pm_matrixMulTri, only: inversion, nothing
    use pm_matrixMulTri, only: setMatMulTri
    use pm_io, only: display_type
    use pm_io, only: getFormat

    implicit none

    type(display_type) :: disp
    character(:, SK), allocatable :: cform, gform
    integer(IK) :: nrow, ncol, ndim, roffA, coffA, roffB, coffB, incB
    cform = getFormat([cmplx(0., 0., CKC)], ed = SK_'f', signed = .true.)
    gform = getFormat([real(0., RKC)], ed = SK_'f', signed = .true.)

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! BLAS 3 TRSV / TRMV: triangular matrix-vector multiplication: complex.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        complex(CKC) :: alpha, beta
        complex(CKC), parameter :: COMPLEXDUM = cmplx(huge(0._CKC), huge(0._CKC), CKC)
        complex(CKC), allocatable :: triMat(:,:), genMat(:), genRef(:), solMat(:)

        genRef = [complex(CKC) :: COMPLEXDUM, COMPLEXDUM, (5.0, 5.0), (24.0, 4.0), (49.0, 3.0), (80.0, 2.0), COMPLEXDUM]
        triMat = reshape([complex(CKC) :: COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM &
                                        , COMPLEXDUM,  COMPLEXDUM,  (2.0, 2.0), (3.0,  3.0),   (2.0, 2.0),   COMPLEXDUM,   COMPLEXDUM &
                                        , COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, (2.0,  2.0),   (5.0, 5.0),   COMPLEXDUM,   COMPLEXDUM &
                                        , COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM,   (3.0, 3.0),   COMPLEXDUM,   COMPLEXDUM &
                                        , COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM &
                                        , COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM &
                                        ], shape = [6, 7], order = [2, 1])
        solMat = [complex(CKC) :: COMPLEXDUM, COMPLEXDUM, (5.0, 5.0), (4.0, 4.0), (3.0, 3.0), (2.0, 2.0), COMPLEXDUM]
        genMat = genRef
        call disp%skip()
        call disp%show("triMat")
        call disp%show( triMat , format = cform )
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("ndim = 4; roffA = 1; coffA = 1; incB = 1;")
                        ndim = 4; roffA = 1; coffA = 1; incB = 1;
        call disp%show("call setMatMulTri(triMat, upperUnit, transUnit, genMat(3:6), ndim, roffA, coffA, incB) ! blas trsv contiguous interface.")
                        call setMatMulTri(triMat, upperUnit, transUnit, genMat(3:6), ndim, roffA, coffA, incB) ! blas trsv contiguous interface.
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("genMat - solMat")
        call disp%show( genMat - solMat , format = cform )
        call disp%skip()
        call disp%show("genMat = genRef")
                        genMat = genRef
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("call setMatMulTri(triMat(2:5, 2:5), upperUnit, transUnit, genMat(3:6)) ! blas trsv simplified interface.")
                        call setMatMulTri(triMat(2:5, 2:5), upperUnit, transUnit, genMat(3:6)) ! blas trsv simplified interface.
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("genMat - solMat")
        call disp%show( genMat - solMat , format = cform )
        call disp%skip()
        call disp%show("call setMatMulTri(triMat, upperUnit, transHerm, genMat(3:6), ndim, roffA, coffA, incB) ! blas trmv contiguous interface.")
                        call setMatMulTri(triMat, upperUnit, transHerm, genMat(3:6), ndim, roffA, coffA, incB) ! blas trmv contiguous interface.
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("genMat - genRef")
        call disp%show( genMat - genRef , format = cform )
        call disp%skip()
        call disp%show("genMat = solMat")
                        genMat = solMat
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("call setMatMulTri(triMat(2:5, 2:5), upperUnit, transHerm, genMat(3:6)) ! blas trmv simplified interface.")
                        call setMatMulTri(triMat(2:5, 2:5), upperUnit, transHerm, genMat(3:6)) ! blas trmv simplified interface.
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("genMat - genRef")
        call disp%show( genMat - genRef , format = cform )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! BLAS 3 TRSV / TRMV: triangular matrix-vector multiplication:    real.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        real(RKC), parameter :: DUM = huge(DUM)
        real(RKC), allocatable :: triMat(:,:), genMat(:), genRef(:), solMat(:)

        triMat = reshape([real(RKC):: DUM, DUM, DUM,  DUM,  DUM,  DUM &
                                    , DUM, DUM, DUM,  DUM,  DUM,  DUM &
                                    , DUM, DUM, DUM,  DUM,  DUM,  DUM &
                                    , DUM, 1.0, DUM,  DUM,  DUM,  DUM &
                                    , DUM, 2.0, 3.0,  DUM,  DUM,  DUM &
                                    , DUM, 3.0, 4.0,  3.0,  DUM,  DUM &
                                    , DUM, DUM, DUM,  DUM,  DUM,  DUM &
                                    , DUM, DUM, DUM,  DUM,  DUM,  DUM &
                                    , DUM, DUM, DUM,  DUM,  DUM,  DUM &
                                    , DUM, DUM, DUM,  DUM,  DUM,  DUM &
                                    ], shape = [10, 6], order = [2, 1])
        genRef = [real(RKC):: DUM, DUM, 1.0, DUM, DUM, 3.0, DUM, DUM, 11.0, DUM, DUM, 24.0, DUM]
        solMat = [real(RKC):: DUM, DUM, 1.0, DUM, DUM, 2.0, DUM, DUM,  3.0, DUM, DUM,  4.0, DUM]
        genMat = genRef
        call disp%skip()
        call disp%show("triMat")
        call disp%show( triMat , format = gform )
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("ndim = 4; roffA = 2; coffA = 1; incB = 3;")
                        ndim = 4; roffA = 2; coffA = 1; incB = 3;
        call disp%show("call setMatMulTri(triMat, lowerUnit, inversion, genMat(3:), ndim, roffA, coffA, incB) ! blas trsv contiguous interface.")
                        call setMatMulTri(triMat, lowerUnit, inversion, genMat(3:), ndim, roffA, coffA, incB) ! blas trsv contiguous interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - solMat")
        call disp%show( genMat - solMat , format = gform )
        call disp%skip()
        call disp%show("genMat = genRef")
                        genMat = genRef
        call disp%show("call setMatMulTri(triMat(3:6, 2:5), lowerUnit, inversion, genMat(3:size(genMat)-1:3)) ! blas trsv simplified interface.")
                        call setMatMulTri(triMat(3:6, 2:5), lowerUnit, inversion, genMat(3:size(genMat)-1:3)) ! blas trsv simplified interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - solMat")
        call disp%show( genMat - solMat , format = gform )
        call disp%skip()
        call disp%show("call setMatMulTri(triMat, lowerUnit, nothing, genMat(3:), ndim, roffA, coffA, incB) ! blas trmv contiguous interface.")
                        call setMatMulTri(triMat, lowerUnit, nothing, genMat(3:), ndim, roffA, coffA, incB) ! blas trmv contiguous interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - genRef")
        call disp%show( genMat - genRef , format = gform )
        call disp%skip()
        call disp%show("genMat = solMat")
                        genMat = solMat
        call disp%show("call setMatMulTri(triMat(3:6, 2:5), lowerUnit, nothing, genMat(3:size(genMat)-1:3)) ! blas trmv simplified interface.")
                        call setMatMulTri(triMat(3:6, 2:5), lowerUnit, nothing, genMat(3:size(genMat)-1:3)) ! blas trmv simplified interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - genRef")
        call disp%show( genMat - genRef , format = gform )
        call disp%skip()

        triMat = reshape([real(RKC):: 1.0,  2.0,  3.0,  2.0,  DUM &
                                    , DUM,  2.0,  2.0,  5.0,  DUM &
                                    , DUM,  DUM,  3.0,  3.0,  DUM &
                                    , DUM,  DUM,  DUM,  1.0,  DUM &
                                    , DUM,  DUM,  DUM,  DUM,  DUM &
                                    , DUM,  DUM,  DUM,  DUM,  DUM &
                                    , DUM,  DUM,  DUM,  DUM,  DUM &
                                    ], shape = [7, 5], order = [2, 1])
        genRef = [real(RKC) :: 5.0, 18.0, 32.0, 41.0]
        solMat = [real(RKC) :: 5.0, 4.0, 3.0, 2.0]
        genMat = genRef
        call disp%skip()
        call disp%show("triMat")
        call disp%show( triMat , format = gform )
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("ndim = 4; roffA = 0; coffA = 0; incB = 1;")
                        ndim = 4; roffA = 0; coffA = 0; incB = 1;
        call disp%show("call setMatMulTri(triMat, upperDiag, transOrth, genMat, ndim, roffA, coffA, incB) ! blas trsm contiguous interface.")
                        call setMatMulTri(triMat, upperDiag, transOrth, genMat, ndim, roffA, coffA, incB) ! blas trsm contiguous interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - solMat")
        call disp%show( genMat - solMat , format = gform )
        call disp%skip()
        call disp%show("genMat = genRef")
                        genMat = genRef
        call disp%show("call setMatMulTri(triMat(1:4, 1:4), upperDiag, transOrth, genMat) ! blas trsv simplified interface.")
                        call setMatMulTri(triMat(1:4, 1:4), upperDiag, transOrth, genMat) ! blas trsv simplified interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - solMat")
        call disp%show( genMat - solMat , format = gform )
        call disp%skip()
        call disp%show("call setMatMulTri(triMat, upperDiag, transSymm, genMat, ndim, roffA, coffA, incB) ! blas trmv contiguous interface.")
                        call setMatMulTri(triMat, upperDiag, transSymm, genMat, ndim, roffA, coffA, incB) ! blas trmv contiguous interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - genRef")
        call disp%show( genMat - genRef , format = gform )
        call disp%skip()
        call disp%show("genMat = solMat")
                        genMat = solMat
        call disp%show("call setMatMulTri(triMat(1:4, 1:4), upperDiag, transSymm, genMat) ! blas trmv simplified interface.")
                        call setMatMulTri(triMat(1:4, 1:4), upperDiag, transSymm, genMat) ! blas trmv simplified interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - genRef")
        call disp%show( genMat - genRef , format = gform )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! BLAS 3 TRSM / TRMM: triangular-general / general-triangular matrix multiplication: complex.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        complex(CKC) :: alpha, beta
        complex(CKC), parameter :: COMPLEXDUM = cmplx(huge(0._CKC), huge(0._CKC), CKC)
        complex(CKC), allocatable, dimension(:,:) :: triMat, genMat, genRef, solMat

        genRef = reshape([complex(CKC) :: COMPLEXDUM,    COMPLEXDUM,    COMPLEXDUM,  COMPLEXDUM,     COMPLEXDUM,    COMPLEXDUM &
                                        , COMPLEXDUM, (22.0, -41.0),  (7.0, -26.0),  (9.0, 0.0),  (-15.0, -3.0),  (-15.0, 8.0) &
                                        , COMPLEXDUM, (29.0, -18.0), (24.0, -10.0),  (9.0, 6.0), (-12.0, -24.0), (-19.0, -8.0) &
                                        , COMPLEXDUM,  (-15.0, 2.0), (-3.0, -21.0), (-2.0, 4.0),  (-4.0, -12.0), (-10.0, -6.0) &
                                        , COMPLEXDUM,    COMPLEXDUM,    COMPLEXDUM,  COMPLEXDUM,     COMPLEXDUM,    COMPLEXDUM &
                                        , COMPLEXDUM,    COMPLEXDUM,    COMPLEXDUM,  COMPLEXDUM,     COMPLEXDUM,    COMPLEXDUM &
                                        ], shape = [6, 6], order = [2, 1])
        triMat = reshape([complex(CKC) :: COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM &
                                        , COMPLEXDUM,  COMPLEXDUM, (2.0, -3.0),  COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM &
                                        , COMPLEXDUM,  COMPLEXDUM, (2.0, -4.0), (3.0, -1.0),   COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM &
                                        , COMPLEXDUM,  COMPLEXDUM, (2.0,  2.0), (1.0,  2.0),  (1.0,  1.0),   COMPLEXDUM,   COMPLEXDUM &
                                        , COMPLEXDUM,  COMPLEXDUM, (0.0,  0.0), (3.0, -1.0),  (0.0, -1.0), (-2.0,  1.0),   COMPLEXDUM &
                                        , COMPLEXDUM,  COMPLEXDUM, (2.0,  2.0), (4.0,  0.0), (-1.0,  2.0),  (2.0, -4.0), (-1.0, -4.0) &
                                        , COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM &
                                        , COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM &
                                        ], shape = [8, 7], order = [2, 1])
        solMat = reshape([complex(CKC) :: COMPLEXDUM,    COMPLEXDUM,    COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM &
                                        , COMPLEXDUM,    (3.0, 0.0),    (4.0, 0.0), (-1.0, -2.0), (-1.0, -1.0), (-1.0, -4.0) &
                                        , COMPLEXDUM,   (2.0, -1.0),    (1.0, 2.0), (-1.0, -3.0),   (0.0, 2.0),  (3.0, -4.0) &
                                        , COMPLEXDUM,   (-2.0, 1.0),  (-1.0, -3.0),  (-3.0, 1.0),   (0.0, 0.0),  (2.0, -2.0) &
                                        , COMPLEXDUM,    COMPLEXDUM,    COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM &
                                        , COMPLEXDUM,    COMPLEXDUM,    COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM,   COMPLEXDUM &
                                        ], shape = [6, 6], order = [2, 1])
        genMat = genRef
        call disp%skip()
        call disp%show("triMat")
        call disp%show( triMat , format = cform )
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("alpha = (1._CKC, 0._CKC); nrow = 3; ncol = 5; roffA = 1; coffA = 1; roffB = 1; coffB = 2;")
                        alpha = (1._CKC, 0._CKC); nrow = 3; ncol = 5; roffA = 1; coffA = 1; roffB = 1; coffB = 2;
        call disp%show("call setMatMulTri(genMat, triMat, lowerDiag, inversion, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trsm contiguous interface.")
                        call setMatMulTri(genMat, triMat, lowerDiag, inversion, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trsm contiguous interface.
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("genMat - solMat")
        call disp%show( genMat - solMat , format = cform )
        call disp%skip()
        call disp%show("genMat = genRef")
                        genMat = genRef
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("call setMatMulTri(genMat(2:4, 2:), triMat(2:6, 3:), lowerDiag, inversion) ! blas trsm simplified interface.")
                        call setMatMulTri(genMat(2:4, 2:), triMat(2:6, 3:), lowerDiag, inversion) ! blas trsm simplified interface.
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("genMat - solMat")
        call disp%show( genMat - solMat , format = cform )
        call disp%skip()
        call disp%show("call setMatMulTri(genMat, triMat, lowerDiag, nothing, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trmm contiguous interface.")
                        call setMatMulTri(genMat, triMat, lowerDiag, nothing, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trmm contiguous interface.
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("genMat - genRef")
        call disp%show( genMat - genRef , format = cform )
        call disp%skip()
        call disp%show("genMat = solMat")
                        genMat = solMat
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("call setMatMulTri(genMat(2:4, 2:), triMat(2:6, 3:), lowerDiag, nothing) ! blas trmm simplified interface.")
                        call setMatMulTri(genMat(2:4, 2:), triMat(2:6, 3:), lowerDiag, nothing) ! blas trmm simplified interface.
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("genMat - genRef")
        call disp%show( genMat - genRef , format = cform )
        call disp%skip()


        genMat = genRef
        call disp%skip()
        call disp%show("triMat")
        call disp%show( triMat , format = cform )
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("alpha = (1._CKC, 0._CKC); nrow = 3; ncol = 5; roffA = 1; coffA = 1; roffB = 2; coffB = 1;")
                        alpha = (1._CKC, 0._CKC); nrow = 3; ncol = 5; roffA = 1; coffA = 1; roffB = 2; coffB = 1;
        call disp%show("call setMatMulTri(genMat, transpose(triMat), upperDiag, inversion, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trsm contiguous interface.")
                        call setMatMulTri(genMat, transpose(triMat), upperDiag, inversion, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trsm contiguous interface.
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("genMat - solMat")
        call disp%show( genMat - solMat , format = cform )
        call disp%skip()
        call disp%show("genMat = genRef")
                        genMat = genRef
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("call setMatMulTri(genMat(2:4, 2:), transpose(triMat(2:6, 3:)), upperDiag, inversion) ! blas trsm simplified interface.")
                        call setMatMulTri(genMat(2:4, 2:), transpose(triMat(2:6, 3:)), upperDiag, inversion) ! blas trsm simplified interface.
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("genMat - solMat")
        call disp%show( genMat - solMat , format = cform )
        call disp%skip()
        call disp%show("call setMatMulTri(genMat, transpose(triMat), upperDiag, nothing, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trmm contiguous interface.")
                        call setMatMulTri(genMat, transpose(triMat), upperDiag, nothing, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trmm contiguous interface.
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("genMat - genRef")
        call disp%show( genMat - genRef , format = cform )
        call disp%skip()
        call disp%show("genMat = solMat")
                        genMat = solMat
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("call setMatMulTri(genMat(2:4, 2:), transpose(triMat(2:6, 3:)), upperDiag, nothing) ! blas trmm simplified interface.")
                        call setMatMulTri(genMat(2:4, 2:), transpose(triMat(2:6, 3:)), upperDiag, nothing) ! blas trmm simplified interface.
        call disp%show("genMat")
        call disp%show( genMat , format = cform )
        call disp%show("genMat - genRef")
        call disp%show( genMat - genRef , format = cform )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! BLAS 3 TRSM / TRMM: triangular-general / general-triangular matrix multiplication:    real.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        real(RKC):: alpha
        real(RKC), parameter :: DUM = huge(DUM), ONE_THIRD = 1._RKC / 3._RKC
        real(RKC), allocatable, dimension(:,:) :: triMat, genMat, genRef, solMat

        triMat = reshape([real(RKC):: 3.0, -1.0,  2.0,  2.0,  1.0 &
                                    , DUM, -2.0,  4.0, -1.0,  3.0 &
                                    , DUM,  DUM, -3.0,  0.0,  2.0 &
                                    , DUM,  DUM,  DUM,  4.0, -2.0 &
                                    , DUM,  DUM,  DUM,  DUM,  1.0 &
                                    , DUM,  DUM,  DUM,  DUM,  DUM &
                                    , DUM,  DUM,  DUM,  DUM,  DUM &
                                    ], shape = [7, 5], order = [2, 1])
        genRef = reshape([real(RKC)::   6.0, 10.0,  -2.0 &
                                    , -16.0, -1.0,   6.0 &
                                    ,  -2.0,  1.0,  -4.0 &
                                    ,  14.0,  0.0, -14.0 &
                                    ,  -1.0,  2.0,   1.0 &
                                    ,   DUM,  DUM,   DUM &
                                    ], shape = [6, 3], order = [2, 1])
        solMat = reshape([real(RKC)::  2.0, 3.0,  1.0 &
                                    ,  5.0, 5.0,  4.0 &
                                    ,  0.0, 1.0,  2.0 &
                                    ,  3.0, 1.0, -3.0 &
                                    , -1.0, 2.0,  1.0 &
                                    ,  DUM, DUM,  DUM &
                                    ], shape = [6, 3], order = [2, 1])
        genMat = genRef
        call disp%skip()
        call disp%show("triMat")
        call disp%show( triMat , format = gform )
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("alpha = 1._RKC; nrow = 5; ncol = 3; roffA = 0; coffA = 0; roffB = 0; coffB = 0;")
                        alpha = 1._RKC; nrow = 5; ncol = 3; roffA = 0; coffA = 0; roffB = 0; coffB = 0;
        call disp%show("call setMatMulTri(triMat, upperDiag, inversion, genMat, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trsm contiguous interface.")
                        call setMatMulTri(triMat, upperDiag, inversion, genMat, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trsm contiguous interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - solMat")
        call disp%show( genMat - solMat , format = gform )
        call disp%skip()
        call disp%show("genMat = genRef")
                        genMat = genRef
        call disp%show("call setMatMulTri(triMat(1:5, 1:5), upperDiag, inversion, genMat(1:5, 1:3), alpha) ! blas trsm simplified interface.")
                        call setMatMulTri(triMat(1:5, 1:5), upperDiag, inversion, genMat(1:5, 1:3), alpha) ! blas trsm simplified interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - solMat")
        call disp%show( genMat - solMat , format = gform )
        call disp%skip()
        call disp%show("call setMatMulTri(triMat, upperDiag, nothing, genMat, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trmm contiguous interface.")
                        call setMatMulTri(triMat, upperDiag, nothing, genMat, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trmm contiguous interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - genRef")
        call disp%show( genMat - genRef , format = gform )
        call disp%skip()
        call disp%show("genMat = solMat")
                        genMat = solMat
        call disp%show("call setMatMulTri(triMat(1:5, 1:5), upperDiag, nothing, genMat(1:5, 1:3), alpha) ! blas trmm simplified interface.")
                        call setMatMulTri(triMat(1:5, 1:5), upperDiag, nothing, genMat(1:5, 1:3), alpha) ! blas trmm simplified interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - genRef")
        call disp%show( genMat - genRef , format = gform )
        call disp%skip()

        triMat = reshape([real(RKC):: -1.0, -4.0, -2.0,  2.0,  3.0 &
                                    ,  DUM, -2.0,  2.0,  2.0,  2.0 &
                                    ,  DUM,  DUM, -3.0, -1.0,  4.0 &
                                    ,  DUM,  DUM,  DUM,  1.0,  0.0 &
                                    ,  DUM,  DUM,  DUM,  DUM, -2.0 &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ], shape = [7, 5], order = [2, 1])
        genRef = reshape([real(RKC):: -1.0, -2.0,  -3.0,  -4.0 &
                                    ,  2.0, -2.0, -14.0, -12.0 &
                                    , 10.0,  5.0,  -8.0,  -7.0 &
                                    , 14.0, 15.0,   1.0,   8.0 &
                                    , -3.0,  4.0,   3.0,  16.0 &
                                    ,  DUM,  DUM,   DUM,   DUM &
                                    ], shape = [6, 4], order = [2, 1])
        solMat = reshape([real(RKC)::  1._RKC,   2._RKC,              3._RKC,  4._RKC             &
                                    , -3._RKC,  -3._RKC,              1._RKC, -2._RKC             &
                                    , -6._RKC,  -5._RKC,  1._RKC + ONE_THIRD, -2._RKC + ONE_THIRD &
                                    , 12._RKC,  12._RKC, -6._RKC + ONE_THIRD,  2._RKC + ONE_THIRD &
                                    ,-12._RKC, -12._RKC, +7._RKC - ONE_THIRD, -7._RKC - ONE_THIRD &
                                    ,     DUM,      DUM,                 DUM,                 DUM &
                                    ], shape = [6, 4], order = [2, 1])
        genMat = genRef
        call disp%skip()
        call disp%show("triMat")
        call disp%show( triMat , format = gform )
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("alpha = 1._RKC; nrow = 5; ncol = 4; roffA = 0; coffA = 0; roffB = 0; coffB = 0;")
                        alpha = 1._RKC; nrow = 5; ncol = 4; roffA = 0; coffA = 0; roffB = 0; coffB = 0;
        call disp%show("call setMatMulTri(triMat, upperDiag, transOrth, genMat, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trsm contiguous interface.")
                        call setMatMulTri(triMat, upperDiag, transOrth, genMat, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trsm contiguous interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - solMat")
        call disp%show( genMat - solMat , format = gform )
        call disp%skip()
        call disp%show("genMat = genRef")
                        genMat = genRef
        call disp%show("call setMatMulTri(triMat(1:5, 1:5), upperDiag, transOrth, genMat(1:5, 1:4)) ! trsm simplified interface.")
                        call setMatMulTri(triMat(1:5, 1:5), upperDiag, transOrth, genMat(1:5, 1:4)) ! trsm simplified interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - solMat")
        call disp%show( genMat - solMat , format = gform )
        call disp%skip()
        call disp%show("call setMatMulTri(triMat, upperDiag, transSymm, genMat, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trmm contiguous interface.")
                        call setMatMulTri(triMat, upperDiag, transSymm, genMat, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trmm contiguous interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - genRef")
        call disp%show( genMat - genRef , format = gform )
        call disp%skip()
        call disp%show("genMat = solMat")
                        genMat = solMat
        call disp%show("call setMatMulTri(triMat(1:5, 1:5), upperDiag, transSymm, genMat(1:5, 1:4)) ! trmm simplified interface.")
                        call setMatMulTri(triMat(1:5, 1:5), upperDiag, transSymm, genMat(1:5, 1:4)) ! trmm simplified interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - genRef")
        call disp%show( genMat - genRef , format = gform )
        call disp%skip()

        genRef = reshape([real(RKC):: DUM,  DUM,  DUM,  DUM, DUM,  DUM &
                                    , DUM, 10.0,  4.0,  0.0, 0.0,  1.0 &
                                    , DUM, 10.0, 14.0, -4.0, 6.0, -3.0 &
                                    , DUM, -8.0,  2.0, -5.0, 4.0, -2.0 &
                                    , DUM,  DUM,  DUM,  DUM, DUM,  DUM &
                                    ], shape = [5, 6], order = [2, 1])
        triMat = reshape([real(RKC):: DUM, DUM, DUM,  DUM,  DUM,  DUM &
                                    , DUM, DUM, DUM,  DUM,  DUM,  DUM &
                                    , DUM, DUM, DUM,  DUM,  DUM,  DUM &
                                    , DUM, 2.0, DUM,  DUM,  DUM,  DUM &
                                    , DUM, 2.0, 3.0,  DUM,  DUM,  DUM &
                                    , DUM, 2.0, 1.0,  1.0,  DUM,  DUM &
                                    , DUM, 0.0, 3.0,  0.0, -2.0,  DUM &
                                    , DUM, 2.0, 4.0, -1.0,  2.0, -1.0 &
                                    , DUM, DUM, DUM,  DUM,  DUM,  DUM &
                                    , DUM, DUM, DUM,  DUM,  DUM,  DUM &
                                    ], shape = [10, 6], order = [2, 1])
        solMat = reshape([real(RKC):: DUM,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    , DUM,  3.0,  4.0, -1.0, -1.0, -1.0 &
                                    , DUM,  2.0,  1.0, -1.0,  0.0,  3.0 &
                                    , DUM, -2.0, -1.0, -3.0,  0.0,  2.0 &
                                    , DUM,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ], shape = [5, 6], order = [2, 1])
        genMat = genRef
        call disp%skip()
        call disp%show("triMat")
        call disp%show( triMat , format = gform )
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("alpha = 1._RKC; nrow = 3; ncol = 5; roffA = 1; coffA = 1; roffB = 3; coffB = 1;")
                        alpha = 1._RKC; nrow = 3; ncol = 5; roffA = 1; coffA = 1; roffB = 3; coffB = 1;
        call disp%show("call setMatMulTri(genMat, triMat, lowerDiag, inversion, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trsm contiguous interface.")
                        call setMatMulTri(genMat, triMat, lowerDiag, inversion, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trsm contiguous interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - solMat")
        call disp%show( genMat - solMat , format = gform )
        call disp%skip()
        call disp%show("genMat = genRef")
                        genMat = genRef
        call disp%show("call setMatMulTri(genMat(2:4, 2:6), triMat(4:8, 2:6), lowerDiag, inversion) ! trsm simplified interface.")
                        call setMatMulTri(genMat(2:4, 2:6), triMat(4:8, 2:6), lowerDiag, inversion) ! trsm simplified interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - solMat")
        call disp%show( genMat - solMat , format = gform )
        call disp%skip()
        call disp%show("call setMatMulTri(genMat, triMat, lowerDiag, nothing, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trmm contiguous interface.")
                        call setMatMulTri(genMat, triMat, lowerDiag, nothing, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trmm contiguous interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - genRef")
        call disp%show( genMat - genRef , format = gform )
        call disp%skip()
        call disp%show("genMat = solMat")
                        genMat = solMat
        call disp%show("call setMatMulTri(genMat(2:4, 2:6), triMat(4:8, 2:6), lowerDiag, nothing) ! trmm simplified interface.")
                        call setMatMulTri(genMat(2:4, 2:6), triMat(4:8, 2:6), lowerDiag, nothing) ! trmm simplified interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - genRef")
        call disp%show( genMat - genRef , format = gform )
        call disp%skip()

        genRef = reshape([real(RKC):: DUM,  DUM,  DUM,  DUM, DUM,  DUM &
                                    , DUM,  DUM,  DUM,  DUM, DUM,  DUM &
                                    , 1.0,  4.0, -2.0, 10.0, 2.0, -6.0 &
                                    , DUM,  DUM,  DUM,  DUM, DUM,  DUM &
                                    ], shape = [4, 6], order = [2, 1])
        triMat = reshape([real(RKC):: DUM, DUM, DUM,  DUM,  DUM,  DUM &
                                    , DUM, 2.0, -3.0, 1.0,  2.0,  4.0 &
                                    , DUM, DUM,  0.0, 1.0,  1.0, -2.0 &
                                    , DUM, DUM,  DUM, 4.0, -1.0,  1.0 &
                                    , DUM, DUM,  DUM, DUM,  0.0, -1.0 &
                                    , DUM, DUM,  DUM, DUM,  DUM,  2.0 &
                                    , DUM, DUM,  DUM, DUM,  DUM,  DUM &
                                    , DUM, DUM,  DUM, DUM,  DUM,  DUM &
                                    ], shape = [8, 6], order = [2, 1])
        solMat = reshape([real(RKC):: DUM,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    , DUM,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    , 1.0,  2.0,  1.0,  3.0, -1.0, -2.0 &
                                    , DUM,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ], shape = [4, 6], order = [2, 1])
        genMat = genRef
        call disp%skip()
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("triMat")
        call disp%show( triMat , format = gform )
        call disp%show("alpha = 1._RKC; nrow = 1; ncol = 6; roffA = 2; coffA = 0; roffB = 1; coffB = 0;")
                        alpha = 1._RKC; nrow = 1; ncol = 6; roffA = 2; coffA = 0; roffB = 1; coffB = 0;
        call disp%show("call setMatMulTri(genMat, triMat, upperUnit, inversion, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trsm contiguous interface.")
                        call setMatMulTri(genMat, triMat, upperUnit, inversion, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trsm contiguous interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - solMat")
        call disp%show( genMat - solMat , format = gform )
        call disp%skip()
        call disp%show("genMat = genRef")
                        genMat = genRef
        call disp%show("call setMatMulTri(genMat(3:3, 1:), triMat(2:7, 1:), upperUnit, inversion) ! trsm simplified interface.")
                        call setMatMulTri(genMat(3:3, 1:), triMat(2:7, 1:), upperUnit, inversion) ! trsm simplified interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - solMat")
        call disp%show( genMat - solMat , format = gform )
        call disp%skip()
        call disp%show("call setMatMulTri(genMat, triMat, upperUnit, nothing, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trmm contiguous interface.")
                        call setMatMulTri(genMat, triMat, upperUnit, nothing, alpha, nrow, ncol, roffA, coffA, roffB, coffB) ! blas trmm contiguous interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - genRef")
        call disp%show( genMat - genRef , format = gform )
        call disp%skip()
        call disp%show("genMat = solMat")
                        genMat = solMat
        call disp%show("call setMatMulTri(genMat(3:3, 1:), triMat(2:7, 1:), upperUnit, nothing) ! trmm simplified interface.")
                        call setMatMulTri(genMat(3:3, 1:), triMat(2:7, 1:), upperUnit, nothing) ! trmm simplified interface.
        call disp%show("genMat")
        call disp%show( genMat , format = gform )
        call disp%show("genMat - genRef")
        call disp%show( genMat - genRef , format = gform )
        call disp%skip()

    end block

end program example