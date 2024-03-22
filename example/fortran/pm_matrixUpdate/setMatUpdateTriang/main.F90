program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKC => RKS ! all processor type kinds are supported.
    use pm_kind, only: CKC => CKS ! all processor type kinds are supported.
    use pm_matrixUpdate, only: lowDia, uppDia
    use pm_matrixUpdate, only: transSymm, transHerm
    use pm_matrixUpdate, only: setMatUpdateTriang
    use pm_io, only: getFormat
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp

    real(RKC)       , parameter     :: DUM = -huge(0._RKC)
    complex(CKC)    , parameter     :: CMPLX_DUMM = cmplx(-huge(0._CKC), -huge(0._CKC), CKC)
    integer(IK)                     :: ndim, ndum, roff, coff, roffA, coffA
    character(:, SK), allocatable   :: cform
    cform = getFormat([cmplx(0., 0., CKC)], ed = SK_'f', signed = .true.)

    disp = display_type(file = "main.out.F90")

    block

        real(RKC) :: alpha, beta
        real(RKC), allocatable :: ref(:,:), mat(:,:), matA(:,:), VecX(:), VecY(:)

        matA = reshape( [ real(RKC) :: 0.0,  8.0 &
                                    ,  1.0,  9.0 &
                                    ,  2.0, 10.0 &
                                    ,  3.0, 11.0 &
                                    ,  4.0, 12.0 &
                                    ,  5.0, 13.0 &
                                    ,  6.0, 14.0 &
                                    ,  7.0, 15.0 &
                                    ,  DUM,  DUM &
                                    ], shape = [9, 2], order = [2, 1])
        mat = reshape( [ real(RKC)  :: 0.0,  1.0,  3.0,  6.0,  10.0,  15.0,  21.0,  28.0 &
                                    ,  DUM,  2.0,  4.0,  7.0,  11.0,  16.0,  22.0,  29.0 &
                                    ,  DUM,  DUM,  5.0,  8.0,  12.0,  17.0,  23.0,  30.0 &
                                    ,  DUM,  DUM,  DUM,  9.0,  13.0,  18.0,  24.0,  31.0 &
                                    ,  DUM,  DUM,  DUM,  DUM,  14.0,  19.0,  25.0,  32.0 &
                                    ,  DUM,  DUM,  DUM,  DUM,   DUM,  20.0,  26.0,  33.0 &
                                    ,  DUM,  DUM,  DUM,  DUM,   DUM,   DUM,  27.0,  34.0 &
                                    ,  DUM,  DUM,  DUM,  DUM,   DUM,   DUM,   DUM,  35.0 &
                                    ,  DUM,  DUM,  DUM,  DUM,   DUM,   DUM,   DUM,   DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,   DUM,   DUM,   DUM,   DUM &
                                    ], shape = [10, 8], order = [2, 1])
        ref = reshape( [ real(RKC) :: 64.0, 73.0,  83.0,  94.0, 106.0, 119.0, 133.0, 148.0 &
                                    ,  DUM, 84.0,  96.0, 109.0, 123.0, 138.0, 154.0, 171.0 &
                                    ,  DUM,  DUM, 109.0, 124.0, 140.0, 157.0, 175.0, 194.0 &
                                    ,  DUM,  DUM,   DUM, 139.0, 157.0, 176.0, 196.0, 217.0 &
                                    ,  DUM,  DUM,   DUM,   DUM, 174.0, 195.0, 217.0, 240.0 &
                                    ,  DUM,  DUM,   DUM,   DUM,   DUM, 214.0, 238.0, 263.0 &
                                    ,  DUM,  DUM,   DUM,   DUM,   DUM,   DUM, 259.0, 286.0 &
                                    ,  DUM,  DUM,   DUM,   DUM,   DUM,   DUM,   DUM, 309.0 &
                                    ,  DUM,  DUM,   DUM,   DUM,   DUM,   DUM,   DUM,   DUM &
                                    ,  DUM,  DUM,   DUM,   DUM,   DUM,   DUM,   DUM,   DUM &
                                    ], shape = [10, 8], order = [2, 1])
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA )
        call disp%show("mat")
        call disp%show( mat )
        call disp%show("alpha = 1._RKC; beta = 1._RKC; ndim = 8; ndum = 2; roff = 0; coff = 0; roffA = 0; coffA = 0;")
                        alpha = 1._RKC; beta = 1._RKC; ndim = 8; ndum = 2; roff = 0; coff = 0; roffA = 0; coffA = 0;
        call disp%show("call setMatUpdateTriang(mat, uppDia, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)")
                        call setMatUpdateTriang(mat, uppDia, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
        call disp%show("mat")
        call disp%show( mat )
        call disp%show("mat - ref")
        call disp%show( mat - ref )
        call disp%skip()

        matA = reshape( [ real(RKC) :: 0.0,  3.0,  6.0,  9.0, 12.0, 15.0, 18.0, 21.0 &
                                    ,  1.0,  4.0,  7.0, 10.0, 13.0, 16.0, 19.0, 22.0 &
                                    ,  2.0,  5.0,  8.0, 11.0, 14.0, 17.0, 20.0, 23.0 &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ], shape = [4, 8], order = [2, 1])
        mat = reshape( [ real(RKC) :: 0.0,  DUM,  DUM,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    , 1.0,  8.0,  DUM,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    , 2.0,  9.0, 15.0,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    , 3.0, 10.0, 16.0, 21.0,  DUM,  DUM,  DUM,  DUM &
                                    , 4.0, 11.0, 17.0, 22.0, 26.0,  DUM,  DUM,  DUM &
                                    , 5.0, 12.0, 18.0, 23.0, 27.0, 30.0,  DUM,  DUM &
                                    , 6.0, 13.0, 19.0, 24.0, 28.0, 31.0, 33.0,  DUM &
                                    , 7.0, 14.0, 20.0, 25.0, 29.0, 32.0, 34.0, 35.0 &
                                    ], shape = [8, 8], order = [2, 1])
        ref = reshape( [ real(RKC) ::  5.0,   DUM,   DUM,   DUM,   DUM,    DUM,    DUM,    DUM &
                                    , 15.0,  58.0,   DUM,   DUM,   DUM,    DUM,    DUM,    DUM &
                                    , 25.0,  95.0, 164.0,   DUM,   DUM,    DUM,    DUM,    DUM &
                                    , 35.0, 132.0, 228.0, 323.0,   DUM,    DUM,    DUM,    DUM &
                                    , 45.0, 169.0, 292.0, 414.0, 535.0,    DUM,    DUM,    DUM &
                                    , 55.0, 206.0, 356.0, 505.0, 653.0,  800.0,    DUM,    DUM &
                                    , 65.0, 243.0, 420.0, 596.0, 771.0,  945.0, 1118.0,    DUM &
                                    , 75.0, 280.0, 484.0, 687.0, 889.0, 1090.0, 1290.0, 1489.0 &
                                    ], shape = [8, 8], order = [2, 1])
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA )
        call disp%show("mat")
        call disp%show( mat )
        call disp%show("alpha = 1._RKC; beta = 1._RKC; ndim = 8; ndum = 3; roff = 0; coff = 0; roffA = 0; coffA = 0;")
                        alpha = 1._RKC; beta = 1._RKC; ndim = 8; ndum = 3; roff = 0; coff = 0; roffA = 0; coffA = 0;
        call disp%show("call setMatUpdateTriang(mat, lowDia, matA, transSymm, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)")
                        call setMatUpdateTriang(mat, lowDia, matA, transSymm, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
        call disp%show("mat")
        call disp%show( mat )
        call disp%show("mat - ref")
        call disp%show( mat - ref )
        call disp%skip()

    end block

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Complex Symmetric matrix update.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    block

        complex(CKC) :: alpha, beta
        complex(CKC), allocatable :: ref(:,:), mat(:,:), matA(:,:), VecX(:), VecY(:)

        matA = reshape( [ complex(CKC) :: (2.0, 0.0), (3.0, 2.0), (4.0, 1.0), (1.0, 7.0), (0.0, 0.0) &
                                        , (3.0, 3.0), (8.0, 0.0), (2.0, 5.0), (2.0, 4.0), (1.0, 2.0) &
                                        , (1.0, 3.0), (2.0, 1.0), (6.0, 0.0), (3.0, 2.0), (2.0, 2.0) &
                                        ], shape = [3, 5], order = [2, 1])
        mat = reshape( [ complex(CKC) :: (2.0, 1.0), (1.0, 9.0), (4.0, 5.0) &
                                        , CMPLX_DUMM, (3.0, 1.0), (6.0, 7.0) &
                                        , CMPLX_DUMM, CMPLX_DUMM, (8.0, 1.0) &
                                        , CMPLX_DUMM, CMPLX_DUMM, CMPLX_DUMM &
                                        ], shape = [4, 3], order = [2, 1])
        ref = reshape( [ complex(CKC) :: (-57.0, 13.0), (-63.0, 79.0), (-24.0,  70.0) &
                                        ,    CMPLX_DUMM, (-28.0, 90.0), (-55.0, 103.0) &
                                        ,    CMPLX_DUMM,    CMPLX_DUMM,  (13.0,  75.0) &
                                        ,    CMPLX_DUMM,    CMPLX_DUMM,     CMPLX_DUMM &
                                        ], shape = [4, 3], order = [2, 1])
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = cform )
        call disp%show("mat")
        call disp%show( mat , format = cform )
        call disp%show("alpha = (1._CKC, 1._CKC); beta = (1._CKC, 1._CKC); ndim = 3; ndum = 5; roff = 0; coff = 0; roffA = 0; coffA = 0;")
                        alpha = (1._CKC, 1._CKC); beta = (1._CKC, 1._CKC); ndim = 3; ndum = 5; roff = 0; coff = 0; roffA = 0; coffA = 0;
        call disp%show("call setMatUpdateTriang(mat, uppDia, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)")
                        call setMatUpdateTriang(mat, uppDia, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
        call disp%show("mat")
        call disp%show( mat , format = cform )
        call disp%show("mat - ref")
        call disp%show( mat - ref , format = cform )
        call disp%skip()

    end block

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Complex hermitian matrix update.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    block

        real(CKC) :: alpha, beta
        complex(CKC), allocatable :: ref(:,:), mat(:,:), matA(:,:), VecX(:), VecY(:)

        matA = reshape( [ complex(CKC) :: (2.0, 0.0), (3.0, 2.0), (4.0, 1.0) &
                                        , (3.0, 3.0), (8.0, 0.0), (2.0, 5.0) &
                                        , (1.0, 3.0), (2.0, 1.0), (6.0, 0.0) &
                                        , (3.0, 3.0), (8.0, 0.0), (2.0, 5.0) &
                                        , (1.0, 9.0), (3.0, 0.0), (6.0, 7.0) &
                                        ], shape = [5, 3], order = [2, 1])
        mat = reshape( [ complex(CKC) :: (6.0, DUM),  CMPLX_DUMM, CMPLX_DUMM &
                                        , (3.0, 4.0), (10.0, DUM), CMPLX_DUMM &
                                        , (9.0, 1.0), (12.0, 2.0), (3.0, DUM) &
                                        , CMPLX_DUMM,  CMPLX_DUMM, CMPLX_DUMM &
                                        ], shape = [4, 3], order = [2, 1])
        ref = reshape( [ complex(CKC) :: (138.0,  0.0),     CMPLX_DUMM,   CMPLX_DUMM &
                                        ,  (65.0, 80.0), (165.0,   0.0),   CMPLX_DUMM &
                                        , (134.0, 46.0),  (88.0, -88.0), (199.0, 0.0) &
                                        ,    CMPLX_DUMM,     CMPLX_DUMM,   CMPLX_DUMM &
                                        ], shape = [4, 3], order = [2, 1])
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = cform )
        call disp%show("mat")
        call disp%show( mat , format = cform )
        call disp%show("alpha = (1._CKC, 1._CKC); beta = (1._CKC, 1._CKC); ndim = 3; ndum = 5; roff = 0; coff = 0; roffA = 0; coffA = 0;")
                        alpha = (1._CKC, 1._CKC); beta = (1._CKC, 1._CKC); ndim = 3; ndum = 5; roff = 0; coff = 0; roffA = 0; coffA = 0;
        call disp%show("call setMatUpdateTriang(mat, lowDia, matA, transHerm, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)")
                        call setMatUpdateTriang(mat, lowDia, matA, transHerm, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
        call disp%show("mat")
        call disp%show( mat , format = cform )
        call disp%show("mat - ref ! Note the imaginary parts of the diagonals are set to zero on output.")
        call disp%show( mat - ref , format = cform )
        call disp%skip()


        matA = reshape( [ complex(CKC) :: (2.0, 0.0), (3.0, 2.0), (4.0, 1.0) &
                                        , (3.0, 3.0), (8.0, 0.0), (2.0, 5.0) &
                                        , (1.0, 3.0), (2.0, 1.0), (6.0, 0.0) &
                                        , (3.0, 3.0), (8.0, 0.0), (2.0, 5.0) &
                                        , (1.0, 9.0), (3.0, 0.0), (6.0, 7.0) &
                                        ], shape = [3, 5], order = [1, 2])
        matA = conjg(matA)
        mat = reshape( [ complex(CKC)  :: (6.0, DUM),  CMPLX_DUMM, CMPLX_DUMM &
                                        , (3.0, 4.0), (10.0, DUM), CMPLX_DUMM &
                                        , (9.0, 1.0), (12.0, 2.0), (3.0, DUM) &
                                        , CMPLX_DUMM,  CMPLX_DUMM, CMPLX_DUMM &
                                        ], shape = [4, 3], order = [2, 1])
        ref = reshape( [ complex(CKC) :: (138.0,  0.0),     CMPLX_DUMM,   CMPLX_DUMM &
                                        ,  (65.0, 80.0), (165.0,   0.0),   CMPLX_DUMM &
                                        , (134.0, 46.0),  (88.0, -88.0), (199.0, 0.0) &
                                        ,    CMPLX_DUMM,     CMPLX_DUMM,   CMPLX_DUMM &
                                        ], shape = [4, 3], order = [2, 1])
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = cform )
        call disp%show("mat")
        call disp%show( mat , format = cform )
        call disp%show("alpha = 1._CKC; beta = 1._CKC; ndim = 3; ndum = 5; roff = 0; coff = 0; roffA = 0; coffA = 0;")
                        alpha = 1._CKC; beta = 1._CKC; ndim = 3; ndum = 5; roff = 0; coff = 0; roffA = 0; coffA = 0;
        call disp%show("call setMatUpdateTriang(mat, lowDia, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)")
                        call setMatUpdateTriang(mat, lowDia, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
        call disp%show("mat")
        call disp%show( mat , format = cform )
        call disp%show("mat - ref ! Note the imaginary parts of the diagonals are set to zero on output.")
        call disp%show( mat - ref , format = cform )
        call disp%skip()


        matA = reshape( [ complex(CKC) :: (2.0, 0.0), (3.0, 2.0), (4.0, 1.0) &
                                        , (3.0, 3.0), (8.0, 0.0), (2.0, 5.0) &
                                        , (1.0, 3.0), (2.0, 1.0), (6.0, 0.0) &
                                        , (3.0, 3.0), (8.0, 0.0), (2.0, 5.0) &
                                        , (1.0, 9.0), (3.0, 0.0), (6.0, 7.0) &
                                        ], shape = [3, 5])
        matA = conjg(matA)
        mat = reshape( [ complex(CKC)  :: (6.0, DUM),  CMPLX_DUMM, CMPLX_DUMM &
                                        , (3.0, 4.0), (10.0, DUM), CMPLX_DUMM &
                                        , (9.0, 1.0), (12.0, 2.0), (3.0, DUM) &
                                        , CMPLX_DUMM,  CMPLX_DUMM, CMPLX_DUMM &
                                        ], shape = [3, 4])
        ref = reshape( [ complex(CKC)  :: (138., .0), (65., -72.), (134., -44.), CMPLX_DUMM &
                                        , CMPLX_DUMM,   (165.,0.),   (88., 92.), CMPLX_DUMM &
                                        , CMPLX_DUMM,  CMPLX_DUMM,    (199.,0.), CMPLX_DUMM &
                                        ], shape = [3, 4], order = [2, 1])
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = cform )
        call disp%show("mat")
        call disp%show( mat , format = cform )
        call disp%show("alpha = 1._CKC; beta = 1._CKC; ndim = 3; ndum = 5; roff = 0; coff = 0; roffA = 0; coffA = 0;")
                        alpha = 1._CKC; beta = 1._CKC; ndim = 3; ndum = 5; roff = 0; coff = 0; roffA = 0; coffA = 0;
        call disp%show("call setMatUpdateTriang(mat, uppDia, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)")
                        call setMatUpdateTriang(mat, uppDia, matA, alpha, beta, ndim, ndum, roff, coff, roffA, coffA)
        call disp%show("mat")
        call disp%show( mat , format = cform )
        call disp%show("mat - ref ! Note the imaginary parts of the diagonals are set to zero on output.")
        call disp%show( mat - ref , format = cform )
        call disp%skip()

    end block

end program example