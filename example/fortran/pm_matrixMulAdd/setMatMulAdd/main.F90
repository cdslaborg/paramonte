program example

    use pm_kind, only: SK, IK, LK, TKG => RKS
    use pm_matrixMulAdd, only: symmetric, hermitian
    use pm_matrixMulAdd, only: transSymm, transHerm
    use pm_matrixMulAdd, only: uppDia, lowDia
    use pm_matrixMulAdd, only: setMatMulAdd
    use pm_arrayReverse, only: getReversed
    use pm_io, only: display_type
    use pm_io, only: getFormat

    implicit none

    type(display_type) :: disp
    character(:, SK), allocatable :: cform, rform, iform
    integer(IK) :: nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC, incB, incC
    cform = getFormat([cmplx(0., 0., TKG)], ed = SK_'f', signed = .true.)
    rform = getFormat([real(0., TKG)], ed = SK_'f', signed = .true.)
    iform = getFormat([0_IK], ed = SK_'i', signed = .true.)

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! BLAS 2 - GEMV: General matrix-vector multiplications - integer")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        use pm_kind, only: TKG => IKS
        integer(TKG) :: alpha, beta
        integer(TKG), parameter :: DUM = huge(DUM)
        integer(TKG), allocatable :: matA(:,:), matB(:), matC(:), refC(:)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape([ integer(TKG) ::  DUM,  DUM,  DUM,  DUM,  DUM &
                                       ,  DUM,  DUM,  1.0,  2.0,  3.0 &
                                       ,  DUM,  DUM,  2.0,  2.0,  4.0 &
                                       ,  DUM,  DUM,  3.0,  2.0,  2.0 &
                                       ,  DUM,  DUM,  4.0,  2.0,  1.0 &
                                       ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                       ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                       ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                       ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                       ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                       ,  DUM,  DUM,  DUM,  DUM,  DUM ], shape = [11, 5], order = [2, 1])
        matB = [ integer(TKG) ::  DUM,  DUM,  3.0,  2.0,  1.0 ]
        matC = [ integer(TKG) ::  4.0,  DUM,  5.0,  DUM,  2.0, DUM,  3.0 ]
        refC = [ integer(TKG) :: 14.0,  DUM, 19.0,  DUM, 17.0, DUM, 20.0 ]
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = iform )
        call disp%show("matB")
        call disp%show( matB , format = iform )
        call disp%show("matC")
        call disp%show( matC , format = iform )
        call disp%show("alpha = 1._TKG; beta = 1._TKG; nrow = 4; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 2;")
                        alpha = 1._TKG; beta = 1._TKG; nrow = 4; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 2;
        call disp%skip()
        call disp%show("call setMatMulAdd(matA, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC) ! full contiguous interface.")
                        call setMatMulAdd(matA, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
        call disp%show("matC")
        call disp%show( matC , format = iform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = iform )
        call disp%skip()
        call disp%show("matC = [ integer(TKG) ::  4.0,  DUM,  5.0,  DUM,  2.0, DUM,  3.0 ]")
                        matC = [ integer(TKG) ::  4.0,  DUM,  5.0,  DUM,  2.0, DUM,  3.0 ]
        call disp%show("call setMatMulAdd(matA(2:5, 3:5), matB(3:), matC(1::incC)) ! simplified interface.")
                        call setMatMulAdd(matA(2:5, 3:5), matB(3:), matC(1::incC))
        call disp%show("matC")
        call disp%show( matC , format = iform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = iform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape([ integer(TKG)  ::  DUM,  DUM,  DUM,  DUM,  DUM &
                                        ,  DUM,  DUM,  1.0,  2.0,  3.0 &
                                        ,  DUM,  DUM,  2.0,  2.0,  4.0 &
                                        ,  DUM,  DUM,  3.0,  2.0,  2.0 &
                                        ,  DUM,  DUM,  4.0,  2.0,  1.0 &
                                        ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                        ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                        ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                        ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                        ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                        ,  DUM,  DUM,  DUM,  DUM,  DUM ], shape = [11, 5], order = [2, 1])
        matB = [ integer(TKG) ::  DUM,  DUM,  3.0,  2.0,  1.0,  4.0 ]
        matC = [ integer(TKG) ::  1.0,  DUM,  2.0,  DUM,  3.0 ]
        refC = [ integer(TKG) :: 28.0,  DUM, 24.0,  DUM, 29.0 ]
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = iform )
        call disp%show("matB")
        call disp%show( matB , format = iform )
        call disp%show("matC")
        call disp%show( matC , format = iform )
        call disp%show("alpha = 1._TKG; beta = 2._TKG; nrow = 4; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 2;")
                        alpha = 1._TKG; beta = 2._TKG; nrow = 4; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 2;
        call disp%skip()
        call disp%show("call setMatMulAdd(matA, transSymm, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC) ! full contiguous interface.")
                        call setMatMulAdd(matA, transSymm, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
        call disp%show("matC")
        call disp%show( matC , format = iform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = iform )
        call disp%skip()
        call disp%show("matC = [ integer(TKG) ::  1.0,  DUM,  2.0,  DUM,  3.0 ]")
                        matC = [ integer(TKG) ::  1.0,  DUM,  2.0,  DUM,  3.0 ]
        call disp%show("call setMatMulAdd(matA(2:5, 3:5), transSymm, matB(3:), matC(1::incC), beta = beta) ! simplified interface.")
                        call setMatMulAdd(matA(2:5, 3:5), transSymm, matB(3:), matC(1::incC), beta = beta)
        call disp%show("matC")
        call disp%show( matC , format = iform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = iform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape([ integer(TKG)  ::  DUM,  DUM,  DUM,  DUM,  DUM &
                                        ,  DUM,  DUM,  1.0,  2.0,  3.0 &
                                        ,  DUM,  DUM,  2.0,  2.0,  4.0 &
                                        ,  DUM,  DUM,  3.0,  2.0,  2.0 &
                                        ,  DUM,  DUM,  4.0,  2.0,  1.0 &
                                        ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                        ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                        ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                        ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                        ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                        ,  DUM,  DUM,  DUM,  DUM,  DUM ], shape = [11, 5], order = [2, 1])
        matB = [ integer(TKG) ::  DUM,  DUM,  3.0,  2.0,  1.0 ]
        matC = [ integer(TKG) ::  4.0,  DUM,  5.0,  DUM,  2.0,  DUM,  3.0 ]
        refC = [ integer(TKG) :: 14.0,  DUM, 19.0,  DUM, 17.0,  DUM, 20.0 ]
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = iform )
        call disp%show("matB")
        call disp%show( matB , format = iform )
        call disp%show("matC")
        call disp%show( matC , format = iform )
        call disp%show("alpha = 1._TKG; beta = 1._TKG; nrow = 4; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 2;")
                        alpha = 1._TKG; beta = 1._TKG; nrow = 4; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 2;
        call disp%skip()
        call disp%show("call setMatMulAdd(matA, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC) ! full contiguous interface.")
                        call setMatMulAdd(matA, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
        call disp%show("matC")
        call disp%show( matC , format = iform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = iform )
        call disp%skip()
        call disp%show("matC = [ integer(TKG) ::  4.0,  DUM,  5.0,  DUM,  2.0,  DUM,  3.0 ]")
                        matC = [ integer(TKG) ::  4.0,  DUM,  5.0,  DUM,  2.0,  DUM,  3.0 ]
        call disp%show("call setMatMulAdd(matA(2:5, 3:5), matB(3:), matC(1::incC)) ! simplified interface.")
                        call setMatMulAdd(matA(2:5, 3:5), matB(3:), matC(1::incC))
        call disp%show("matC")
        call disp%show( matC , format = iform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = iform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! BLAS 2 - GEMV: General matrix-vector multiplications - complex")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        use pm_kind, only: TKG => RKS
        complex(TKG) :: alpha, beta
        complex(TKG), parameter :: COMPLEXDUM = cmplx(huge(0._TKG), huge(0._TKG), TKG)
        complex(TKG), allocatable :: matA(:,:), matB(:), matC(:), refC(:)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape([ complex(TKG)  ::  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  (1.0, 2.0), (3.0, 5.0),  (2.0, 0.0) &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  (2.0, 3.0), (7.0, 9.0),  (4.0, 8.0) &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  (7.0, 4.0), (1.0, 4.0),  (6.0, 0.0) &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  (8.0, 2.0), (2.0, 5.0),  (8.0, 0.0) &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  (9.0, 1.0), (3.0, 6.0),  (1.0, 0.0) &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM ], shape = [11, 5], order = [2, 1])
        matB = [ complex(TKG) ::  COMPLEXDUM,  COMPLEXDUM, (1.0, 2.0), (4.0, 0.0), (1.0, 1.0) ]
        matC = [ complex(TKG) ::  (1.0, 2.0), (4.0, 0.0), (1.0, -1.0), (3.0, 4.0), (2.0, 0.0) ]
        refC = [ complex(TKG) ::  (12.0, 28.0), (24.0, 55.0), (10.0, 39.0), (23.0, 50.0), (22.0, 44.0) ]
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = cform )
        call disp%show("matB")
        call disp%show( matB , format = cform )
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("alpha = 1._TKG; beta = 1._TKG; nrow = 5; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 1;")
                        alpha = 1._TKG; beta = 1._TKG; nrow = 5; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 1;
        call disp%skip()
        call disp%show("call setMatMulAdd(matA, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC) ! full contiguous interface.")
                        call setMatMulAdd(matA, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = cform )
        call disp%skip()
        call disp%show("matC = [ complex(TKG) ::  (1.0, 2.0), (4.0, 0.0), (1.0, -1.0), (3.0, 4.0), (2.0, 0.0) ]")
                        matC = [ complex(TKG) ::  (1.0, 2.0), (4.0, 0.0), (1.0, -1.0), (3.0, 4.0), (2.0, 0.0) ]
        call disp%show("call setMatMulAdd(matA(2:5, 3:5), matB(3:), matC(1::incC)) ! simplified interface.")
                        call setMatMulAdd(matA(2:6, 3:5), matB(3:), matC(1::incC))
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = cform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape([ complex(TKG)  ::  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  (1.0, 2.0), (3.0, 5.0),  (2.0, 0.0) &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  (2.0, 3.0), (7.0, 9.0),  (4.0, 8.0) &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  (7.0, 4.0), (1.0, 4.0),  (6.0, 0.0) &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  (8.0, 2.0), (2.0, 5.0),  (8.0, 0.0) &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  (9.0, 1.0), (3.0, 6.0),  (1.0, 0.0) &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM ], shape = [11, 5], order = [2, 1])
        matB = [ complex(TKG) ::  COMPLEXDUM,  COMPLEXDUM, (1.0, 2.0), (4.0, 0.0), (1.0, 1.0), (3.0, 4.0), (2.0, 0.0) ]
        matC = [ complex(TKG) ::  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM ]
        refC = [ complex(TKG) ::  (42.0, 67.0), (10.0, 87.0), (50.0, 74.0) ]
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = cform )
        call disp%show("matB")
        call disp%show( matB , format = cform )
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("alpha = 1._TKG; beta = 0._TKG; nrow = 5; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 1;")
                        alpha = 1._TKG; beta = 0._TKG; nrow = 5; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 1;
        call disp%skip()
        call disp%show("call setMatMulAdd(matA, transSymm, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC) ! full contiguous interface.")
                        call setMatMulAdd(matA, transSymm, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = cform )
        call disp%skip()
        call disp%show("matC = [ complex(TKG) ::  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM ]")
                        matC = [ complex(TKG) ::  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM ]
        call disp%show("call setMatMulAdd(matA(2:5, 3:5), transSymm, matB(3:), matC(1::incC), beta = beta) ! simplified interface.")
                        call setMatMulAdd(matA(2:6, 3:5), transSymm, matB(3:), matC(1::incC), beta = beta)
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = cform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape([ complex(TKG)  ::  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  (1.0, 2.0), (3.0, 5.0),  (2.0, 0.0) &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  (2.0, 3.0), (7.0, 9.0),  (4.0, 8.0) &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  (7.0, 4.0), (1.0, 4.0),  (6.0, 0.0) &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  (8.0, 2.0), (2.0, 5.0),  (8.0, 0.0) &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  (9.0, 1.0), (3.0, 6.0),  (1.0, 0.0) &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM &
                                        ,  COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM ], shape = [11, 5], order = [2, 1])
        matB = [ complex(TKG) ::  COMPLEXDUM,  COMPLEXDUM, (1.0, 2.0), (4.0, 0.0), (1.0, 1.0), (3.0, 4.0), (2.0, 0.0) ]
        matC = [ complex(TKG) ::  (1.0, 2.0), (4.0, 0.0), (1.0, -1.0) ]
        refC = [ complex(TKG) ::  (-73.0, -13.0), (-74.0, 57.0), (-49.0, -11.0) ]
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = cform )
        call disp%show("matB")
        call disp%show( matB , format = cform )
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("alpha = -1._TKG; beta = 1._TKG; nrow = 5; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 1;")
                        alpha = -1._TKG; beta = 1._TKG; nrow = 5; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 1;
        call disp%skip()
        call disp%show("call setMatMulAdd(matA, transHerm, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC) ! full contiguous interface.")
                        call setMatMulAdd(matA, transHerm, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = cform )
        call disp%skip()
        call disp%show("matC = [ complex(TKG) ::  (1.0, 2.0), (4.0, 0.0), (1.0, -1.0) ]")
                        matC = [ complex(TKG) ::  (1.0, 2.0), (4.0, 0.0), (1.0, -1.0) ]
        call disp%show("call setMatMulAdd(matA(2:5, 3:5), transHerm, matB(3:), matC(1::incC), alpha) ! simplified interface.")
                        call setMatMulAdd(matA(2:6, 3:5), transHerm, matB(3:), matC(1::incC), alpha)
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = cform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! BLAS 2 - GEMV: General matrix-vector multiplications -    real")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        use pm_kind, only: TKG => RKS
        real(TKG) :: alpha, beta
        real(TKG), parameter :: DUM = huge(DUM)
        real(TKG), allocatable :: matA(:,:), matB(:), matC(:), refC(:)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape([ real(TKG) ::  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  1.0,  2.0,  3.0 &
                                    ,  DUM,  DUM,  2.0,  2.0,  4.0 &
                                    ,  DUM,  DUM,  3.0,  2.0,  2.0 &
                                    ,  DUM,  DUM,  4.0,  2.0,  1.0 &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM ], shape = [11, 5], order = [2, 1])
        matB = [ real(TKG) ::  DUM,  DUM,  3.0,  2.0,  1.0 ]
        matC = [ real(TKG) ::  4.0,  DUM,  5.0,  DUM,  2.0, DUM,  3.0 ]
        refC = [ real(TKG) :: 14.0,  DUM, 19.0,  DUM, 17.0, DUM, 20.0 ]
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = rform )
        call disp%show("matB")
        call disp%show( matB , format = rform )
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("alpha = 1._TKG; beta = 1._TKG; nrow = 4; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 2;")
                        alpha = 1._TKG; beta = 1._TKG; nrow = 4; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 2;
        call disp%skip()
        call disp%show("call setMatMulAdd(matA, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC) ! full contiguous interface.")
                        call setMatMulAdd(matA, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()
        call disp%show("matC = [ real(TKG) ::  4.0,  DUM,  5.0,  DUM,  2.0, DUM,  3.0 ]")
                        matC = [ real(TKG) ::  4.0,  DUM,  5.0,  DUM,  2.0, DUM,  3.0 ]
        call disp%show("call setMatMulAdd(matA(2:5, 3:5), matB(3:), matC(1::incC)) ! simplified interface.")
                        call setMatMulAdd(matA(2:5, 3:5), matB(3:), matC(1::incC))
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape([ real(TKG) ::  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  1.0,  2.0,  3.0 &
                                    ,  DUM,  DUM,  2.0,  2.0,  4.0 &
                                    ,  DUM,  DUM,  3.0,  2.0,  2.0 &
                                    ,  DUM,  DUM,  4.0,  2.0,  1.0 &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM ], shape = [11, 5], order = [2, 1])
        matB = [ real(TKG) ::  DUM,  DUM,  3.0,  2.0,  1.0,  4.0 ]
        matC = [ real(TKG) ::  1.0,  DUM,  2.0,  DUM,  3.0 ]
        refC = [ real(TKG) :: 28.0,  DUM, 24.0,  DUM, 29.0 ]
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = rform )
        call disp%show("matB")
        call disp%show( matB , format = rform )
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("alpha = 1._TKG; beta = 2._TKG; nrow = 4; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 2;")
                        alpha = 1._TKG; beta = 2._TKG; nrow = 4; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 2;
        call disp%skip()
        call disp%show("call setMatMulAdd(matA, transSymm, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC) ! full contiguous interface.")
                        call setMatMulAdd(matA, transSymm, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()
        call disp%show("matC = [ real(TKG) ::  1.0,  DUM,  2.0,  DUM,  3.0 ]")
                        matC = [ real(TKG) ::  1.0,  DUM,  2.0,  DUM,  3.0 ]
        call disp%show("call setMatMulAdd(matA(2:5, 3:5), transSymm, matB(3:), matC(1::incC), beta = beta) ! simplified interface.")
                        call setMatMulAdd(matA(2:5, 3:5), transSymm, matB(3:), matC(1::incC), beta = beta)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape([ real(TKG) ::  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  1.0,  2.0,  3.0 &
                                    ,  DUM,  DUM,  2.0,  2.0,  4.0 &
                                    ,  DUM,  DUM,  3.0,  2.0,  2.0 &
                                    ,  DUM,  DUM,  4.0,  2.0,  1.0 &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM ], shape = [11, 5], order = [2, 1])
        matB = [ real(TKG) ::  DUM,  DUM,  3.0,  2.0,  1.0 ]
        matC = [ real(TKG) ::  4.0,  DUM,  5.0,  DUM,  2.0,  DUM,  3.0 ]
        refC = [ real(TKG) :: 14.0,  DUM, 19.0,  DUM, 17.0,  DUM, 20.0 ]
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = rform )
        call disp%show("matB")
        call disp%show( matB , format = rform )
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("alpha = 1._TKG; beta = 1._TKG; nrow = 4; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 2;")
                        alpha = 1._TKG; beta = 1._TKG; nrow = 4; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 2;
        call disp%skip()
        call disp%show("call setMatMulAdd(matA, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC) ! full contiguous interface.")
                        call setMatMulAdd(matA, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()
        call disp%show("matC = [ real(TKG) ::  4.0,  DUM,  5.0,  DUM,  2.0,  DUM,  3.0 ]")
                        matC = [ real(TKG) ::  4.0,  DUM,  5.0,  DUM,  2.0,  DUM,  3.0 ]
        call disp%show("call setMatMulAdd(matA(2:5, 3:5), matB(3:), matC(1::incC)) ! simplified interface.")
                        call setMatMulAdd(matA(2:5, 3:5), matB(3:), matC(1::incC))
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end block

    block

        use pm_kind, only: TKG => RKH
        real(TKG) :: alpha, beta
        real(TKG), parameter :: DUM = huge(DUM)
        real(TKG), allocatable :: matA(:,:), matB(:), matC(:), refC(:)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape([ real(TKG) ::  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  1.0,  2.0,  3.0 &
                                    ,  DUM,  DUM,  2.0,  2.0,  4.0 &
                                    ,  DUM,  DUM,  3.0,  2.0,  2.0 &
                                    ,  DUM,  DUM,  4.0,  2.0,  1.0 &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM ], shape = [11, 5], order = [2, 1])
        matB = [ real(TKG) ::  DUM,  DUM,  3.0,  2.0,  1.0 ]
        matC = [ real(TKG) ::  4.0,  DUM,  5.0,  DUM,  2.0, DUM,  3.0 ]
        refC = [ real(TKG) :: 14.0,  DUM, 19.0,  DUM, 17.0, DUM, 20.0 ]
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = rform )
        call disp%show("matB")
        call disp%show( matB , format = rform )
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("alpha = 1._TKG; beta = 1._TKG; nrow = 4; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 2;")
                        alpha = 1._TKG; beta = 1._TKG; nrow = 4; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 2;
        call disp%skip()
        call disp%show("call setMatMulAdd(matA, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC) ! full contiguous interface.")
                        call setMatMulAdd(matA, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()
        call disp%show("matC = [ real(TKG) ::  4.0,  DUM,  5.0,  DUM,  2.0, DUM,  3.0 ]")
                        matC = [ real(TKG) ::  4.0,  DUM,  5.0,  DUM,  2.0, DUM,  3.0 ]
        call disp%show("call setMatMulAdd(matA(2:5, 3:5), matB(3:), matC(1::incC)) ! simplified interface.")
                        call setMatMulAdd(matA(2:5, 3:5), matB(3:), matC(1::incC))
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape([ real(TKG) ::  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  1.0,  2.0,  3.0 &
                                    ,  DUM,  DUM,  2.0,  2.0,  4.0 &
                                    ,  DUM,  DUM,  3.0,  2.0,  2.0 &
                                    ,  DUM,  DUM,  4.0,  2.0,  1.0 &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM ], shape = [11, 5], order = [2, 1])
        matB = [ real(TKG) ::  DUM,  DUM,  3.0,  2.0,  1.0,  4.0 ]
        matC = [ real(TKG) ::  1.0,  DUM,  2.0,  DUM,  3.0 ]
        refC = [ real(TKG) :: 28.0,  DUM, 24.0,  DUM, 29.0 ]
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = rform )
        call disp%show("matB")
        call disp%show( matB , format = rform )
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("alpha = 1._TKG; beta = 2._TKG; nrow = 4; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 2;")
                        alpha = 1._TKG; beta = 2._TKG; nrow = 4; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 2;
        call disp%skip()
        call disp%show("call setMatMulAdd(matA, transSymm, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC) ! full contiguous interface.")
                        call setMatMulAdd(matA, transSymm, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()
        call disp%show("matC = [ real(TKG) ::  1.0,  DUM,  2.0,  DUM,  3.0 ]")
                        matC = [ real(TKG) ::  1.0,  DUM,  2.0,  DUM,  3.0 ]
        call disp%show("call setMatMulAdd(matA(2:5, 3:5), transSymm, matB(3:), matC(1::incC), beta = beta) ! simplified interface.")
                        call setMatMulAdd(matA(2:5, 3:5), transSymm, matB(3:), matC(1::incC), beta = beta)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape([ real(TKG) ::  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  1.0,  2.0,  3.0 &
                                    ,  DUM,  DUM,  2.0,  2.0,  4.0 &
                                    ,  DUM,  DUM,  3.0,  2.0,  2.0 &
                                    ,  DUM,  DUM,  4.0,  2.0,  1.0 &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM ], shape = [11, 5], order = [2, 1])
        matB = [ real(TKG) ::  DUM,  DUM,  3.0,  2.0,  1.0 ]
        matC = [ real(TKG) ::  4.0,  DUM,  5.0,  DUM,  2.0,  DUM,  3.0 ]
        refC = [ real(TKG) :: 14.0,  DUM, 19.0,  DUM, 17.0,  DUM, 20.0 ]
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = rform )
        call disp%show("matB")
        call disp%show( matB , format = rform )
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("alpha = 1._TKG; beta = 1._TKG; nrow = 4; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 2;")
                        alpha = 1._TKG; beta = 1._TKG; nrow = 4; ncol = 3; roffA = 1; coffA = 2; incB = 1; incC = 2;
        call disp%skip()
        call disp%show("call setMatMulAdd(matA, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC) ! full contiguous interface.")
                        call setMatMulAdd(matA, matB(3:), matC, alpha, beta, nrow, ncol, roffA, coffA, incB, incC)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()
        call disp%show("matC = [ real(TKG) ::  4.0,  DUM,  5.0,  DUM,  2.0,  DUM,  3.0 ]")
                        matC = [ real(TKG) ::  4.0,  DUM,  5.0,  DUM,  2.0,  DUM,  3.0 ]
        call disp%show("call setMatMulAdd(matA(2:5, 3:5), matB(3:), matC(1::incC)) ! simplified interface.")
                        call setMatMulAdd(matA(2:5, 3:5), matB(3:), matC(1::incC))
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Complete general integer matrix-matrix multiplications.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        use pm_kind, only: TKG => IK
        use pm_distUnif, only: setUnifRand
        integer(TKG), allocatable, dimension(:,:) :: matA, matB, matC, refC
        integer(TKG) :: alpha, beta
        integer(IK) :: nrow, ncol

        alpha = 2_TKG; beta = 3_TKG
        nrow = 2; ncol = 2; ndum = 3
        allocate(matA(nrow, ndum), matB(ndum, ncol), matC(nrow, ncol))

        call disp%skip()
        call disp%show("call setUnifRand(matA, lb = -10_TKG, ub = +10_TKG)")
                        call setUnifRand(matA, lb = -10_TKG, ub = +10_TKG)
        call disp%show("call setUnifRand(matB, lb = -10_TKG, ub = +10_TKG)")
                        call setUnifRand(matB, lb = -10_TKG, ub = +10_TKG)
        call disp%show("call setUnifRand(matC, lb = -10_TKG, ub = +10_TKG)")
                        call setUnifRand(matC, lb = -10_TKG, ub = +10_TKG)
        call disp%show("matA")
        call disp%show( matA , format = iform )
        call disp%show("matB")
        call disp%show( matB , format = iform )
        call disp%show("matC")
        call disp%show( matC , format = iform )
        call disp%show("[alpha, beta]")
        call disp%show( [alpha, beta] )
        call disp%show("[nrow, ncol, ndum]")
        call disp%show( [nrow, ncol, ndum] )
        call disp%show("refC = matmul(alpha * matA, matB) + beta * matC ! reference value.")
                        refC = matmul(alpha * matA, matB) + beta * matC
        call disp%show("call setMatMulAdd(matA, matB, matC, alpha, beta)")
                        call setMatMulAdd(matA, matB, matC, alpha, beta)
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = iform )
        call disp%skip()
        call disp%show("call setUnifRand(matC, lb = -10_TKG, ub = +10_TKG) ! reset for new multiplication.")
                        call setUnifRand(matC, lb = -10_TKG, ub = +10_TKG)
        call disp%show("matA = transpose(matA)")
                        matA = transpose(matA)
        call disp%show("matA")
        call disp%show( matA , format = iform )
        call disp%show("matC")
        call disp%show( matC , format = iform )
        call disp%show("refC = matmul(alpha * transpose(matA), matB) + beta * matC ! reference value.")
                        refC = matmul(alpha * transpose(matA), matB) + beta * matC
        call disp%show("call setMatMulAdd(matA, transSymm, matB, matC, alpha, beta)")
                        call setMatMulAdd(matA, transSymm, matB, matC, alpha, beta)
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = iform )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Complete general complex matrix-matrix multiplications.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        use pm_kind, only: TKG => CKS
        use pm_distUnif, only: setUnifRand
        complex(TKG), allocatable, dimension(:,:) :: matA, matB, matC, refC
        complex(TKG) :: alpha, beta
        integer(IK) :: nrow, ncol

        nrow = 2; ncol = 2; ndum = 3
        alpha = (1._TKG, 0._TKG); beta = (0._TKG, 0._TKG)
        allocate(matA(nrow, ndum), matB(ndum, ncol), matC(nrow, ncol))

        call disp%skip()
        call disp%show("call setUnifRand(matA)")
                        call setUnifRand(matA)
        call disp%show("call setUnifRand(matB)")
                        call setUnifRand(matB)
        call disp%show("call setUnifRand(matC)")
                        call setUnifRand(matC)
        call disp%show("matA")
        call disp%show( matA , format = cform )
        call disp%show("matB")
        call disp%show( matB , format = cform )
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("[alpha, beta]")
        call disp%show( [alpha, beta] )
        call disp%show("[nrow, ncol, ndum]")
        call disp%show( [nrow, ncol, ndum] )
        call disp%show("refC = matmul(alpha * matA, matB) + beta * matC ! reference value.")
                        refC = matmul(alpha * matA, matB) + beta * matC
        call disp%show("call setMatMulAdd(matA, matB, matC, alpha, beta)")
                        call setMatMulAdd(matA, matB, matC, alpha, beta)
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = cform )
        call disp%skip()
        call disp%show("call setUnifRand(matC) ! reset for new multiplication.")
                        call setUnifRand(matC)
        call disp%show("matA = conjg(transpose(matA))")
                        matA = conjg(transpose(matA))
        call disp%show("matA")
        call disp%show( matA , format = cform )
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("refC = matmul(alpha * conjg(transpose(matA)), matB) + beta * matC ! reference value.")
                        refC = matmul(alpha * conjg(transpose(matA)), matB) + beta * matC
        call disp%show("call setMatMulAdd(matA, transHerm, matB, matC, alpha, beta)")
                        call setMatMulAdd(matA, transHerm, matB, matC, alpha, beta)
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = cform )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Subset general integer matrix-matrix multiplication.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        use pm_kind, only: IKG => IK
        integer(TKG) :: alpha, beta
        integer(TKG), parameter :: DUM = huge(DUM)
        integer(TKG), allocatable, dimension(:,:) :: matA, matB, matC, refC

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape( [ integer(TKG) ::  1.0,  2.0,  -1.0,  -1.0,  4.0 &
                                        ,  2.0,  0.0,   1.0,   1.0, -1.0 &
                                        ,  1.0, -1.0,  -1.0,   1.0,  2.0 &
                                        , -3.0,  2.0,   2.0,   2.0,  0.0 &
                                        ,  4.0,  0.0,  -2.0,   1.0, -1.0 &
                                        , -1.0, -1.0,   1.0,  -3.0,  2.0 &
                                        ,  DUM,  DUM,   DUM,   DUM,  DUM &
                                        ,  DUM,  DUM,   DUM,   DUM,  DUM ], shape = [8, 5], order = [2, 1])
        matB = reshape( [ integer(TKG) ::  1.0, -1.0,   0.0,   2.0 &
                                        ,  2.0,  2.0,  -1.0,  -2.0 &
                                        ,  1.0,  0.0,  -1.0,   1.0 &
                                        , -3.0, -1.0,   1.0,  -1.0 &
                                        ,  4.0,  2.0,  -1.0,   1.0 &
                                        ,  DUM,  DUM,   DUM,   DUM ], shape = [6, 4], order = [2, 1])
        matC = reshape( [ integer(TKG) ::  1.0,  1.0,   1.0,   1.0 &
                                        ,  1.0,  1.0,   1.0,   1.0 &
                                        ,  1.0,  1.0,   1.0,   1.0 &
                                        ,  1.0,  1.0,   1.0,   1.0 &
                                        ,  1.0,  1.0,   1.0,   1.0 &
                                        ,  1.0,  1.0,   1.0,   1.0 &
                                        ,  DUM,  DUM,   DUM,   DUM ], shape = [6, 4], order = [2, 1])
        refC = reshape( [ integer(TKG) ::  24.0, 13.0, -5.0,  3.0 &
                                        ,  -3.0, -4.0,  2.0,  4.0 &
                                        ,   4.0,  1.0,  2.0,  5.0 &
                                        ,  -2.0,  6.0, -1.0, -9.0 &
                                        ,  -4.0, -6.0,  5.0,  5.0 &
                                        ,  16.0,  7.0, -4.0,  7.0 &
                                        ,   DUM,  DUM,  DUM,  DUM ], shape = [6, 4], order = [2, 1])
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = iform )
        call disp%show("matB")
        call disp%show( matB , format = iform )
        call disp%show("matC")
        call disp%show( matC , format = iform )
        call disp%show("alpha = 1._TKG; beta = 1._TKG; nrow = 6; ncol = 4; ndum = 5; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;")
                        alpha = 1._TKG; beta = 1._TKG; nrow = 6; ncol = 4; ndum = 5; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;
        call disp%show("call setMatMulAdd(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)")
                        call setMatMulAdd(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
        call disp%show("matC")
        call disp%show( matC , format = iform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = iform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape( [ integer(TKG) ::  1.0, -3.0 &
                                        ,  2.0,  4.0 &
                                        ,  1.0, -1.0 &
                                        ,  DUM,  DUM ], shape = [4, 2], order = [2, 1])
        matB = reshape( [ integer(TKG) ::  1.0, -3.0 &
                                        ,  2.0,  4.0 &
                                        ,  1.0, -1.0 ], shape = [3, 2], order = [2, 1])
        matC = reshape( [ integer(TKG) ::  1.0,  1.0,  1.0 &
                                        ,  1.0,  1.0,  1.0 &
                                        ,  1.0,  1.0,  1.0 &
                                        ,  DUM,  DUM,  DUM &
                                        ,  DUM,  DUM,  DUM ], shape = [5, 3], order = [2, 1])
        refC = reshape( [ integer(TKG) ::  11.0, -9.,  5.0 &
                                        ,  -9.0, 21., -1.0 &
                                        ,   5.0, -1.,  3.0 &
                                        ,   DUM, DUM,  DUM &
                                        ,   DUM, DUM,  DUM ], shape = [5, 3], order = [2, 1])
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = iform )
        call disp%show("matB")
        call disp%show( matB , format = iform )
        call disp%show("matC")
        call disp%show( matC , format = iform )
        call disp%show("alpha = 1._TKG; beta = 1._TKG; nrow = 3; ncol = 3; ndum = 2; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;")
                        alpha = 1._TKG; beta = 1._TKG; nrow = 3; ncol = 3; ndum = 2; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;
        call disp%show("call setMatMulAdd(matA, matB, transSymm, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)")
                        call setMatMulAdd(matA, matB, transSymm, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
        call disp%show("matC")
        call disp%show( matC , format = iform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = iform )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Subset general complex matrix-matrix multiplication.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        use pm_kind, only: TKG => CKS
        complex(TKG) :: alpha, beta
        complex(TKG), parameter :: COMPLEXDUM = cmplx(huge(0._TKG), huge(0._TKG), TKG)
        complex(TKG), allocatable, dimension(:,:) :: matA, matB, matC, refC

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape( [ complex(TKG) :: (1.0, 5.0), (9.0, 2.0), (1.0, 9.0) &
                                        , (2.0, 4.0), (8.0, 3.0), (1.0, 8.0) &
                                        , (3.0, 3.0), (7.0, 5.0), (1.0, 7.0) &
                                        , (4.0, 2.0), (4.0, 7.0), (1.0, 5.0) &
                                        , (5.0, 1.0), (5.0, 1.0), (1.0, 6.0) &
                                        , (6.0, 6.0), (3.0, 6.0), (1.0, 4.0) &
                                        , COMPLEXDUM, COMPLEXDUM, COMPLEXDUM &
                                        , COMPLEXDUM, COMPLEXDUM, COMPLEXDUM ], shape = [8, 3], order = [2, 1])
        matB = reshape( [ complex(TKG) :: (1.0, 8.0), (2.0, 7.0) &
                                        , (4.0, 4.0), (6.0, 8.0) &
                                        , (6.0, 2.0), (4.0, 5.0) &
                                        , COMPLEXDUM, COMPLEXDUM ], shape = [4, 2], order = [2, 1])
        matC = reshape( [ complex(TKG) :: (0.5, 0.0), (0.5, 0.0) &
                                        , (0.5, 0.0), (0.5, 0.0) &
                                        , (0.5, 0.0), (0.5, 0.0) &
                                        , (0.5, 0.0), (0.5, 0.0) &
                                        , (0.5, 0.0), (0.5, 0.0) &
                                        , (0.5, 0.0), (0.5, 0.0) &
                                        , COMPLEXDUM, COMPLEXDUM &
                                        , COMPLEXDUM, COMPLEXDUM ], shape = [8, 2], order = [2, 1])
        refC = reshape( [ complex(TKG) :: (-22.0, 113.0), (-35.0, 142.0) &
                                        , (-19.0, 114.0), (-35.0, 141.0) &
                                        , (-20.0, 119.0), (-43.0, 146.0) &
                                        , (-27.0, 110.0), (-58.0, 131.0) &
                                        ,   (8.0, 103.0),   (0.0, 112.0) &
                                        , (-55.0, 116.0), (-75.0, 135.0) &
                                        ,     COMPLEXDUM,     COMPLEXDUM &
                                        ,     COMPLEXDUM,     COMPLEXDUM ], shape = [8, 2], order = [2, 1])
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = cform )
        call disp%show("matB")
        call disp%show( matB , format = cform )
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("alpha = (1._TKG, 0._TKG); beta = (2._TKG, 0._TKG); nrow = 6; ncol = 2; ndum = 3; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;")
                        alpha = (1._TKG, 0._TKG); beta = (2._TKG, 0._TKG); nrow = 6; ncol = 2; ndum = 3; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;
        call disp%show("call setMatMulAdd(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)")
                        call setMatMulAdd(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = cform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape( [ complex(TKG) :: (1.0, 3.0), (-3.0, 2.0) &
                                        , (2.0, 5.0),  (4.0, 6.0) &
                                        , (1.0, 1.0), (-1.0, 9.0) ], shape = [3, 2], order = [2, 1])
        matB = reshape( [ complex(TKG) :: (1.0, 2.0), (-3.0, 2.0) &
                                        , (2.0, 6.0),  (4.0, 5.0) &
                                        , (1.0, 2.0), (-1.0, 8.0) &
                                        , COMPLEXDUM,  COMPLEXDUM ], shape = [4, 2], order = [2, 1])
        matC = reshape( [ complex(TKG) :: COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM &
                                        , COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM &
                                        , COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM &
                                        , COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM ], shape = [4, 3], order = [2, 1])
        refC = reshape( [ complex(TKG) :: (20.0,   1.0), (18.0, 23.0), (26.0,  23.0) &
                                        , (12.0, -25.0), (80.0,  2.0), (56.0, -37.0) &
                                        , (24.0, -26.0), (49.0, 37.0), (76.0,  -2.0) &
                                        ,    COMPLEXDUM,   COMPLEXDUM,    COMPLEXDUM ], shape = [4, 3], order = [2, 1])
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = cform )
        call disp%show("matB")
        call disp%show( matB , format = cform )
        call disp%show("matC ! Note that the initialization is irrelevant because `beta = (0., 0.)`.")
        call disp%show( matC , format = cform )
        call disp%show("alpha = (1._TKG, 0._TKG); beta = (0._TKG, 0._TKG); nrow = 3; ncol = 3; ndum = 2; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;")
                        alpha = (1._TKG, 0._TKG); beta = (0._TKG, 0._TKG); nrow = 3; ncol = 3; ndum = 2; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;
        call disp%show("call setMatMulAdd(matA, matB, transHerm, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)")
                        call setMatMulAdd(matA, matB, transHerm, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = cform )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Subset general    real matrix-matrix multiplication.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        real(TKG) :: alpha, beta
        real(TKG), parameter :: DUM = huge(DUM)
        real(TKG), allocatable, dimension(:,:) :: matA, matB, matC, refC

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape([ real(TKG) ::  1.0,  2.0,  -1.0,  -1.0,  4.0 &
                                    ,  2.0,  0.0,   1.0,   1.0, -1.0 &
                                    ,  1.0, -1.0,  -1.0,   1.0,  2.0 &
                                    , -3.0,  2.0,   2.0,   2.0,  0.0 &
                                    ,  4.0,  0.0,  -2.0,   1.0, -1.0 &
                                    , -1.0, -1.0,   1.0,  -3.0,  2.0 &
                                    ,  DUM,  DUM,   DUM,   DUM,  DUM &
                                    ,  DUM,  DUM,   DUM,   DUM,  DUM ], shape = [8, 5], order = [2, 1])
        matB = reshape([ real(TKG) ::  1.0, -1.0,   0.0,   2.0 &
                                    ,  2.0,  2.0,  -1.0,  -2.0 &
                                    ,  1.0,  0.0,  -1.0,   1.0 &
                                    , -3.0, -1.0,   1.0,  -1.0 &
                                    ,  4.0,  2.0,  -1.0,   1.0 &
                                    ,  DUM,  DUM,   DUM,   DUM ], shape = [6, 4], order = [2, 1])
        matC = reshape([ real(TKG) ::  0.5,  0.5,   0.5,   0.5 &
                                    ,  0.5,  0.5,   0.5,   0.5 &
                                    ,  0.5,  0.5,   0.5,   0.5 &
                                    ,  0.5,  0.5,   0.5,   0.5 &
                                    ,  0.5,  0.5,   0.5,   0.5 &
                                    ,  0.5,  0.5,   0.5,   0.5 &
                                    ,  DUM,  DUM,   DUM,   DUM ], shape = [6, 4], order = [2, 1])
        refC = reshape([ real(TKG) ::  24.0, 13.0, -5.0,  3.0 &
                                    ,  -3.0, -4.0,  2.0,  4.0 &
                                    ,   4.0,  1.0,  2.0,  5.0 &
                                    ,  -2.0,  6.0, -1.0, -9.0 &
                                    ,  -4.0, -6.0,  5.0,  5.0 &
                                    ,  16.0,  7.0, -4.0,  7.0 &
                                    ,   DUM,  DUM,  DUM,  DUM ], shape = [6, 4], order = [2, 1])
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = rform )
        call disp%show("matB")
        call disp%show( matB , format = rform )
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("alpha = 1._TKG; beta = 2._TKG; nrow = 6; ncol = 4; ndum = 5; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;")
                        alpha = 1._TKG; beta = 2._TKG; nrow = 6; ncol = 4; ndum = 5; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;
        call disp%show("call setMatMulAdd(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)")
                        call setMatMulAdd(matA, matB, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape([ real(TKG) ::  1.0, -3.0 &
                                    ,  2.0,  4.0 &
                                    ,  1.0, -1.0 &
                                    ,  DUM,  DUM ], shape = [4, 2], order = [2, 1])
        matB = reshape([ real(TKG) ::  1.0, -3.0 &
                                    ,  2.0,  4.0 &
                                    ,  1.0, -1.0 ], shape = [3, 2], order = [2, 1])
        matC = reshape([ real(TKG) ::  0.5,  0.5,  0.5 &
                                    ,  0.5,  0.5,  0.5 &
                                    ,  0.5,  0.5,  0.5 &
                                    ,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM ], shape = [5, 3], order = [2, 1])
        refC = reshape([ real(TKG) ::  11.0, -9.,  5.0 &
                                    ,  -9.0, 21., -1.0 &
                                    ,   5.0, -1.,  3.0 &
                                    ,   DUM, DUM,  DUM &
                                    ,   DUM, DUM,  DUM ], shape = [5, 3], order = [2, 1])
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = rform )
        call disp%show("matB")
        call disp%show( matB , format = rform )
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("alpha = 1._TKG; beta = 2._TKG; nrow = 3; ncol = 3; ndum = 2; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;")
                        alpha = 1._TKG; beta = 2._TKG; nrow = 3; ncol = 3; ndum = 2; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;
        call disp%show("call setMatMulAdd(matA, matB, transSymm, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)")
                        call setMatMulAdd(matA, matB, transSymm, matC, alpha, beta, nrow, ncol, ndum, roffA, coffA, roffB, coffB, roffC, coffC)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Subset symmetric complex matrix-matrix multiplication.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        use pm_kind, only: TKG => CKS
        complex(TKG) :: alpha, beta
        real(TKG), parameter :: DUM = huge(DUM)
        complex(TKG), parameter :: COMPLEXDUM = cmplx(huge(0._TKG), huge(0._TKG), TKG)
        complex(TKG), allocatable, dimension(:,:) :: matA, matB, matC, refC

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matB = reshape( [ complex(TKG) :: (1.0, 5.0), (-3.0, 2.0),  (1.0, 6.0) &
                                        , COMPLEXDUM,  (4.0, 5.0), (-1.0, 4.0) &
                                        , COMPLEXDUM,  COMPLEXDUM,  (2.0, 5.0) &
                                        , COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM &
                                        ], shape = [4, 3], order = [2, 1])
        matA = reshape( [ complex(TKG) :: (1.0, 1.0), (-3.0, 2.0),  (3.0, 3.0) &
                                        , (2.0, 6.0),  (4.0, 5.0), (-1.0, 4.0) &
                                        , COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM &
                                        ], shape = [3, 3], order = [2, 1])
        matC = reshape( [ complex(TKG) ::  (13.0, 6.0), (-18.0, 6.0), (10.0, 7.0) &
                                        , (-11.0, 8.0),  (11.0, 1.0), (-4.0, 2.0) &
                                        ,   COMPLEXDUM,   COMPLEXDUM,  COMPLEXDUM &
                                        ,   COMPLEXDUM,   COMPLEXDUM,  COMPLEXDUM &
                                        ,   COMPLEXDUM,   COMPLEXDUM,  COMPLEXDUM &
                                        ], shape = [5, 3], order = [2, 1])
        refC = reshape( [ complex(TKG) ::  (-96.0,   72.0), (-141.0, -226.0), (-112.0,   38.0) &
                                        , (-230.0, -269.0), (-133.0,  -23.0), (-272.0, -198.0) &
                                        ,       COMPLEXDUM,       COMPLEXDUM,       COMPLEXDUM &
                                        ,       COMPLEXDUM,       COMPLEXDUM,       COMPLEXDUM &
                                        ,       COMPLEXDUM,       COMPLEXDUM,       COMPLEXDUM &
                                        ], shape = [5, 3], order = [2, 1])
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = cform )
        call disp%show("matB")
        call disp%show( matB , format = cform )
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("alpha = (2._TKG, 3._TKG); beta = (1._TKG, 6._TKG); nrow = 2; ncol = 3; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;")
                        alpha = (2._TKG, 3._TKG); beta = (1._TKG, 6._TKG); nrow = 2; ncol = 3; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;
        call disp%show("call setMatMulAdd(matA, matB, symmetric, uppDia, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)")
                        call setMatMulAdd(matA, matB, symmetric, uppDia, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = cform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matB = reshape( [ complex(TKG) ::  (1.0, DUM), COMPLEXDUM, COMPLEXDUM &
                                        ,  (3.0, 2.0), (4.0, DUM), COMPLEXDUM &
                                        , (-1.0, 6.0), (1.0, 4.0), (2.0, DUM) &
                                        ,  COMPLEXDUM, COMPLEXDUM, COMPLEXDUM &
                                        ], shape = [4, 3], order = [2, 1])
        matA = reshape( [ complex(TKG) :: (1.0, 1.0), (-3.0, 2.0),  (3.0, 3.0) &
                                        , (2.0, 6.0),  (4.0, 5.0), (-1.0, 4.0) &
                                        , COMPLEXDUM,  COMPLEXDUM,  COMPLEXDUM &
                                        ], shape = [3, 3], order = [2, 1])
        matC = reshape( [ complex(TKG) ::  (13.0, 6.0), (-18.0, 6.0), (10.0, 7.0) &
                                        , (-11.0, 8.0),  (11.0, 1.0), (-4.0, 2.0) &
                                        ,   COMPLEXDUM,   COMPLEXDUM,  COMPLEXDUM &
                                        ,   COMPLEXDUM,   COMPLEXDUM,  COMPLEXDUM &
                                        ,   COMPLEXDUM,   COMPLEXDUM,  COMPLEXDUM &
                                        ], shape = [5, 3], order = [2, 1])
        refC = reshape( [ complex(TKG) :: (-137.0,  17.0), (-158.0, -102.0), (-39.0, 141.0) &
                                        , (-154.0, -77.0),  (-63.0,  186.0), (159.0, 104.0) &
                                        ,      COMPLEXDUM,       COMPLEXDUM,     COMPLEXDUM &
                                        ,      COMPLEXDUM,       COMPLEXDUM,     COMPLEXDUM &
                                        ,      COMPLEXDUM,       COMPLEXDUM,     COMPLEXDUM &
                                        ], shape = [5, 3], order = [2, 1])
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = cform )
        call disp%show("matB")
        call disp%show( matB , format = cform )
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("alpha = (2._TKG, 3._TKG); beta = (1._TKG, 6._TKG); nrow = 2; ncol = 3; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;")
                        alpha = (2._TKG, 3._TKG); beta = (1._TKG, 6._TKG); nrow = 2; ncol = 3; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;
        call disp%show("call setMatMulAdd(matA, matB, hermitian, lowDia, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)")
                        call setMatMulAdd(matA, matB, hermitian, lowDia, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = cform )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Subset symmetric    real matrix-matrix multiplication.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        real(TKG) :: alpha, beta
        real(TKG), parameter :: DUM = huge(DUM)
        real(TKG), allocatable, dimension(:,:) :: matA, matB, matC, refC

        matA = reshape([ real(TKG) ::  1.0,  2.0, -1.0, -1.0,  4.0 &
                                    ,  DUM,  0.0,  1.0,  1.0, -1.0 &
                                    ,  DUM,  DUM, -1.0,  1.0,  2.0 &
                                    ,  DUM,  DUM,  DUM,  2.0,  0.0 &
                                    ,  DUM,  DUM,  DUM,  DUM, -1.0 &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ,  DUM,  DUM,  DUM,  DUM,  DUM &
                                    ], shape = [8, 5], order = [2, 1])
        matB = reshape([ real(TKG) ::  1.0, -1.0,  0.0,  2.0 &
                                    ,  2.0,  2.0, -1.0, -2.0 &
                                    ,  1.0,  0.0, -1.0,  1.0 &
                                    , -3.0, -1.0,  1.0, -1.0 &
                                    ,  4.0,  2.0, -1.0,  1.0 &
                                    ,  DUM,  DUM,  DUM,  DUM &
                                    ], shape = [6, 4], order = [2, 1])
        matC = reshape([ real(TKG) :: 23.0, 12.0, -6.0,  2.0 &
                                    , -4.0, -5.0,  1.0,  3.0 &
                                    ,  5.0,  6.0, -1.0, -4.0 &
                                    , -4.0,  1.0,  0.0, -5.0 &
                                    ,  8.0, -4.0, -2.0, 13.0 &
                                    ], shape = [5, 4], order = [2, 1])
        refC = reshape([ real(TKG) ::  69.0,  36.0, -18.0,   6.0 &
                                    , -12.0, -15.0,   3.0,   9.0 &
                                    ,  15.0,  18.0,  -3.0, -12.0 &
                                    , -12.0,   3.0,   0.0, -15.0 &
                                    ,   8.0, -20.0,  -2.0,  35.0 &
                                    ], shape = [5, 4], order = [2, 1])
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = rform )
        call disp%show("matB")
        call disp%show( matB , format = rform )
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("alpha = 2._TKG; beta = 1._TKG; nrow = 5; ncol = 4; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;")
                        alpha = 2._TKG; beta = 1._TKG; nrow = 5; ncol = 4; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;
        call disp%show("call setMatMulAdd(matA, symmetric, uppDia, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)")
                        call setMatMulAdd(matA, symmetric, uppDia, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

        matC = reshape([ real(TKG) :: 23.0, 12.0, -6.0,  2.0 &
                                    , -4.0, -5.0,  1.0,  3.0 &
                                    ,  5.0,  6.0, -1.0, -4.0 &
                                    , -4.0,  1.0,  0.0, -5.0 &
                                    ,  8.0, -4.0, -2.0, 13.0 &
                                    ], shape = [5, 4], order = [2, 1])
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("call setMatMulAdd(matA(1:5, 1:5), symmetric, uppDia, matB(1:5, 1:4), matC, alpha, beta)")
                        call setMatMulAdd(matA(1:5, 1:5), symmetric, uppDia, matB(1:5, 1:4), matC, alpha, beta)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape([ real(TKG) :: 1.0,  DUM,  DUM &
                                    , 2.0,  4.0,  DUM &
                                    , 1.0, -1.0, -1.0 &
                                    , DUM,  DUM,  DUM &
                                    ], shape = [4, 3], order = [2, 1])
        matB = reshape([ real(TKG) :: 1.0, -3.0,  2.0,  2.0, -1.0,  2.0 &
                                    , 2.0,  4.0,  0.0,  0.0,  1.0, -2.0 &
                                    , 1.0, -1.0, -1.0, -1.0, -1.0,  1.0 &
                                    ], shape = [3, 6], order = [2, 1])
        matC = reshape([ real(TKG) ::  6.0,  4.0, 1.0, 1.0,  0.0, -1.0 &
                                    ,  9.0, 11.0, 5.0, 5.0,  3.0, -5.0 &
                                    , -2.0, -6.0, 3.0, 3.0, -1.0, 32.0 &
                                    ,  DUM,  DUM, DUM, DUM,  DUM,  DUM &
                                    ,  DUM,  DUM, DUM, DUM,  DUM,  DUM &
                                    ], shape = [5, 6], order = [2, 1])
        refC = reshape([ real(TKG) :: 24.0,  16.0,  4.0,  4.0,  0.0,  -4.0 &
                                    , 36.0,  44.0, 20.0, 20.0, 12.0, -20.0 &
                                    , -8.0, -24.0, 12.0, 12.0, -4.0,  70.0 &
                                    ,  DUM,   DUM,  DUM,  DUM,  DUM,   DUM &
                                    ,  DUM,   DUM,  DUM,  DUM,  DUM,   DUM &
                                    ], shape = [5, 6], order = [2, 1])
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = rform )
        call disp%show("matB")
        call disp%show( matB , format = rform )
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("alpha = 2._TKG; beta = 2._TKG; nrow = 3; ncol = 6; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;")
                        alpha = 2._TKG; beta = 2._TKG; nrow = 3; ncol = 6; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;
        call disp%show("call setMatMulAdd(matA, symmetric, lowDia, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)")
                        call setMatMulAdd(matA, symmetric, lowDia, matB, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

        matC = reshape([ real(TKG) ::  6.0,  4.0, 1.0, 1.0,  0.0, -1.0 &
                                    ,  9.0, 11.0, 5.0, 5.0,  3.0, -5.0 &
                                    , -2.0, -6.0, 3.0, 3.0, -1.0, 32.0 &
                                    ,  DUM,  DUM, DUM, DUM,  DUM,  DUM &
                                    ,  DUM,  DUM, DUM, DUM,  DUM,  DUM &
                                    ], shape = [5, 6], order = [2, 1])
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("call setMatMulAdd(matA(1:3, 1:3), symmetric, lowDia, matB(1:3, 1:6), matC(1:3, 1:6), alpha, beta)")
                        call setMatMulAdd(matA(1:3, 1:3), symmetric, lowDia, matB(1:3, 1:6), matC(1:3, 1:6), alpha, beta)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matB = reshape([ real(TKG) :: 1.0, -3.0,  1.0 &
                                    , DUM,  4.0, -1.0 &
                                    , DUM,  DUM,  2.0 &
                                    , DUM,  DUM,  DUM &
                                    ], shape = [4, 3], order = [2, 1])
        matA = reshape([ real(TKG) :: 1.0, -3.0,  3.0 &
                                    , 2.0,  4.0, -1.0 &
                                    , DUM,  DUM,  DUM &
                                    ], shape = [3, 3], order = [2, 1])
        matC = reshape([ real(TKG) ::  13.0, -18.0, 10.0 &
                                    , -11.0,  11.0, -4.0 &
                                    ,   DUM,   DUM,  DUM &
                                    ,   DUM,   DUM,  DUM &
                                    ,   DUM,   DUM,  DUM &
                                    ], shape = [5, 3], order = [2, 1])
        refC = reshape([ real(TKG) ::  39.0, -54.0,  30.0 &
                                    , -33.0,  33.0, -12.0 &
                                    ,   DUM,   DUM,   DUM &
                                    ,   DUM,   DUM,   DUM &
                                    ,   DUM,   DUM,   DUM &
                                    ], shape = [5, 3], order = [2, 1])
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = rform )
        call disp%show("matB")
        call disp%show( matB , format = rform )
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("alpha = 2._TKG; beta = 1._TKG; nrow = 2; ncol = 3; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;")
                        alpha = 2._TKG; beta = 1._TKG; nrow = 2; ncol = 3; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;
        call disp%show("call setMatMulAdd(matA, matB, symmetric, uppDia, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)")
                        call setMatMulAdd(matA, matB, symmetric, uppDia, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

        matC = reshape([ real(TKG) ::  13.0, -18.0, 10.0 &
                                    , -11.0,  11.0, -4.0 &
                                    ,   DUM,   DUM,  DUM &
                                    ,   DUM,   DUM,  DUM &
                                    ,   DUM,   DUM,  DUM &
                                    ], shape = [5, 3], order = [2, 1])
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("call setMatMulAdd(matA(1:2, 1:3), matB(1:3, 1:3), symmetric, uppDia, matC(1:2, 1:3), alpha, beta)")
                        call setMatMulAdd(matA(1:2, 1:3), matB(1:3, 1:3), symmetric, uppDia, matC(1:2, 1:3), alpha, beta)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matB = reshape([ real(TKG) :: 1.0,  DUM, DUM &
                                    , 2.0, 10.0, DUM &
                                    , 1.0, 11.0, 4.0 &
                                    ], shape = [3, 3], order = [2, 1])
        matA = reshape([ real(TKG) :: 1.0, -3.0,  2.0 &
                                    , 2.0,  4.0,  0.0 &
                                    , 1.0, -1.0, -1.0 &
                                    ], shape = [3, 3], order = [2, 1])
        matC = reshape([ real(TKG) ::  1.0,  5.0, -9.0 &
                                    , -3.0, 10.0, -2.0 &
                                    , -2.0,  8.0,  0.0 &
                                    ], shape = [3, 3], order = [2, 1])
        refC = reshape([ real(TKG) ::   4.0,  11.0,  15.0 &
                                    , -13.0, -34.0, -48.0 &
                                    ,   0.0,  27.0,  14.0 &
                                    ], shape = [3, 3], order = [2, 1])
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = rform )
        call disp%show("matB")
        call disp%show( matB , format = rform )
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("alpha = -1._TKG; beta = 1._TKG; nrow = 3; ncol = 3; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;")
                        alpha = -1._TKG; beta = 1._TKG; nrow = 3; ncol = 3; roffA = 0; coffA = 0; roffB = 0; coffB = 0; roffC = 0; coffC = 0;
        call disp%show("call setMatMulAdd(matA, matB, symmetric, lowDia, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)")
                        call setMatMulAdd(matA, matB, symmetric, lowDia, matC, alpha, beta, nrow, ncol, roffA, coffA, roffB, coffB, roffC, coffC)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

        matC = reshape([ real(TKG) ::  1.0,  5.0, -9.0 &
                                    , -3.0, 10.0, -2.0 &
                                    , -2.0,  8.0,  0.0 &
                                    ], shape = [3, 3], order = [2, 1])
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("call setMatMulAdd(matA(1:3, 1:3), matB(1:3, 1:3), symmetric, lowDia, matC(1:3, 1:3), alpha, beta)")
                        call setMatMulAdd(matA(1:3, 1:3), matB(1:3, 1:3), symmetric, lowDia, matC(1:3, 1:3), alpha, beta)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Subset symmetric    real matrix-vector multiplication.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        real(TKG)   :: alpha, beta
        real(TKG)   , parameter :: DUM = huge(DUM)
        real(TKG)   , allocatable :: matA(:,:), matB(:), matC(:), refC(:)
        integer(IK) :: ndim, incB, incC

        matA = reshape([ real(TKG) :: 8.0, DUM, DUM &
                                    , 4.0, 6.0, DUM &
                                    , 2.0, 7.0, 3.0 &
                                    ], shape = [3, 3], order = [2, 1])
        matB = [3.00, 2.0, 1.00]
        matC = [5.00, DUM, 3.00, DUM, 2.00]
        refC = [39.0, DUM, 34.0, DUM, 25.0]
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = rform )
        call disp%show("matB")
        call disp%show( matB , format = rform )
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("alpha = 1._TKG; beta = 1._TKG; ndim = 3; roffA = 0; coffA = 0; incB = 1; incC = 2;")
                        alpha = 1._TKG; beta = 1._TKG; ndim = 3; roffA = 0; coffA = 0; incB = 1; incC = 2;
        call disp%show("call setMatMulAdd(matA, symmetric, lowDia, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)")
                        call setMatMulAdd(matA, symmetric, lowDia, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

        matB = [3.00, 2.0, 1.00]
        matC = [5.00, 3.00, 2.00]
        refC = [39.0, 34.0, 25.0]
        call disp%show("matB")
        call disp%show( matB , format = rform )
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("alpha = 1._TKG; beta = 1._TKG;")
                        alpha = 1._TKG; beta = 1._TKG;
        call disp%show("call setMatMulAdd(matA, symmetric, lowDia, matB, matC, alpha, beta)")
                        call setMatMulAdd(matA, symmetric, lowDia, matB, matC, alpha, beta)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape([ real(TKG) :: 8.0, 4.0, 2.0 &
                                    , DUM, 6.0, 7.0 &
                                    , DUM, DUM, 3.0 &
                                    , DUM, DUM, DUM &
                                    ], shape = [3, 3], order = [2, 1])
        matB = [4.0, DUM, 2.0, DUM, 1.0]
        matC = [6.0, 5.0, 4.0]
        refC = [36.0, 54.0, 36.0]
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = rform )
        call disp%show("matB")
        call disp%show( matB , format = rform )
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("alpha = 1._TKG; beta = 2._TKG; ndim = 3; roffA = 0; coffA = 0; incB = -2; incC = 1;")
                        alpha = 1._TKG; beta = 2._TKG; ndim = 3; roffA = 0; coffA = 0; incB = -2; incC = 1;
        call disp%show("call setMatMulAdd(matA, symmetric, uppDia, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)")
                        call setMatMulAdd(matA, symmetric, uppDia, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

        matB = matB(size(matB):1:incB)
        matC = [6.0, 5.0, 4.0]
        call disp%show("matB")
        call disp%show( matB , format = rform )
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("alpha = 1._TKG; beta = 2._TKG;")
                        alpha = 1._TKG; beta = 2._TKG;
        call disp%show("call setMatMulAdd(matA, symmetric, uppDia, matB, matC, alpha, beta)")
                        call setMatMulAdd(matA, symmetric, uppDia, matB, matC, alpha, beta)
        call disp%show("matC")
        call disp%show( matC , format = rform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = rform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Subset Hermitian complex matrix-vector multiplication.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        use pm_kind, only: TKG => CKS
        complex(TKG)    :: alpha, beta
        real(TKG)       , parameter :: DUM = huge(DUM)
        complex(TKG)    , parameter :: COMPLEXDUM = cmplx(huge(0._TKG), huge(0._TKG), TKG)
        complex(TKG)    , allocatable :: matA(:,:), matB(:), matC(:), refC(:)
        integer(IK)     :: ndim, incB, incC

        matA = reshape([ real(TKG) :: (1.0,  DUM), COMPLEXDUM, COMPLEXDUM &
                                    , (3.0, -5.0), (7.0, DUM), COMPLEXDUM &
                                    , (2.0,  3.0), (4.0, 8.0), (6.0, DUM) &
                                    ], shape = [3, 3], order = [2, 1])
        matB = [(1.0, 2.0), (4.0, 0.0), (3.0, 4.0)]
        matC = [(1.0, 0.0), COMPLEXDUM, (2.0, -1.0), COMPLEXDUM, (2.0, 1.0)]
        refC = [(20., 10.), COMPLEXDUM, (45., 21.0), COMPLEXDUM, (38., 29.)]
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = cform )
        call disp%show("matB")
        call disp%show( matB , format = cform )
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("alpha = (1._TKG, 0._TKG); beta = (1._TKG, 0._TKG); ndim = 3; roffA = 0; coffA = 0; incB = 1; incC = 2;")
                        alpha = (1._TKG, 0._TKG); beta = (1._TKG, 0._TKG); ndim = 3; roffA = 0; coffA = 0; incB = 1; incC = 2;
        call disp%show("call setMatMulAdd(matA, hermitian, lowDia, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)")
                        call setMatMulAdd(matA, hermitian, lowDia, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = cform )
        call disp%skip()

        matC = [(1.0, 0.0), (2.0, -1.0), (2.0, 1.0)]
        refC = refC(1::incC)
        call disp%show("matB")
        call disp%show( matB , format = cform )
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("call setMatMulAdd(matA, hermitian, lowDia, matB, matC, alpha, beta)")
                        call setMatMulAdd(matA, hermitian, lowDia, matB, matC, alpha, beta)
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = cform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        matA = reshape([ real(TKG) :: COMPLEXDUM, COMPLEXDUM, COMPLEXDUM, COMPLEXDUM,  COMPLEXDUM &
                                    , COMPLEXDUM, COMPLEXDUM, (1.0, DUM), (3.0, 5.0), (2.0, -3.0) &
                                    , COMPLEXDUM, COMPLEXDUM, COMPLEXDUM, (7.0, DUM), (4.0, -8.0) &
                                    , COMPLEXDUM, COMPLEXDUM, COMPLEXDUM, COMPLEXDUM, (6.0,  DUM) &
                                    ], shape = [4, 5], order = [2, 1])
        matB = [(3.0, 4.0), COMPLEXDUM, (4.0, 0.0), COMPLEXDUM, (1.0, 2.0)]
        matC = [COMPLEXDUM, COMPLEXDUM, COMPLEXDUM, COMPLEXDUM, COMPLEXDUM]
        refC = [(19., 10.), COMPLEXDUM, (43., 22.), COMPLEXDUM, (36., 28.)]
        call disp%skip()
        call disp%show("matA")
        call disp%show( matA , format = cform )
        call disp%show("matB")
        call disp%show( matB , format = cform )
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("alpha = (1._TKG, 0._TKG); beta = (0._TKG, 0._TKG); ndim = 3; roffA = 1; coffA = 2; incB = -2; incC = 2;")
                        alpha = (1._TKG, 0._TKG); beta = (0._TKG, 0._TKG); ndim = 3; roffA = 1; coffA = 2; incB = -2; incC = 2;
        call disp%show("call setMatMulAdd(matA, hermitian, uppDia, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)")
                        call setMatMulAdd(matA, hermitian, uppDia, matB, matC, alpha, beta, ndim, roffA, coffA, incB, incC)
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = cform )
        call disp%skip()

        matA = matA(2:, 3:)
        matB = matB(size(matB):1:incB)
        matC = [COMPLEXDUM, COMPLEXDUM, COMPLEXDUM]
        refC = refC(1::-incB)
        call disp%show("matA")
        call disp%show( matA , format = cform )
        call disp%show("matB")
        call disp%show( matB , format = cform )
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("call setMatMulAdd(matA, hermitian, uppDia, matB, matC, alpha, beta)")
                        call setMatMulAdd(matA, hermitian, uppDia, matB, matC, alpha, beta)
        call disp%show("matC")
        call disp%show( matC , format = cform )
        call disp%show("matC - refC")
        call disp%show( matC - refC , format = cform )
        call disp%skip()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end block

end program example