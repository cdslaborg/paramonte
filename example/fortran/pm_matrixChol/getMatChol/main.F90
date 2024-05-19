program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: TKG => RKS ! all processor type kinds are supported.
    use pm_matrixChol, only: getMatChol, lowDia, uppDia
    use pm_io, only: display_type
    use pm_io, only: getFormat

    implicit none

    type(display_type) :: disp

    character(:, SK), allocatable   :: cform, gform
    real(TKG)       , parameter     :: DUM = -huge(0._TKG)
    complex(TKG)    , parameter     :: CMPLX_DUMM = cmplx(-huge(0._TKG), -huge(0._TKG), TKG)
    cform = getFormat([CMPLX_DUMM], ed = SK_'f', signed = .true.)
    gform = getFormat([DUM], ed = SK_'f', signed = .true.)

    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Cholesky factorization of    real positive-definite square matrix update.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    block
        real(TKG), allocatable :: cholref(:,:), cholmat(:,:)

        cholmat = reshape([real(TKG)   :: 1.0, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM &
                                        , 1.0, 2.0, DUM, DUM, DUM, DUM, DUM, DUM, DUM &
                                        , 1.0, 2.0, 3.0, DUM, DUM, DUM, DUM, DUM, DUM &
                                        , 1.0, 2.0, 3.0, 4.0, DUM, DUM, DUM, DUM, DUM &
                                        , 1.0, 2.0, 3.0, 4.0, 5.0, DUM, DUM, DUM, DUM &
                                        , 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, DUM, DUM, DUM &
                                        , 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, DUM, DUM &
                                        , 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, DUM &
                                        , 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 &
                                        ], shape = [9, 9], order = [2, 1])
        cholref = reshape([real(TKG)       :: 1.0, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM &
                                            , 1.0, 1.0, DUM, DUM, DUM, DUM, DUM, DUM, DUM &
                                            , 1.0, 1.0, 1.0, DUM, DUM, DUM, DUM, DUM, DUM &
                                            , 1.0, 1.0, 1.0, 1.0, DUM, DUM, DUM, DUM, DUM &
                                            , 1.0, 1.0, 1.0, 1.0, 1.0, DUM, DUM, DUM, DUM &
                                            , 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, DUM, DUM, DUM &
                                            , 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, DUM, DUM &
                                            , 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, DUM &
                                            , 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 &
                                            ], shape = [9, 9], order = [2, 1])
        call disp%skip()
        call disp%show("cholmat")
        call disp%show( cholmat , format = gform )
        call disp%show("cholmat = getMatChol(cholmat, lowDia)")
                        cholmat = getMatChol(cholmat, lowDia)
        call disp%show("cholmat")
        call disp%show( cholmat , format = gform )
        call disp%show("cholref")
        call disp%show( cholref , format = gform )
        call disp%skip()


        cholmat = reshape([real(TKG)   :: DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM &
                                        , DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM &
                                        , DUM, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 &
                                        , DUM, DUM, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0 &
                                        , DUM, DUM, DUM, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0 &
                                        , DUM, DUM, DUM, DUM, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0 &
                                        , DUM, DUM, DUM, DUM, DUM, 5.0, 5.0, 5.0, 5.0, 5.0 &
                                        , DUM, DUM, DUM, DUM, DUM, DUM, 6.0, 6.0, 6.0, 6.0 &
                                        , DUM, DUM, DUM, DUM, DUM, DUM, DUM, 7.0, 7.0, 7.0 &
                                        , DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, 8.0, 8.0 &
                                        , DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, 9.0 &
                                        ], shape = [11, 10], order = [2, 1])
        cholref = reshape([real(TKG)       :: DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM &
                                            , DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM &
                                            , DUM, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 &
                                            , DUM, DUM, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 &
                                            , DUM, DUM, DUM, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 &
                                            , DUM, DUM, DUM, DUM, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 &
                                            , DUM, DUM, DUM, DUM, DUM, 1.0, 1.0, 1.0, 1.0, 1.0 &
                                            , DUM, DUM, DUM, DUM, DUM, DUM, 1.0, 1.0, 1.0, 1.0 &
                                            , DUM, DUM, DUM, DUM, DUM, DUM, DUM, 1.0, 1.0, 1.0 &
                                            , DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, 1.0, 1.0 &
                                            , DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, 1.0 &
                                            ], shape = [11, 10], order = [2, 1])
        call disp%skip()
        call disp%show("cholmat")
        call disp%show( cholmat , format = gform )
        call disp%show("cholmat = getMatChol(cholmat(3:,2:), uppDia)")
                        cholmat = getMatChol(cholmat(3:,2:), uppDia)
        call disp%show("cholmat")
        call disp%show( cholmat , format = gform )
        call disp%show("cholref")
        call disp%show( cholref , format = gform )
        call disp%skip()

    end block

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Cholesky factorization of complex positive-definite square matrix update.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    block

        complex(TKG), allocatable :: cholref(:,:), cholmat(:,:)

        cholmat = reshape( [complex(TKG)   :: CMPLX_DUMM, CMPLX_DUMM,   CMPLX_DUMM,  CMPLX_DUMM,  CMPLX_DUMM &
                                            , CMPLX_DUMM, CMPLX_DUMM, (25.0,  DUM),  CMPLX_DUMM,  CMPLX_DUMM &
                                            , CMPLX_DUMM, CMPLX_DUMM, (-5.0,  5.0), (51.0, DUM),  CMPLX_DUMM &
                                            , CMPLX_DUMM, CMPLX_DUMM, (10.0, -5.0), ( 4.0, 6.0), (71.0, DUM) &
                                            ], shape = [4, 5], order = [2, 1])
        cholref = reshape( [complex(TKG)       :: CMPLX_DUMM, CMPLX_DUMM,   CMPLX_DUMM, CMPLX_DUMM,  CMPLX_DUMM &
                                                , CMPLX_DUMM, CMPLX_DUMM, ( 5.0,  0.0), CMPLX_DUMM,  CMPLX_DUMM &
                                                , CMPLX_DUMM, CMPLX_DUMM, (-1.0,  1.0), (7.0, 0.0),  CMPLX_DUMM &
                                                , CMPLX_DUMM, CMPLX_DUMM, ( 2.0, -1.0), (1.0, 1.0),  (8.0, 0.0) &
                                                ], shape = [4, 5], order = [2, 1])
        call disp%skip()
        call disp%show("cholmat")
        call disp%show( cholmat , format = cform )
        call disp%show("cholmat = getMatChol(cholmat(2:,3:), lowDia)")
                        cholmat = getMatChol(cholmat(2:,3:), lowDia)
        call disp%show("cholmat")
        call disp%show( cholmat , format = cform )
        call disp%show("cholref")
        call disp%show( cholref , format = cform )
        call disp%skip()

        cholmat = reshape( [complex(TKG)   :: CMPLX_DUMM, CMPLX_DUMM,   CMPLX_DUMM, CMPLX_DUMM &
                                            , (9.0, DUM),  ( 3.0,3.0), ( 3.0,-3.0), CMPLX_DUMM &
                                            , CMPLX_DUMM,  (18.0,DUM), ( 8.0,-6.0), CMPLX_DUMM &
                                            , CMPLX_DUMM,  CMPLX_DUMM, (43.0, DUM), CMPLX_DUMM &
                                            ], shape = [4, 4], order = [2, 1])
        cholref = reshape( [complex(TKG)   :: CMPLX_DUMM, CMPLX_DUMM,  CMPLX_DUMM, CMPLX_DUMM &
                                            , (3.0, 0.0), (1.0, 1.0), (1.0, -1.0), CMPLX_DUMM &
                                            , CMPLX_DUMM, (4.0, 0.0), (2.0, -1.0), CMPLX_DUMM &
                                            , CMPLX_DUMM, CMPLX_DUMM, (6.0,  0.0), CMPLX_DUMM &
                                            ], shape = [4, 4], order = [2, 1])
        call disp%skip()
        call disp%show("cholmat")
        call disp%show( cholmat , format = cform )
        call disp%show("cholmat = getMatChol(cholmat(2:,:3), uppDia)")
                        cholmat = getMatChol(cholmat(2:,:3), uppDia)
        call disp%show("cholmat")
        call disp%show( cholmat , format = cform )
        call disp%show("cholref")
        call disp%show( cholref , format = cform )
        call disp%skip()

    end block

end program example