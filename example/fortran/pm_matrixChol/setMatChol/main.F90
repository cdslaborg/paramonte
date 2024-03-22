program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: TKC => RKS ! all processor type kinds are supported.
    use pm_kind, only: CKC => CKS ! all processor type kinds are supported.
    use pm_matrixChol, only: setMatChol, lowDia, uppDia
    use pm_matrixChol, only: iteration, recursion
    use pm_matrixChol, only: nothing, transHerm
    use pm_arrayResize, only: setResized
    use pm_io, only: display_type
    use pm_io, only: getFormat

    implicit none

    type(display_type) :: disp

    character(:, SK), allocatable   :: cform, gform
    real(TKC)       , parameter     :: DUM = -huge(0._TKC)
    complex(CKC)    , parameter     :: CMPLX_DUMM = cmplx(-huge(0._CKC), -huge(0._CKC), CKC)
    integer(IK)                     :: info, ndim, roff, coff
    cform = getFormat([CMPLX_DUMM], ed = SK_'f', signed = .true.)
    gform = getFormat([DUM], ed = SK_'f', signed = .true.)

    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKC => RKS
        real(TKC), allocatable :: mat(:,:), cholow(:,:), choupp(:,:)
        mat = reshape(  [ 1._TKC, 0._TKC, 2._TKC &
                        , 0._TKC, 4._TKC, 0._TKC &
                        , 2._TKC, 0._TKC, 8._TKC &
                        ], shape = [3,3], order = [2, 1])
        call disp%skip
        call disp%show("mat")
        call disp%show( mat )
        call setResized(choupp, shape(mat, IK))
        call disp%show("choupp = 0")
                        choupp = 0
        call disp%show("call setMatChol(mat, uppDia, info, choupp, nothing)")
                        call setMatChol(mat, uppDia, info, choupp, nothing)
        call disp%show("if (info /= 0) error stop")
                        if (info /= 0) error stop
        call disp%show("choupp")
        call disp%show( choupp )
        call setResized(cholow, shape(mat, IK))
        call disp%show("cholow = 0")
                        cholow = 0
        call disp%show("call setMatChol(mat, uppDia, info, cholow, transHerm)")
                        call setMatChol(mat, uppDia, info, cholow, transHerm)
        call disp%show("if (info /= 0) error stop")
                        if (info /= 0) error stop
        call disp%show("cholow")
        call disp%show( cholow )
        call disp%show("matmul(cholow, choupp) - mat ! must be all zero-valued elements.")
        call disp%show( matmul(cholow, choupp) - mat )
        call disp%skip
    end block

    block
        use pm_kind, only: TKC => RKS
        complex(TKC), allocatable :: mat(:,:), cholow(:,:), choupp(:,:)
        mat = reshape(  [ (9.0,  0.0), (3.0, 3.0), (3.0, -3.0) &
                        , (3.0, -3.0),(18.0, 0.0), (8.0, -6.0) &
                        , (3.0,  3.0), (8.0, 6.0),(43.0,  0.0) &
                        ], shape = [3,3], order = [2, 1])
        call disp%skip
        call disp%show("mat")
        call disp%show( mat )
        call setResized(choupp, shape(mat, IK))
        call disp%show("choupp = 0")
                        choupp = 0
        call disp%show("call setMatChol(mat, uppDia, info, choupp, nothing)")
                        call setMatChol(mat, uppDia, info, choupp, nothing)
        call disp%show("if (info /= 0) error stop")
                        if (info /= 0) error stop
        call disp%show("choupp")
        call disp%show( choupp )
        call setResized(cholow, shape(mat, IK))
        call disp%show("cholow = 0")
                        cholow = 0
        call disp%show("call setMatChol(mat, uppDia, info, cholow, transHerm)")
                        call setMatChol(mat, uppDia, info, cholow, transHerm)
        call disp%show("if (info /= 0) error stop")
                        if (info /= 0) error stop
        call disp%show("cholow")
        call disp%show( cholow )
        call disp%show("matmul(cholow, choupp) - mat ! must be all zero-valued elements.")
        call disp%show( matmul(cholow, choupp) - mat )
        call disp%skip
    end block

    block
        use pm_kind, only: TKC => RKS
        complex(TKC), allocatable :: mat(:,:), cholow(:,:), choupp(:,:)
        mat = reshape(  [  (25.0, 0.0), (-5.0, -5.0), (10.0, 5.0) &
                        ,  (-5.0, 5.0),  (51.0, 0.0), (4.0, -6.0) &
                        , (10.0, -5.0),   (4.0, 6.0), (71.0, 0.0) &
                        ], shape = [3,3], order = [2, 1])
        call disp%skip
        call disp%show("mat")
        call disp%show( mat )
        call setResized(choupp, shape(mat, IK))
        call disp%show("choupp = 0")
                        choupp = 0
        call disp%show("call setMatChol(mat, uppDia, info, choupp, nothing)")
                        call setMatChol(mat, uppDia, info, choupp, nothing)
        call disp%show("if (info /= 0) error stop")
                        if (info /= 0) error stop
        call disp%show("choupp")
        call disp%show( choupp )
        call setResized(cholow, shape(mat, IK))
        call disp%show("cholow = 0")
                        cholow = 0
        call disp%show("call setMatChol(mat, uppDia, info, cholow, transHerm)")
                        call setMatChol(mat, uppDia, info, cholow, transHerm)
        call disp%show("if (info /= 0) error stop")
                        if (info /= 0) error stop
        call disp%show("cholow")
        call disp%show( cholow )
        call disp%show("matmul(cholow, choupp) - mat ! must be all zero-valued elements.")
        call disp%show( matmul(cholow, choupp) - mat )
        call disp%skip
    end block

    block
        use pm_kind, only: TKC => RKS
        real(TKC), allocatable :: mat(:,:), choupp(:,:), cholow(:,:)
        mat = reshape(  [ 1._TKC, 0._TKC, 2._TKC &
                        , 0._TKC, 4._TKC, 0._TKC &
                        , 2._TKC, 0._TKC, 8._TKC &
                        ], shape = [3,3], order = [2, 1])
        call disp%skip
        call disp%show("mat")
        call disp%show( mat )
        call setResized(cholow, shape(mat, IK))
        call disp%show("cholow = 0")
                        cholow = 0
        call disp%show("call setMatChol(mat, lowDia, info, cholow, nothing)")
                        call setMatChol(mat, lowDia, info, cholow, nothing)
        call disp%show("if (info /= 0) error stop")
                        if (info /= 0) error stop
        call disp%show("cholow")
        call disp%show( cholow )
        call setResized(choupp, shape(mat, IK))
        call disp%show("choupp = 0")
        call disp%show("choupp = 0")
                        choupp = 0
        call disp%show("call setMatChol(mat, lowDia, info, choupp, transHerm)")
                        call setMatChol(mat, lowDia, info, choupp, transHerm)
        call disp%show("if (info /= 0) error stop")
                        if (info /= 0) error stop
        call disp%show("choupp")
        call disp%show( choupp )
        call disp%show("matmul(cholow, choupp) - mat ! must be all zero-valued elements.")
        call disp%show( matmul(cholow, choupp) - mat )
        call disp%skip
    end block

    block
        use pm_kind, only: TKC => RKS
        complex(TKC), allocatable :: mat(:,:), choupp(:,:), cholow(:,:)
        mat = reshape(  [ (9.0,  0.0), (3.0, 3.0), (3.0, -3.0) &
                        , (3.0, -3.0),(18.0, 0.0), (8.0, -6.0) &
                        , (3.0,  3.0), (8.0, 6.0),(43.0,  0.0) &
                        ], shape = [3,3], order = [2, 1])
        call disp%skip
        call disp%show("mat")
        call disp%show( mat )
        call setResized(cholow, shape(mat, IK))
        call disp%show("cholow = 0")
                        cholow = 0
        call disp%show("call setMatChol(mat, lowDia, info, cholow, nothing)")
                        call setMatChol(mat, lowDia, info, cholow, nothing)
        call disp%show("if (info /= 0) error stop")
                        if (info /= 0) error stop
        call disp%show("cholow")
        call disp%show( cholow )
        call setResized(choupp, shape(mat, IK))
        call disp%show("choupp = 0")
        call disp%show("choupp = 0")
                        choupp = 0
        call disp%show("call setMatChol(mat, lowDia, info, choupp, transHerm)")
                        call setMatChol(mat, lowDia, info, choupp, transHerm)
        call disp%show("if (info /= 0) error stop")
                        if (info /= 0) error stop
        call disp%show("choupp")
        call disp%show( choupp )
        call disp%show("matmul(cholow, choupp) - mat ! must be all zero-valued elements.")
        call disp%show( matmul(cholow, choupp) - mat )
        call disp%skip
    end block

    block
        use pm_kind, only: TKC => RKS
        complex(TKC), allocatable :: mat(:,:), choupp(:,:), cholow(:,:)
        mat = reshape(  [  (25.0, 0.0), (-5.0, -5.0), (10.0, 5.0) &
                        ,  (-5.0, 5.0),  (51.0, 0.0), (4.0, -6.0) &
                        , (10.0, -5.0),   (4.0, 6.0), (71.0, 0.0) &
                        ], shape = [3,3], order = [2, 1])
        call disp%skip
        call disp%show("mat")
        call disp%show( mat )
        call setResized(cholow, shape(mat, IK))
        call disp%show("cholow = 0")
                        cholow = 0
        call disp%show("call setMatChol(mat, lowDia, info, cholow, nothing)")
                        call setMatChol(mat, lowDia, info, cholow, nothing)
        call disp%show("if (info /= 0) error stop")
                        if (info /= 0) error stop
        call disp%show("cholow")
        call disp%show( cholow )
        call setResized(choupp, shape(mat, IK))
        call disp%show("choupp = 0")
        call disp%show("choupp = 0")
                        choupp = 0
        call disp%show("call setMatChol(mat, lowDia, info, choupp, transHerm)")
                        call setMatChol(mat, lowDia, info, choupp, transHerm)
        call disp%show("if (info /= 0) error stop")
                        if (info /= 0) error stop
        call disp%show("choupp")
        call disp%show( choupp )
        call disp%show("matmul(cholow, choupp) - mat ! must be all zero-valued elements.")
        call disp%show( matmul(cholow, choupp) - mat )
        call disp%skip
    end block

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!In-place and out-of-place factorization within the same matrix")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    block
        integer(IK), parameter :: ndim = 9
        real(TKC) :: cholref(ndim, ndim), mat(ndim, ndim), matchol(ndim, ndim + 1)
        mat = reshape([real(TKC)    ::1.0, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM &
                                    , 1.0, 2.0, DUM, DUM, DUM, DUM, DUM, DUM, DUM &
                                    , 1.0, 2.0, 3.0, DUM, DUM, DUM, DUM, DUM, DUM &
                                    , 1.0, 2.0, 3.0, 4.0, DUM, DUM, DUM, DUM, DUM &
                                    , 1.0, 2.0, 3.0, 4.0, 5.0, DUM, DUM, DUM, DUM &
                                    , 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, DUM, DUM, DUM &
                                    , 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, DUM, DUM &
                                    , 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, DUM &
                                    , 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 &
                                    ], shape = [ndim, ndim], order = [2, 1])
        cholref = reshape([real(TKC)::1.0, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM &
                                    , 1.0, 1.0, DUM, DUM, DUM, DUM, DUM, DUM, DUM &
                                    , 1.0, 1.0, 1.0, DUM, DUM, DUM, DUM, DUM, DUM &
                                    , 1.0, 1.0, 1.0, 1.0, DUM, DUM, DUM, DUM, DUM &
                                    , 1.0, 1.0, 1.0, 1.0, 1.0, DUM, DUM, DUM, DUM &
                                    , 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, DUM, DUM, DUM &
                                    , 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, DUM, DUM &
                                    , 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, DUM &
                                    , 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 &
                                    ], shape = [ndim, ndim], order = [2, 1])
        call disp%skip()
        call disp%show("mat")
        call disp%show( mat , format = gform )
        matchol = DUM
        call disp%show("matchol(:,1:ndim) = mat")
                        matchol(:,1:ndim) = mat
        call disp%show("matchol(:,1:ndim)")
        call disp%show( matchol(:,1:ndim) , format = gform )
        call disp%show("call setMatChol(matchol(:,1:ndim), lowDia, info, matchol(:,1:ndim), nothing)")
                        call setMatChol(matchol(:,1:ndim), lowDia, info, matchol(:,1:ndim), nothing)
        call disp%show("matchol(:,1:ndim)")
        call disp%show( matchol(:,1:ndim) , format = gform )
        call disp%show("matchol(:,1:ndim) - cholref")
        call disp%show( matchol(:,1:ndim) - cholref , format = gform )
        call disp%show("if (info /= 0) error stop 'Cholesky factorization failed.'")
                        if (info /= 0) error stop 'Cholesky factorization failed.'
        call disp%skip()
        matchol = DUM
        call disp%show("matchol(:,1:ndim) = mat")
                        matchol(:,1:ndim) = mat
        call disp%show("matchol(:,1:ndim)")
        call disp%show( matchol(:,1:ndim) , format = gform )
        call disp%show("call setMatChol(matchol(:,1:ndim), lowDia, info, matchol(:,2:ndim+1), transHerm)")
                        call setMatChol(matchol(:,1:ndim), lowDia, info, matchol(:,2:ndim+1), transHerm)
        call disp%show("matchol(:,2:ndim+1)")
        call disp%show( matchol(:,2:ndim+1) , format = gform )
        call disp%show("matchol(:,2:ndim+1) - transpose(cholref)")
        call disp%show( matchol(:,2:ndim+1) - transpose(cholref) , format = gform )
        call disp%show("if (info /= 0) error stop 'Cholesky factorization failed.'")
                        if (info /= 0) error stop 'Cholesky factorization failed.'
        call disp%skip()

    end block

    block
        integer(IK), parameter :: ndim = 3
        complex(TKC) :: cholref(ndim, ndim), mat(ndim, ndim), matchol(ndim, ndim + 1)
        mat = reshape( [complex(CKC)    ::(25.0,  0.0),  CMPLX_DUMM,  CMPLX_DUMM &
                                        , (-5.0,  5.0), (51.0, 0.0),  CMPLX_DUMM &
                                        , (10.0, -5.0), ( 4.0, 6.0), (71.0, 0.0) &
                                        ], shape = [ndim, ndim], order = [2, 1])
        cholref = reshape( [complex(CKC)::( 5.0,  0.0), CMPLX_DUMM,  CMPLX_DUMM &
                                        , (-1.0,  1.0), (7.0, 0.0),  CMPLX_DUMM &
                                        , ( 2.0, -1.0), (1.0, 1.0),  (8.0, 0.0) &
                                        ], shape = [ndim, ndim], order = [2, 1])
        call disp%skip()
        call disp%show("mat")
        call disp%show( mat , format = cform )
        matchol = CMPLX_DUMM
        call disp%show("matchol(:,1:ndim) = mat")
                        matchol(:,1:ndim) = mat
        call disp%show("matchol(:,1:ndim)")
        call disp%show( matchol(:,1:ndim) , format = cform )
        call disp%show("call setMatChol(matchol(:,1:ndim), lowDia, info, matchol(:,1:ndim), nothing)")
                        call setMatChol(matchol(:,1:ndim), lowDia, info, matchol(:,1:ndim), nothing)
        call disp%show("matchol(:,1:ndim)")
        call disp%show( matchol(:,1:ndim) , format = cform )
        call disp%show("cholref")
        call disp%show( cholref , format = cform )
        call disp%show("matchol(:,1:ndim) - cholref")
        call disp%show( matchol(:,1:ndim) - cholref , format = cform )
        call disp%show("if (info /= 0) error stop 'Cholesky factorization failed.'")
                        if (info /= 0) error stop 'Cholesky factorization failed.'
        call disp%skip()
        matchol = CMPLX_DUMM
        call disp%show("matchol(:,1:ndim) = mat")
                        matchol(:,1:ndim) = mat
        call disp%show("matchol(:,1:ndim)")
        call disp%show( matchol(:,1:ndim) , format = cform )
        call disp%show("call setMatChol(matchol(:,1:ndim), lowDia, info, matchol(:,2:ndim+1), transHerm)")
                        call setMatChol(matchol(:,1:ndim), lowDia, info, matchol(:,2:ndim+1), transHerm)
        call disp%show("matchol(:,2:ndim+1)")
        call disp%show( matchol(:,2:ndim+1) , format = cform )
        call disp%show("transpose(conjg(cholref))")
        call disp%show( transpose(conjg(cholref)) , format = cform )
        call disp%show("matchol(:,2:ndim+1) - transpose(conjg(cholref))")
        call disp%show( matchol(:,2:ndim+1) - transpose(conjg(cholref)) , format = cform )
        call disp%show("if (info /= 0) error stop 'Cholesky factorization failed.'")
                        if (info /= 0) error stop 'Cholesky factorization failed.'
        call disp%skip()

    end block

    call runExamples()
    call runExamples(bdim = 2_IK)

contains

    subroutine runExamples(bdim)
        integer(IK), optional :: bdim

        call disp%skip
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!Cholesky factorization of    real positive-definite square matrix update.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip

        block

            real(TKC), allocatable :: cholref(:,:), cholmat(:,:)

            cholmat = reshape([real(TKC)   :: 1.0, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM &
                                            , 1.0, 2.0, DUM, DUM, DUM, DUM, DUM, DUM, DUM &
                                            , 1.0, 2.0, 3.0, DUM, DUM, DUM, DUM, DUM, DUM &
                                            , 1.0, 2.0, 3.0, 4.0, DUM, DUM, DUM, DUM, DUM &
                                            , 1.0, 2.0, 3.0, 4.0, 5.0, DUM, DUM, DUM, DUM &
                                            , 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, DUM, DUM, DUM &
                                            , 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, DUM, DUM &
                                            , 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, DUM &
                                            , 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 &
                                            ], shape = [9, 9], order = [2, 1])
            cholref = reshape([real(TKC)       :: 1.0, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM &
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
            call disp%show("ndim = 9; roff = 0; coff = 0;")
                            ndim = 9; roff = 0; coff = 0;
            if (present(bdim)) then
                call disp%show("bdim")
                call disp%show( bdim )
                call disp%show("call setMatChol(cholmat, lowDia, info, iteration, ndim, roff, coff, bdim)")
                                call setMatChol(cholmat, lowDia, info, iteration, ndim, roff, coff, bdim)
            else
                call disp%show("call setMatChol(cholmat, lowDia, info, recursion, ndim, roff, coff)")
                                call setMatChol(cholmat, lowDia, info, recursion, ndim, roff, coff)
            end if
            call disp%show("cholmat")
            call disp%show( cholmat , format = gform )
            call disp%show("cholmat - cholref")
            call disp%show( cholmat - cholref , format = gform )
            call disp%show("info")
            call disp%show( info )
            call disp%skip()


            cholmat = reshape([real(TKC)   :: DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM &
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
            cholref = reshape([real(TKC)       :: DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM, DUM &
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
            call disp%show("ndim = 9; roff = 2; coff = 1;")
                            ndim = 9; roff = 2; coff = 1;
            if (present(bdim)) then
                call disp%show("bdim")
                call disp%show( bdim )
                call disp%show("call setMatChol(cholmat, uppDia, info, iteration, ndim, roff, coff, bdim)")
                                call setMatChol(cholmat, uppDia, info, iteration, ndim, roff, coff, bdim)
            else
                call disp%show("call setMatChol(cholmat, uppDia, info, recursion, ndim, roff, coff)")
                                call setMatChol(cholmat, uppDia, info, recursion, ndim, roff, coff)
            end if
            call disp%show("cholmat")
            call disp%show( cholmat , format = gform )
            call disp%show("cholmat - cholref")
            call disp%show( cholmat - cholref , format = gform )
            call disp%show("info")
            call disp%show( info )
            call disp%skip()

        end block

        call disp%skip
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%show("!Cholesky factorization of complex positive-definite square matrix update.")
        call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        call disp%skip

        block

            complex(CKC), allocatable :: cholref(:,:), cholmat(:,:)

            cholmat = reshape( [complex(CKC)   :: CMPLX_DUMM, CMPLX_DUMM,   CMPLX_DUMM,  CMPLX_DUMM,  CMPLX_DUMM &
                                                , CMPLX_DUMM, CMPLX_DUMM, (25.0,  DUM),  CMPLX_DUMM,  CMPLX_DUMM &
                                                , CMPLX_DUMM, CMPLX_DUMM, (-5.0,  5.0), (51.0, DUM),  CMPLX_DUMM &
                                                , CMPLX_DUMM, CMPLX_DUMM, (10.0, -5.0), ( 4.0, 6.0), (71.0, DUM) &
                                                ], shape = [4, 5], order = [2, 1])
            cholref = reshape( [complex(CKC)       :: CMPLX_DUMM, CMPLX_DUMM,   CMPLX_DUMM, CMPLX_DUMM,  CMPLX_DUMM &
                                                    , CMPLX_DUMM, CMPLX_DUMM, ( 5.0,  0.0), CMPLX_DUMM,  CMPLX_DUMM &
                                                    , CMPLX_DUMM, CMPLX_DUMM, (-1.0,  1.0), (7.0, 0.0),  CMPLX_DUMM &
                                                    , CMPLX_DUMM, CMPLX_DUMM, ( 2.0, -1.0), (1.0, 1.0),  (8.0, 0.0) &
                                                    ], shape = [4, 5], order = [2, 1])
            call disp%skip()
            call disp%show("cholmat")
            call disp%show( cholmat , format = cform )
            call disp%show("ndim = 3; roff = 1; coff = 2;")
                            ndim = 3; roff = 1; coff = 2;
            call disp%show("call setMatChol(cholmat, lowDia, info, iteration, ndim, roff, coff)")
                            call setMatChol(cholmat, lowDia, info, iteration, ndim, roff, coff)
            call disp%show("cholmat")
            call disp%show( cholmat , format = cform )
            call disp%show("cholmat - cholref")
            call disp%show( cholmat - cholref , format = cform )
            call disp%show("if (info /= 0) error stop 'Cholesky factorization failed.'")
                            if (info /= 0) error stop 'Cholesky factorization failed.'
            call disp%skip()

            cholmat = reshape( [complex(CKC)   :: CMPLX_DUMM, CMPLX_DUMM,   CMPLX_DUMM, CMPLX_DUMM &
                                                , (9.0, DUM),  ( 3.0,3.0), ( 3.0,-3.0), CMPLX_DUMM &
                                                , CMPLX_DUMM,  (18.0,DUM), ( 8.0,-6.0), CMPLX_DUMM &
                                                , CMPLX_DUMM,  CMPLX_DUMM, (43.0, DUM), CMPLX_DUMM &
                                                ], shape = [4, 4], order = [2, 1])
            cholref = reshape( [complex(CKC)       :: CMPLX_DUMM, CMPLX_DUMM,  CMPLX_DUMM, CMPLX_DUMM &
                                                    , (3.0, 0.0), (1.0, 1.0), (1.0, -1.0), CMPLX_DUMM &
                                                    , CMPLX_DUMM, (4.0, 0.0), (2.0, -1.0), CMPLX_DUMM &
                                                    , CMPLX_DUMM, CMPLX_DUMM, (6.0,  0.0), CMPLX_DUMM &
                                                    ], shape = [4, 4], order = [2, 1])
            call disp%skip()
            call disp%show("cholmat")
            call disp%show( cholmat , format = cform )
            call disp%show("ndim = 3; roff = 1; coff = 0;")
                            ndim = 3; roff = 1; coff = 0;
            if (present(bdim)) then
                call disp%show("bdim")
                call disp%show( bdim )
                call disp%show("call setMatChol(cholmat, uppDia, info, iteration, ndim, roff, coff, bdim)")
                                call setMatChol(cholmat, uppDia, info, iteration, ndim, roff, coff, bdim)
            else
                call disp%show("call setMatChol(cholmat, uppDia, info, recursion, ndim, roff, coff)")
                                call setMatChol(cholmat, uppDia, info, recursion, ndim, roff, coff)
            end if
            call disp%show("cholmat")
            call disp%show( cholmat , format = cform )
            call disp%show("cholmat - cholref")
            call disp%show( cholmat - cholref , format = cform )
            call disp%show("if (info /= 0) error stop 'Cholesky factorization failed.'")
                            if (info /= 0) error stop 'Cholesky factorization failed.'
            call disp%skip()

        end block

    end subroutine

end program example