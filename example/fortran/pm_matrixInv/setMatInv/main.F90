program example

    use pm_io, only: getFormat
    use pm_kind, only: SK, IK, LK, TKC => RKS
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRand
    use pm_matrixInv, only: setMatInv
    use pm_matrixInv, only: choUpp, choLow
    use pm_matrixInv, only: upperDiag, lowerDiag
    use pm_matrixInv, only: upperUnit, lowerUnit
    use pm_matrixChol, only: getMatChol, uppDia, lowDia
    use pm_matrixInit, only: getMatInit, uppLowDia
    use pm_distCov, only: getCovRand
    use pm_matrixCopy, only: rdpack, uppDia
    use pm_arrayResize, only: setResized
    use pm_matrixCopy, only: setMatCopy
    use pm_matrixCopy, only: transHerm
    use pm_matrixInit, only: setMatInit
    use pm_arrayFill, only: getFilled
    use pm_matrixLUP, only: setMatLUP
    use pm_err, only: setAsserted
    use pm_val2str, only: getStr

    implicit none

    integer(IK):: itry, info, ndim, ntry = 10
    character(:, SK), allocatable :: cform, rform

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    cform = getFormat(mold = [(0._TKC, 0._TKC)], ed = SK_"f", ndigit = 2_IK, signed = .true._LK)
    rform = getFormat(mold = [0._TKC], ed = SK_"f", ndigit = 2_IK, signed = .true._LK)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute inverse of an upperDiag triangular    real matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS
        real(TKC), allocatable :: mat(:,:), inv(:,:), mul(:,:)
        do itry = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getUnifRand(1_IK, 9_IK, ndim, ndim)")
                            mat = getUnifRand(1_IK, 9_IK, ndim, ndim)
            call disp%show("call setMatInit(mat(2 : ndim, 1 : ndim - 1), lowDia, 0._TKC, 0._TKC)")
                            call setMatInit(mat(2 : ndim, 1 : ndim - 1), lowDia, 0._TKC, 0._TKC)
            call disp%show("mat")
            call disp%show( mat , format = rform )
            call disp%show("call setResized(inv, shape(mat, IK))")
                            call setResized(inv, shape(mat, IK))
            call disp%show("call setMatInv(inv, mat, upperDiag)")
                            call setMatInv(inv, mat, upperDiag)
            call disp%show("inv")
            call disp%show( inv , format = rform )
            call disp%show("mul = matmul(mat, inv)")
                            mul = matmul(mat, inv)
            call disp%show("mul")
            call disp%show( mul , format = rform )
            call disp%skip
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute inverse of a upperDiag triangular complex matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS
        complex(TKC), allocatable :: mat(:,:), inv(:,:), mul(:,:)
        do itry = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getUnifRand((1., 1.), (2., 2.), ndim, ndim)")
                            mat = getUnifRand((1., 1.), (2., 2.), ndim, ndim)
            call disp%show("call setMatInit(mat(2 : ndim, 1 : ndim - 1), lowDia, (0._TKC, 0._TKC), (0._TKC, 0._TKC))")
                            call setMatInit(mat(2 : ndim, 1 : ndim - 1), lowDia, (0._TKC, 0._TKC), (0._TKC, 0._TKC))
            call disp%show("mat")
            call disp%show( mat , format = cform )
            call disp%show("call setResized(inv, shape(mat, IK))")
                            call setResized(inv, shape(mat, IK))
            call disp%show("call setMatInv(inv, mat, upperDiag)")
                            call setMatInv(inv, mat, upperDiag)
            call disp%show("inv")
            call disp%show( inv , format = cform )
            call disp%show("mul = matmul(mat, inv)")
                            mul = matmul(mat, inv)
            call disp%show("mul")
            call disp%show( mul , format = cform )
            call disp%skip
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute inverse of a lowerDiag triangular    real matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS
        real(TKC), allocatable :: mat(:,:), inv(:,:), mul(:,:)
        do itry = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getUnifRand(1_IK, 9_IK, ndim, ndim)")
                            mat = getUnifRand(1_IK, 9_IK, ndim, ndim)
            call disp%show("call setMatInit(mat(1 : ndim - 1, 2 : ndim), uppDia, 0._TKC, 0._TKC)")
                            call setMatInit(mat(1 : ndim - 1, 2 : ndim), uppDia, 0._TKC, 0._TKC)
            call disp%show("mat")
            call disp%show( mat , format = rform )
            call disp%show("call setResized(inv, shape(mat, IK))")
                            call setResized(inv, shape(mat, IK))
            call disp%show("call setMatInv(inv, mat, lowerDiag)")
                            call setMatInv(inv, mat, lowerDiag)
            call disp%show("inv")
            call disp%show( inv , format = rform )
            call disp%show("mul = matmul(mat, inv)")
                            mul = matmul(mat, inv)
            call disp%show("mul")
            call disp%show( mul , format = rform )
            call disp%skip
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute inverse of a lowerDiag triangular complex matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS
        complex(TKC), allocatable :: mat(:,:), inv(:,:), mul(:,:)
        do itry = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getUnifRand((1., 1.), (2., 2.), ndim, ndim)")
                            mat = getUnifRand((1., 1.), (2., 2.), ndim, ndim)
            call disp%show("call setMatInit(mat(1 : ndim - 1, 2 : ndim), uppDia, (0._TKC, 0._TKC), (0._TKC, 0._TKC))")
                            call setMatInit(mat(1 : ndim - 1, 2 : ndim), uppDia, (0._TKC, 0._TKC), (0._TKC, 0._TKC))
            call disp%show("mat")
            call disp%show( mat , format = cform )
            call disp%show("call setResized(inv, shape(mat, IK))")
                            call setResized(inv, shape(mat, IK))
            call disp%show("call setMatInv(inv, mat, lowerDiag)")
                            call setMatInv(inv, mat, lowerDiag)
            call disp%show("inv")
            call disp%show( inv , format = cform )
            call disp%show("mul = matmul(mat, inv)")
                            mul = matmul(mat, inv)
            call disp%show("mul")
            call disp%show( mul , format = cform )
            call disp%skip
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute inverse of a upperUnit triangular    real matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS
        real(TKC), allocatable :: mat(:,:), inv(:,:), mul(:,:)
        do itry = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getUnifRand(1_IK, 9_IK, ndim, ndim)")
                            mat = getUnifRand(1_IK, 9_IK, ndim, ndim)
            call disp%show("call setMatInit(mat, lowDia, 0._TKC, 1._TKC)")
                            call setMatInit(mat, lowDia, 0._TKC, 1._TKC)
            call disp%show("mat")
            call disp%show( mat , format = rform )
            call disp%show("call setResized(inv, shape(mat, IK))")
                            call setResized(inv, shape(mat, IK))
            call disp%show("call setMatInv(inv, mat, upperUnit)")
                            call setMatInv(inv, mat, upperUnit)
            call disp%show("inv")
            call disp%show( inv , format = rform )
            call disp%show("mul = matmul(mat, inv)")
                            mul = matmul(mat, inv)
            call disp%show("mul")
            call disp%show( mul , format = rform )
            call disp%skip
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute inverse of a upperUnit triangular complex matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS
        complex(TKC), allocatable :: mat(:,:), inv(:,:), mul(:,:)
        do itry = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getUnifRand((-1., -1.), (+1., +1.), ndim, ndim)")
                            mat = getUnifRand((-1., -1.), (+1., +1.), ndim, ndim)
            call disp%show("call setMatInit(mat, lowDia, (0._TKC, 0._TKC), (1._TKC, 0._TKC))")
                            call setMatInit(mat, lowDia, (0._TKC, 0._TKC), (1._TKC, 0._TKC))
            call disp%show("mat")
            call disp%show( mat , format = cform )
            call disp%show("call setResized(inv, shape(mat, IK))")
                            call setResized(inv, shape(mat, IK))
            call disp%show("call setMatInv(inv, mat, upperUnit)")
                            call setMatInv(inv, mat, upperUnit)
            call disp%show("inv")
            call disp%show( inv , format = cform )
            call disp%show("mul = matmul(mat, inv)")
                            mul = matmul(mat, inv)
            call disp%show("mul")
            call disp%show( mul , format = cform )
            call disp%skip
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute inverse of a lowerUnit triangular    real matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS
        real(TKC), allocatable :: mat(:,:), inv(:,:), mul(:,:)
        do itry = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getUnifRand(1_IK, 9_IK, ndim, ndim)")
                            mat = getUnifRand(1_IK, 9_IK, ndim, ndim)
            call disp%show("call setMatInit(mat, uppDia, 0._TKC, 1._TKC)")
                            call setMatInit(mat, uppDia, 0._TKC, 1._TKC)
            call disp%show("mat")
            call disp%show( mat , format = rform )
            call disp%show("call setResized(inv, shape(mat, IK))")
                            call setResized(inv, shape(mat, IK))
            call disp%show("call setMatInv(inv, mat, lowerUnit)")
                            call setMatInv(inv, mat, lowerUnit)
            call disp%show("inv")
            call disp%show( inv , format = rform )
            call disp%show("mul = matmul(mat, inv)")
                            mul = matmul(mat, inv)
            call disp%show("mul")
            call disp%show( mul , format = rform )
            call disp%skip
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute inverse of a lowerUnit triangular complex matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS
        complex(TKC), allocatable :: mat(:,:), inv(:,:), mul(:,:)
        do itry = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getUnifRand((-1., -1.), (+1., +1.), ndim, ndim)")
                            mat = getUnifRand((-1., -1.), (+1., +1.), ndim, ndim)
            call disp%show("call setMatInit(mat, uppDia, (0._TKC, 0._TKC), (1._TKC, 0._TKC))")
                            call setMatInit(mat, uppDia, (0._TKC, 0._TKC), (1._TKC, 0._TKC))
            call disp%show("mat")
            call disp%show( mat , format = cform )
            call disp%show("call setResized(inv, shape(mat, IK))")
                            call setResized(inv, shape(mat, IK))
            call disp%show("call setMatInv(inv, mat, lowerUnit)")
                            call setMatInv(inv, mat, lowerUnit)
            call disp%show("inv")
            call disp%show( inv , format = cform )
            call disp%show("mul = matmul(mat, inv)")
                            mul = matmul(mat, inv)
            call disp%show("mul")
            call disp%show( mul , format = cform )
            call disp%skip
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute inverse of a general    real matrix by passing its LUP factorization.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS
        integer(IK), allocatable :: rperm(:)
        real(TKC), allocatable :: mat(:,:), lup(:,:), inv(:,:), mul(:,:)
        do itry = 1, ntry
            call disp%skip
            call disp%show("mat = reshape([1, 0, 2, -1, 5, 0, 0, 3, -9], shape = [3,3], order = [2, 1])")
                            mat = reshape([1, 0, 2, -1, 5, 0, 0, 3, -9], shape = [3,3], order = [2, 1])
            call disp%show("mat")
            call disp%show( mat , format = rform )
            call disp%show("call setResized(rperm, size(mat, 1, IK))")
                            call setResized(rperm, size(mat, 1, IK))
            call disp%show("lup = mat")
                            lup = mat
            call disp%show("lup")
            call disp%show( lup , format = rform )
            call disp%show("call setMatLUP(lup, rperm, info) ! compute the LUP factorization of the matrix.")
                            call setMatLUP(lup, rperm, info)
            call disp%show("if (info /= 0) error stop 'LUP factorization failed.'")
                            if (info /= 0) error stop 'LUP factorization failed.'
            call disp%show("lup")
            call disp%show( lup , format = rform )
            call disp%show("call setResized(inv, shape(mat, IK))")
                            call setResized(inv, shape(mat, IK))
            call disp%show("call setMatInv(inv, lup, rperm) ! compute the inverse of the matrix by passing its LUP factorization.")
                            call setMatInv(inv, lup, rperm)
            call disp%show("inv")
            call disp%show( inv , format = rform )
            call disp%show("mul = matmul(mat, inv)")
                            mul = matmul(mat, inv)
            call disp%show("mul")
            call disp%show( mul , format = rform )
            call disp%skip
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute inverse of a general complex matrix by passing its LUP factorization.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS
        integer(IK), allocatable :: rperm(:)
        complex(TKC), parameter :: mat(*,*) = reshape( &
        [ (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0), (4.0,-1.0), (4.4,-1.0), (4.8,-1.0), (5.2,-1.0) &
        , (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0), (4.0,-1.0), (4.4,-1.0), (4.8,-1.0) &
        , (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0), (4.0,-1.0), (4.4,-1.0) &
        , (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0), (4.0,-1.0) &
        , (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0), (3.6,-1.0) &
        , (4.0, 1.0), (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0), (3.2,-1.0) &
        , (4.4, 1.0), (4.0, 1.0), (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0), (2.8,-1.0) &
        , (4.8, 1.0), (4.4, 1.0), (4.0, 1.0), (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0), (2.4,-1.0) &
        , (5.2, 1.0), (4.8, 1.0), (4.4, 1.0), (4.0, 1.0), (3.6, 1.0), (3.2, 1.0), (2.8, 1.0), (2.4, 1.0), (2.0, 1.0) &
        ], shape = [9, 9], order = [2, 1])
        complex(TKC), dimension(size(mat,1), size(mat,2)) :: inv, lup, mul
        call disp%skip
        call disp%show("call setResized(rperm, size(mat, 1, IK))")
                        call setResized(rperm, size(mat, 1, IK))
        call disp%show("lup = mat")
                        lup = mat
        call disp%show("lup")
        call disp%show( lup , format = cform )
        call disp%show("call setMatLUP(lup, rperm, info) ! compute the LUP factorization of the matrix.")
                        call setMatLUP(lup, rperm, info)
        call disp%show("if (info /= 0) error stop 'LUP factorization failed.'")
                        if (info /= 0) error stop 'LUP factorization failed.'
        call disp%show("lup")
        call disp%show( lup , format = cform )
        call disp%show("call setMatInv(inv, lup, rperm) ! compute the inverse of the matrix by passing its LUP factorization.")
                        call setMatInv(inv, lup, rperm)
        call disp%show("inv")
        call disp%show( inv , format = cform )
        call disp%show("mul = matmul(mat, inv)")
                        mul = matmul(mat, inv)
        call disp%show("mul")
        call disp%show( mul , format = cform )
        call disp%skip
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the inverse of a positive-definite    real matrix by passing its Cholesky factorization.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS
        real(TKC), allocatable :: mat(:,:), chol(:,:), inv(:,:), mul(:,:)
        do itry = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getCovRand(mold = 1._TKC, ndim = ndim)")
                            mat = getCovRand(mold = 1._TKC, ndim = ndim)
            call disp%show("mat")
            call disp%show( mat , format = rform )
            call disp%show("chol = getMatChol(mat, subset = lowDia)")
                            chol = getMatChol(mat, subset = lowDia)
            call disp%show("chol")
            call disp%show( chol , format = rform )
            call disp%show("call setResized(inv, shape(mat, IK))")
                            call setResized(inv, shape(mat, IK))
            call disp%show("call setMatInv(inv, chol, auxil = choLow)")
                            call setMatInv(inv, chol, auxil = choLow)
            call disp%show("inv")
            call disp%show( inv , format = rform )
            call disp%show("mul = matmul(mat, inv)")
                            mul = matmul(mat, inv)
            call disp%show("mul")
            call disp%show( mul , format = rform )
            call disp%skip
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the inverse of a positive-definite complex matrix by passing its Cholesky factorization.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS
        integer(IK), parameter :: ndim = 3_IK
        complex(TKC), allocatable :: chol(:,:), inv(:,:), mul(:,:)
        complex(TKC), parameter :: mat(*,*) = reshape(  [ (9.0,  0.0), (3.0, 3.0), (3.0, -3.0) &
                                                        , (3.0, -3.0),(18.0, 0.0), (8.0, -6.0) &
                                                        , (3.0,  3.0), (8.0, 6.0),(43.0,  0.0) &
                                                        ], shape = [ndim, ndim], order = [2, 1])
        call disp%skip
        call disp%show("mat")
        call disp%show( mat , format = cform )
        call disp%show("chol = getMatChol(mat, subset = lowDia)")
                        chol = getMatChol(mat, subset = lowDia)
        call disp%show("chol")
        call disp%show( chol , format = cform )
        call disp%show("call setResized(inv, shape(mat, IK))")
                        call setResized(inv, shape(mat, IK))
        call disp%show("call setMatInv(inv, chol, auxil = choLow)")
                        call setMatInv(inv, chol, auxil = choLow)
        call disp%show("inv")
        call disp%show( inv , format = cform )
        call disp%show("mul = matmul(mat, inv)")
                        mul = matmul(mat, inv)
        call disp%show("mul")
        call disp%show( mul , format = cform )
        call disp%skip
    end block

    block
        use pm_kind, only: TKC => RKS
        integer(IK), parameter :: ndim = 3_IK
        complex(TKC), allocatable :: chol(:,:), inv(:,:), mul(:,:)
        complex(TKC), parameter :: mat(*,*) = reshape(  [ (9.0,  0.0), (3.0, 3.0), (3.0, -3.0) &
                                                        , (3.0, -3.0),(18.0, 0.0), (8.0, -6.0) &
                                                        , (3.0,  3.0), (8.0, 6.0),(43.0,  0.0) &
                                                        ], shape = [ndim, ndim], order = [2, 1])
        call disp%skip
        call disp%show("mat")
        call disp%show( mat , format = cform )
        call disp%show("chol = getMatChol(mat, subset = uppDia)")
                        chol = getMatChol(mat, subset = uppDia)
        call disp%show("chol")
        call disp%show( chol , format = cform )
        call disp%show("call setResized(inv, shape(mat, IK))")
                        call setResized(inv, shape(mat, IK))
        call disp%show("call setMatInv(inv, chol, auxil = choUpp)")
                        call setMatInv(inv, chol, auxil = choUpp)
        call disp%show("inv")
        call disp%show( inv , format = cform )
        call disp%show("mul = matmul(mat, inv)")
                        mul = matmul(mat, inv)
        call disp%show("mul")
        call disp%show( mul , format = cform )
        call disp%skip
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute a subset of the inverse of a positive-definite    real matrix by passing its Cholesky factorization.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS
        real(TKC), allocatable :: mat(:,:), chol(:,:), inv(:,:), mul(:,:)
        do itry = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getCovRand(mold = 1._TKC, ndim = ndim)")
                            mat = getCovRand(mold = 1._TKC, ndim = ndim)
            call disp%show("mat")
            call disp%show( mat , format = rform )
            call disp%show("chol = getMatChol(mat, subset = lowDia)")
                            chol = getMatChol(mat, subset = lowDia)
            call disp%show("chol")
            call disp%show( chol , format = rform )
            call disp%show("call setResized(inv, shape(mat, IK))")
                            call setResized(inv, shape(mat, IK))
            call disp%skip
            call disp%show("inv = 0")
                            inv = 0
            call disp%show("call setMatInv(inv, chol, auxil = choLow, subset = uppDia)")
                            call setMatInv(inv, chol, auxil = choLow, subset = uppDia)
            call disp%show("inv")
            call disp%show( inv , format = rform )
            call disp%show("call setMatCopy(inv(2:,1:), rdpack, inv(1:,2:), rdpack, uppDia, transHerm) ! symmetrize `inv` for multiplication.")
                            call setMatCopy(inv(2:,1:), rdpack, inv(1:,2:), rdpack, uppDia, transHerm)
            call disp%show("inv")
            call disp%show( inv , format = rform )
            call disp%show("mul = matmul(mat, inv)")
                            mul = matmul(mat, inv)
            call disp%show("mul")
            call disp%show( mul , format = rform )
            call disp%skip
            call disp%show("inv = 0")
                            inv = 0
            call disp%show("call setMatInv(inv, chol, auxil = choLow, subset = lowDia)")
                            call setMatInv(inv, chol, auxil = choLow, subset = lowDia)
            call disp%show("inv")
            call disp%show( inv , format = rform )
            call disp%show("call setMatCopy(inv(1:,2:), rdpack, inv(2:,1:), rdpack, lowDia, transHerm) ! symmetrize `inv` for multiplication.")
                            call setMatCopy(inv(1:,2:), rdpack, inv(2:,1:), rdpack, lowDia, transHerm)
            call disp%show("inv")
            call disp%show( inv , format = rform )
            call disp%show("mul = matmul(mat, inv)")
                            mul = matmul(mat, inv)
            call disp%show("mul")
            call disp%show( mul , format = rform )
            call disp%skip
        end do
    end block

    block
        use pm_kind, only: TKC => RKS
        real(TKC), allocatable :: mat(:,:), chol(:,:), inv(:,:), mul(:,:)
        do itry = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getCovRand(mold = 1._TKC, ndim = ndim)")
                            mat = getCovRand(mold = 1._TKC, ndim = ndim)
            call disp%show("mat")
            call disp%show( mat , format = rform )
            call disp%show("chol = getMatChol(mat, subset = uppDia)")
                            chol = getMatChol(mat, subset = uppDia)
            call disp%show("chol")
            call disp%show( chol , format = rform )
            call disp%show("call setResized(inv, shape(mat, IK))")
                            call setResized(inv, shape(mat, IK))
            call disp%skip
            call disp%show("inv = 0")
                            inv = 0
            call disp%show("call setMatInv(inv, chol, auxil = choUpp, subset = uppDia)")
                            call setMatInv(inv, chol, auxil = choUpp, subset = uppDia)
            call disp%show("inv")
            call disp%show( inv , format = rform )
            call disp%show("call setMatCopy(inv(2:,1:), rdpack, inv(1:,2:), rdpack, uppDia, transHerm) ! symmetrize `inv` for multiplication.")
                            call setMatCopy(inv(2:,1:), rdpack, inv(1:,2:), rdpack, uppDia, transHerm)
            call disp%show("inv")
            call disp%show( inv , format = rform )
            call disp%show("mul = matmul(mat, inv)")
                            mul = matmul(mat, inv)
            call disp%show("mul")
            call disp%show( mul , format = rform )
            call disp%skip
            call disp%show("inv = 0")
                            inv = 0
            call disp%show("call setMatInv(inv, chol, auxil = choUpp, subset = lowDia)")
                            call setMatInv(inv, chol, auxil = choUpp, subset = lowDia)
            call disp%show("inv")
            call disp%show( inv , format = rform )
            call disp%show("call setMatCopy(inv(1:,2:), rdpack, inv(2:,1:), rdpack, lowDia, transHerm) ! symmetrize `inv` for multiplication.")
                            call setMatCopy(inv(1:,2:), rdpack, inv(2:,1:), rdpack, lowDia, transHerm)
            call disp%show("inv")
            call disp%show( inv , format = rform )
            call disp%show("mul = matmul(mat, inv)")
                            mul = matmul(mat, inv)
            call disp%show("mul")
            call disp%show( mul , format = rform )
            call disp%skip
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute a subset of the inverse of a positive-definite complex matrix by passing its Cholesky factorization.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS
        integer(IK), parameter :: ndim = 3_IK
        complex(TKC), allocatable :: chol(:,:), inv(:,:), mul(:,:)
        complex(TKC), parameter :: mat(*,*) = reshape(  [ (9.0,  0.0), (3.0, 3.0), (3.0, -3.0) &
                                                        , (3.0, -3.0),(18.0, 0.0), (8.0, -6.0) &
                                                        , (3.0,  3.0), (8.0, 6.0),(43.0,  0.0) &
                                                        ], shape = [ndim, ndim], order = [2, 1])
        call disp%skip
        call disp%show("mat")
        call disp%show( mat , format = cform )
        call disp%show("chol = getMatChol(mat, subset = lowDia)")
                        chol = getMatChol(mat, subset = lowDia)
        call disp%show("chol")
        call disp%show( chol , format = cform )
        call disp%show("call setResized(inv, shape(mat, IK))")
                        call setResized(inv, shape(mat, IK))
        call disp%show("inv = 0")
                        inv = 0
        call disp%show("call setMatInv(inv, chol, auxil = choLow, subset = uppDia)")
                        call setMatInv(inv, chol, auxil = choLow, subset = uppDia)
        call disp%show("inv")
        call disp%show( inv , format = cform )
        call disp%show("call setMatCopy(inv(2:,1:), rdpack, inv(1:,2:), rdpack, uppDia, transHerm) ! symmetrize `inv` for multiplication.")
                        call setMatCopy(inv(2:,1:), rdpack, inv(1:,2:), rdpack, uppDia, transHerm)
        call disp%show("inv")
        call disp%show( inv , format = cform )
        call disp%show("mul = matmul(mat, inv)")
                        mul = matmul(mat, inv)
        call disp%show("mul")
        call disp%show( mul , format = cform )
        call disp%skip
        call disp%show("inv = 0")
                        inv = 0
        call disp%show("call setMatInv(inv, chol, auxil = choLow, subset = lowDia)")
                        call setMatInv(inv, chol, auxil = choLow, subset = lowDia)
        call disp%show("inv")
        call disp%show( inv , format = cform )
        call disp%show("call setMatCopy(inv(1:,2:), rdpack, inv(2:,1:), rdpack, lowDia, transHerm) ! symmetrize `inv` for multiplication.")
                        call setMatCopy(inv(1:,2:), rdpack, inv(2:,1:), rdpack, lowDia, transHerm)
        call disp%show("inv")
        call disp%show( inv , format = cform )
        call disp%show("mul = matmul(mat, inv)")
                        mul = matmul(mat, inv)
        call disp%show("mul")
        call disp%show( mul , format = cform )
        call disp%skip
    end block

    block
        use pm_kind, only: TKC => RKS
        integer(IK), parameter :: ndim = 3_IK
        complex(TKC), allocatable :: chol(:,:), inv(:,:), mul(:,:)
        complex(TKC), parameter :: mat(*,*) = reshape(  [ (9.0,  0.0), (3.0, 3.0), (3.0, -3.0) &
                                                        , (3.0, -3.0),(18.0, 0.0), (8.0, -6.0) &
                                                        , (3.0,  3.0), (8.0, 6.0),(43.0,  0.0) &
                                                        ], shape = [ndim, ndim], order = [2, 1])
        call disp%skip
        call disp%show("mat")
        call disp%show( mat , format = cform )
        call disp%show("chol = getMatChol(mat, subset = uppDia)")
                        chol = getMatChol(mat, subset = uppDia)
        call disp%show("chol")
        call disp%show( chol , format = cform )
        call disp%show("call setResized(inv, shape(mat, IK))")
                        call setResized(inv, shape(mat, IK))
        call disp%show("inv = 0")
                        inv = 0
        call disp%show("call setMatInv(inv, chol, auxil = choUpp, subset = lowDia)")
                        call setMatInv(inv, chol, auxil = choUpp, subset = lowDia)
        call disp%show("inv")
        call disp%show( inv , format = cform )
        call disp%show("call setMatCopy(inv(1:,2:), rdpack, inv(2:,1:), rdpack, lowDia, transHerm) ! symmetrize `inv` for multiplication.")
                        call setMatCopy(inv(1:,2:), rdpack, inv(2:,1:), rdpack, lowDia, transHerm)
        call disp%show("inv")
        call disp%show( inv , format = cform )
        call disp%show("mul = matmul(mat, inv)")
                        mul = matmul(mat, inv)
        call disp%show("mul")
        call disp%show( mul , format = cform )
        call disp%skip
        call disp%show("inv = 0")
                        inv = 0
        call disp%show("call setMatInv(inv, chol, auxil = choUpp, subset = uppDia)")
                        call setMatInv(inv, chol, auxil = choUpp, subset = uppDia)
        call disp%show("inv")
        call disp%show( inv , format = cform )
        call disp%show("call setMatCopy(inv(2:,1:), rdpack, inv(1:,2:), rdpack, uppDia, transHerm) ! symmetrize `inv` for multiplication.")
                        call setMatCopy(inv(2:,1:), rdpack, inv(1:,2:), rdpack, uppDia, transHerm)
        call disp%show("inv")
        call disp%show( inv , format = cform )
        call disp%show("mul = matmul(mat, inv)")
                        mul = matmul(mat, inv)
        call disp%show("mul")
        call disp%show( mul , format = cform )
        call disp%skip
    end block

    block
        use pm_kind, only: TKC => RKS
        complex(TKC), allocatable :: low(:,:), upp(:,:), mat(:,:), chol(:,:), inv(:,:), mul(:,:)
        do itry = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            low = getUnifRand((1._TKC, -1._TKC), (+2._TKC, +1._TKC), ndim, ndim)
            call setMatInit(low, uppDia, (0._TKC, 0._TKC), cmplx(getUnifRand(1._TKC, 2._TKC, ndim), 0._TKC, TKC))
            upp = transpose(conjg(low))
            mat = matmul(low, upp)
            call disp%show("mat")
            call disp%show( mat , format = cform )
            call disp%show("chol = getMatChol(mat, subset = uppDia)")
                            chol = getMatChol(mat, subset = uppDia)
            call disp%show("chol")
            call disp%show( chol , format = cform )
            call disp%show("call setResized(inv, shape(mat, IK))")
                            call setResized(inv, shape(mat, IK))
            call disp%show("inv = 0")
                            inv = 0
            call disp%show("call setMatInv(inv, chol, auxil = choUpp, subset = lowDia)")
                            call setMatInv(inv, chol, auxil = choUpp, subset = lowDia)
            call disp%show("inv")
            call disp%show( inv , format = cform )
            call disp%show("call setMatCopy(inv(1:,2:), rdpack, inv(2:,1:), rdpack, lowDia, transHerm) ! symmetrize `inv` for multiplication.")
                            call setMatCopy(inv(1:,2:), rdpack, inv(2:,1:), rdpack, lowDia, transHerm)
            call disp%show("inv")
            call disp%show( inv , format = cform )
            call disp%show("mul = matmul(mat, inv)")
                            mul = matmul(mat, inv)
            call disp%show("mul")
            call disp%show( mul , format = cform )
            call disp%skip
            call disp%show("inv = 0")
                            inv = 0
            call disp%show("call setMatInv(inv, chol, auxil = choUpp, subset = uppDia)")
                            call setMatInv(inv, chol, auxil = choUpp, subset = uppDia)
            call disp%show("inv")
            call disp%show( inv , format = cform )
            call disp%show("call setMatCopy(inv(2:,1:), rdpack, inv(1:,2:), rdpack, uppDia, transHerm) ! symmetrize `inv` for multiplication.")
                            call setMatCopy(inv(2:,1:), rdpack, inv(1:,2:), rdpack, uppDia, transHerm)
            call disp%show("inv")
            call disp%show( inv , format = cform )
            call disp%show("mul = matmul(mat, inv)")
                            mul = matmul(mat, inv)
            call disp%show("mul")
            call disp%show( mul , format = cform )
            call disp%skip
        end do
    end block

end program example