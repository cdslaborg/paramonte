program example

    use pm_io, only: getFormat
    use pm_kind, only: SK, IK, LK, TKG => RKS
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRand
    use pm_matrixInv, only: getMatInv
    use pm_matrixInv, only: choUpp, choLow
    use pm_matrixInv, only: upperDiag, lowerDiag
    use pm_matrixInv, only: upperUnit, lowerUnit
    use pm_matrixChol, only: getMatChol, uppDia, lowDia
    use pm_matrixInit, only: getMatInit, uppLowDia
    use pm_distCov, only: getCovRand
    use pm_arrayResize, only: setResized
    use pm_matrixInit, only: setMatInit
    use pm_arrayFill, only: getFilled
    use pm_matrixLUP, only: setMatLUP
    use pm_err, only: setAsserted
    use pm_val2str, only: getStr

    implicit none

    integer(IK):: i, info, ndim, ntry = 10
    character(:, SK), allocatable :: cform, rform

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    cform = getFormat(mold = [(0._TKG, 0._TKG)], ed = SK_"f", ndigit = 2_IK, signed = .true._LK)
    rform = getFormat(mold = [0._TKG], ed = SK_"f", ndigit = 2_IK, signed = .true._LK)

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute inverse of a general    real matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => RKS
        real(TKG), allocatable :: mat(:,:), inv(:,:), mul(:,:)
        real(TKG) :: invDetSqrt
        do i = 1, 1
            call disp%skip
            call disp%show("mat = reshape([1, 0, 2, -1, 5, 0, 0, 3, -9], shape = [3,3], order = [2, 1])")
                            mat = reshape([1, 0, 2, -1, 5, 0, 0, 3, -9], shape = [3,3], order = [2, 1])
            call disp%show("mat")
            call disp%show( mat , format = rform )
            call disp%show("inv = getMatInv(mat)")
                            inv = getMatInv(mat)
            call disp%show("inv")
            call disp%show( inv , format = rform )
            call disp%show("inv = getMatInv(mat, info) ! gracefully catch inversion errors.")
                            inv = getMatInv(mat, info)
            call disp%show("if (info /= 0) error stop 'inversion failed.'")
                            if (info /= 0) error stop 'inversion failed.'
            call disp%show("inv")
            call disp%show( inv , format = rform )
            call disp%show("inv = getMatInv(mat, invDetSqrt) ! compute the sqrt of the determinant of the inverse of the general matrix along with the inverse.")
                            inv = getMatInv(mat, invDetSqrt)
            call disp%show("invDetSqrt")
            call disp%show( invDetSqrt )
            call disp%show("inv")
            call disp%show( inv , format = rform )
            call disp%show("mul = matmul(mat, inv)")
                            mul = matmul(mat, inv)
            call disp%show("mul")
            call disp%show( mul, format = rform )
            call disp%skip
        end do
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute inverse of a general complex matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => RKS
        complex(TKG) :: invDetSqrt
        complex(TKG), allocatable :: inv(:,:), mul(:,:)
        complex(TKG), parameter :: mat(*,*) = reshape( &
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
        call disp%skip
        call disp%show("mat")
        call disp%show( mat , format = cform )
        call disp%show("inv = getMatInv(mat)")
                        inv = getMatInv(mat)
        call disp%show("inv")
        call disp%show( inv , format = cform )
        call disp%show("inv = getMatInv(mat, info) ! gracefully catch inversion errors.")
                        inv = getMatInv(mat, info)
        call disp%show("if (info /= 0) error stop 'inversion failed.'")
                        if (info /= 0) error stop 'inversion failed.'
        call disp%show("inv")
        call disp%show( inv , format = cform )
        call disp%show("inv = getMatInv(mat, invDetSqrt) ! compute the sqrt of the determinant of the inverse of the general matrix along with the inverse.")
                        inv = getMatInv(mat, invDetSqrt)
        call disp%show("invDetSqrt")
        call disp%show( invDetSqrt )
        call disp%show("inv")
        call disp%show( inv , format = cform )
        call disp%show("mul = matmul(mat, inv)")
                        mul = matmul(mat, inv)
        call disp%show("mul")
        call disp%show( mul, format = cform )
        call disp%skip
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute inverse of an upperDiag triangular    real matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => RKS
        real(TKG), allocatable :: mat(:,:), inv(:,:), mul(:,:)
        do i = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getUnifRand(1_IK, 9_IK, ndim, ndim)")
                            mat = getUnifRand(1_IK, 9_IK, ndim, ndim)
            call disp%show("call setMatInit(mat(2 : ndim, 1 : ndim - 1), lowDia, 0._TKG, 0._TKG)")
                            call setMatInit(mat(2 : ndim, 1 : ndim - 1), lowDia, 0._TKG, 0._TKG)
            call disp%show("mat")
            call disp%show( mat , format = rform )
            call disp%show("inv = getMatInv(mat, upperDiag)")
                            inv = getMatInv(mat, upperDiag)
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
        use pm_kind, only: TKG => RKS
        complex(TKG), allocatable :: mat(:,:), inv(:,:), mul(:,:)
        do i = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getUnifRand((1., 1.), (2., 2.), ndim, ndim)")
                            mat = getUnifRand((1., 1.), (2., 2.), ndim, ndim)
            call disp%show("call setMatInit(mat(2 : ndim, 1 : ndim - 1), lowDia, (0._TKG, 0._TKG), (0._TKG, 0._TKG))")
                            call setMatInit(mat(2 : ndim, 1 : ndim - 1), lowDia, (0._TKG, 0._TKG), (0._TKG, 0._TKG))
            call disp%show("mat")
            call disp%show( mat , format = cform )
            call disp%show("inv = getMatInv(mat, upperDiag)")
                            inv = getMatInv(mat, upperDiag)
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
        use pm_kind, only: TKG => RKS
        real(TKG), allocatable :: mat(:,:), inv(:,:), mul(:,:)
        do i = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getUnifRand(1_IK, 9_IK, ndim, ndim)")
                            mat = getUnifRand(1_IK, 9_IK, ndim, ndim)
            call disp%show("call setMatInit(mat(1 : ndim - 1, 2 : ndim), uppDia, 0._TKG, 0._TKG)")
                            call setMatInit(mat(1 : ndim - 1, 2 : ndim), uppDia, 0._TKG, 0._TKG)
            call disp%show("mat")
            call disp%show( mat , format = rform )
            call disp%show("inv = getMatInv(mat, lowerDiag)")
                            inv = getMatInv(mat, lowerDiag)
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
        use pm_kind, only: TKG => RKS
        complex(TKG), allocatable :: mat(:,:), inv(:,:), mul(:,:)
        do i = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getUnifRand((1., 1.), (2., 2.), ndim, ndim)")
                            mat = getUnifRand((1., 1.), (2., 2.), ndim, ndim)
            call disp%show("call setMatInit(mat(1 : ndim - 1, 2 : ndim), uppDia, (0._TKG, 0._TKG), (0._TKG, 0._TKG))")
                            call setMatInit(mat(1 : ndim - 1, 2 : ndim), uppDia, (0._TKG, 0._TKG), (0._TKG, 0._TKG))
            call disp%show("mat")
            call disp%show( mat , format = cform )
            call disp%show("inv = getMatInv(mat, lowerDiag)")
                            inv = getMatInv(mat, lowerDiag)
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
        use pm_kind, only: TKG => RKS
        real(TKG), allocatable :: mat(:,:), inv(:,:), mul(:,:)
        do i = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getUnifRand(1_IK, 9_IK, ndim, ndim)")
                            mat = getUnifRand(1_IK, 9_IK, ndim, ndim)
            call disp%show("call setMatInit(mat, lowDia, 0._TKG, 1._TKG)")
                            call setMatInit(mat, lowDia, 0._TKG, 1._TKG)
            call disp%show("mat")
            call disp%show( mat , format = rform )
            call disp%show("inv = getMatInv(mat, upperUnit)")
                            inv = getMatInv(mat, upperUnit)
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
        use pm_kind, only: TKG => RKS
        complex(TKG), allocatable :: mat(:,:), inv(:,:), mul(:,:)
        do i = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getUnifRand((-1., -1.), (+1., +1.), ndim, ndim)")
                            mat = getUnifRand((-1., -1.), (+1., +1.), ndim, ndim)
            call disp%show("call setMatInit(mat, lowDia, (0._TKG, 0._TKG), (1._TKG, 0._TKG))")
                            call setMatInit(mat, lowDia, (0._TKG, 0._TKG), (1._TKG, 0._TKG))
            call disp%show("mat")
            call disp%show( mat , format = cform )
            call disp%show("inv = getMatInv(mat, upperUnit)")
                            inv = getMatInv(mat, upperUnit)
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
        use pm_kind, only: TKG => RKS
        real(TKG), allocatable :: mat(:,:), inv(:,:), mul(:,:)
        do i = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getUnifRand(1_IK, 9_IK, ndim, ndim)")
                            mat = getUnifRand(1_IK, 9_IK, ndim, ndim)
            call disp%show("call setMatInit(mat, uppDia, 0._TKG, 1._TKG)")
                            call setMatInit(mat, uppDia, 0._TKG, 1._TKG)
            call disp%show("mat")
            call disp%show( mat , format = rform )
            call disp%show("inv = getMatInv(mat, lowerUnit)")
                            inv = getMatInv(mat, lowerUnit)
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
        use pm_kind, only: TKG => RKS
        complex(TKG), allocatable :: mat(:,:), inv(:,:), mul(:,:)
        do i = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getUnifRand((-1., -1.), (+1., +1.), ndim, ndim)")
                            mat = getUnifRand((-1., -1.), (+1., +1.), ndim, ndim)
            call disp%show("call setMatInit(mat, uppDia, (0._TKG, 0._TKG), (1._TKG, 0._TKG))")
                            call setMatInit(mat, uppDia, (0._TKG, 0._TKG), (1._TKG, 0._TKG))
            call disp%show("mat")
            call disp%show( mat , format = cform )
            call disp%show("inv = getMatInv(mat, lowerUnit)")
                            inv = getMatInv(mat, lowerUnit)
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
        use pm_kind, only: TKG => RKS
        integer(IK), allocatable :: rperm(:)
        real(TKG), allocatable :: mat(:,:), lup(:,:), inv(:,:), mul(:,:)
        do i = 1, ntry
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
            call disp%show("inv = getMatInv(lup, rperm) ! compute the inverse of the matrix by passing its LUP factorization.")
                            inv = getMatInv(lup, rperm)
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
        use pm_kind, only: TKG => RKS
        integer(IK), allocatable :: rperm(:)
        complex(TKG), parameter :: mat(*,*) = reshape( &
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
        complex(TKG), dimension(size(mat,1), size(mat,2)) :: inv, lup, mul
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
        call disp%show("inv = getMatInv(lup, rperm) ! compute the inverse of the matrix by passing its LUP factorization.")
                        inv = getMatInv(lup, rperm)
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
        use pm_kind, only: TKG => RKS
        real(TKG), allocatable :: mat(:,:), chol(:,:), inv(:,:), mul(:,:)
        do i = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getCovRand(mold = 1._TKG, ndim = ndim)")
                            mat = getCovRand(mold = 1._TKG, ndim = ndim)
            call disp%show("mat")
            call disp%show( mat , format = rform )
            call disp%show("chol = getMatChol(mat, subset = lowDia)")
                            chol = getMatChol(mat, subset = lowDia)
            call disp%show("chol")
            call disp%show( chol , format = rform )
            call disp%show("inv = getMatInv(chol, auxil = choLow)")
                            inv = getMatInv(chol, auxil = choLow)
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
        use pm_kind, only: TKG => RKS
        real(TKG), allocatable :: mat(:,:), chol(:,:), inv(:,:), mul(:,:)
        do i = 1, ntry
            call disp%skip
            call disp%show("ndim = getUnifRand(1_IK, 6_IK)")
                            ndim = getUnifRand(1_IK, 6_IK)
            call disp%show("ndim ! matrix rank")
            call disp%show( ndim )
            call disp%show("mat = getCovRand(mold = 1._TKG, ndim = ndim)")
                            mat = getCovRand(mold = 1._TKG, ndim = ndim)
            call disp%show("mat")
            call disp%show( mat , format = rform )
            call disp%show("chol = getMatChol(mat, subset = uppDia)")
                            chol = getMatChol(mat, subset = uppDia)
            call disp%show("chol")
            call disp%show( chol , format = rform )
            call disp%show("inv = getMatInv(chol, auxil = choUpp)")
                            inv = getMatInv(chol, auxil = choUpp)
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
        use pm_kind, only: TKG => RKS
        integer(IK), parameter :: ndim = 3_IK
        complex(TKG), allocatable :: chol(:,:), inv(:,:), mul(:,:)
        complex(TKG), parameter :: mat(*,*) = reshape(  [ (9.0,  0.0), (3.0, 3.0), (3.0, -3.0) &
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
        call disp%show("inv = getMatInv(chol, auxil = choLow)")
                        inv = getMatInv(chol, auxil = choLow)
        call disp%show("inv")
        call disp%show( inv , format = cform )
        call disp%show("mul = matmul(mat, inv)")
                        mul = matmul(mat, inv)
        call disp%show("mul")
        call disp%show( mul , format = cform )
        call disp%skip
    end block

    block
        use pm_kind, only: TKG => RKS
        integer(IK), parameter :: ndim = 3_IK
        complex(TKG), allocatable :: chol(:,:), inv(:,:), mul(:,:)
        complex(TKG), parameter :: mat(*,*) = reshape(  [ (9.0,  0.0), (3.0, 3.0), (3.0, -3.0) &
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
        call disp%show("inv = getMatInv(chol, auxil = choUpp)")
                        inv = getMatInv(chol, auxil = choUpp)
        call disp%show("inv")
        call disp%show( inv , format = cform )
        call disp%show("mul = matmul(mat, inv)")
                        mul = matmul(mat, inv)
        call disp%show("mul")
        call disp%show( mul , format = cform )
        call disp%skip
    end block

end program example