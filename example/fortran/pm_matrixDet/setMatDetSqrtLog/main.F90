program example

    use pm_kind, only: SK, IK, LK
    use pm_matrixCopy, only: rdpack
    use pm_matrixCopy, only: transHerm
    use pm_matrixCopy, only: getMatCopy
    use pm_matrixDet, only: getMatDetSqrtLog
    use pm_matrixDet, only: setMatDetSqrtLog
    use pm_matrixTrace, only: getMatMulTraceLog
    use pm_matrixChol, only: uppDia, lowDia
    use pm_matrixChol, only: getMatChol
    use pm_distCov, only: getCovRand
    use pm_arrayResize, only: setResized
    use pm_distUnif, only: getUnifRand
    use pm_matrixInv, only: getMatInv
    use pm_io, only: display_type
    use pm_io, only: getFormat

    implicit none

    integer(IK) :: info, ndim, itry, ntry = 10
    character(:, SK), allocatable :: format
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the sqrt of the determinant of the positive definite matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKS
        real(TKC), allocatable :: mat(:,:)
        real(TKC) :: detSqrtLog
        format = getFormat(mold = [0._TKC], ed = SK_"e")
        do itry = 1, ntry
            call disp%skip()
            call disp%show("ndim = getUnifRand(1, 5)")
                            ndim = getUnifRand(1, 5)
            call disp%show("ndim")
            call disp%show( ndim )
            call disp%show("call setResized(mat, [ndim, ndim + 1])")
                            call setResized(mat, [ndim, ndim + 1])
            call disp%show("mat(:,1:ndim) = getCovRand(mold = 1._TKC, ndim = ndim)")
                            mat(:,1:ndim) = getCovRand(mold = 1._TKC, ndim = ndim)
            call disp%show("mat(:,1:ndim)")
            call disp%show( mat(:,1:ndim) , format = format )
            call disp%show("mat(:,1:ndim) = getMatCopy(rdpack, mat(:,1:ndim), rdpack, lowDia, init = 0._TKC) ! reset the upper.")
                            mat(:,1:ndim) = getMatCopy(rdpack, mat(:,1:ndim), rdpack, lowDia, init = 0._TKC)
            call disp%show("mat(:,1:ndim)")
            call disp%show( mat(:,1:ndim) , format = format )
            call disp%show("call setMatDetSqrtLog(mat(:,1:ndim), lowDia, detSqrtLog, info, mat(:,2:ndim+1), transHerm)")
                            call setMatDetSqrtLog(mat(:,1:ndim), lowDia, detSqrtLog, info, mat(:,2:ndim+1), transHerm)
            call disp%show("if (info /= 0) error stop 'Cholesky factorization failed.'")
                            if (info /= 0) error stop 'Cholesky factorization failed.'
            call disp%show("detSqrtLog")
            call disp%show( detSqrtLog )
            call disp%show("mat")
            call disp%show( mat , format = format )
            call disp%show("getMatChol(mat(:,1:ndim), lowDia)")
            call disp%show( getMatChol(mat(:,1:ndim), lowDia) , format = format )
            call disp%show("getMatMulTraceLog(getMatChol(mat(:,1:ndim), lowDia)) ! for comparison.")
            call disp%show( getMatMulTraceLog(getMatChol(mat(:,1:ndim), lowDia)) )
            call disp%show("detSqrtLog - getMatMulTraceLog(getMatChol(mat(:,1:ndim), lowDia)) ! must be one.")
            call disp%show( detSqrtLog - getMatMulTraceLog(getMatChol(mat(:,1:ndim), lowDia)) )
            call disp%skip()
            call disp%show("mat(:,2:ndim+1) = getCovRand(mold = 1._TKC, ndim = ndim)")
                            mat(:,2:ndim+1) = getCovRand(mold = 1._TKC, ndim = ndim)
            call disp%show("mat(:,2:ndim+1) = getMatCopy(rdpack, mat(:,2:ndim+1), rdpack, uppDia, init = 0._TKC) ! reset the lower.")
                            mat(:,2:ndim+1) = getMatCopy(rdpack, mat(:,2:ndim+1), rdpack, uppDia, init = 0._TKC)
            call disp%show("mat(:,2:ndim+1)")
            call disp%show( mat(:,2:ndim+1) , format = format )
            call disp%show("call setMatDetSqrtLog(mat(:,2:ndim+1), uppDia, detSqrtLog, info, mat(:,1:ndim), transHerm)")
                            call setMatDetSqrtLog(mat(:,2:ndim+1), uppDia, detSqrtLog, info, mat(:,1:ndim), transHerm)
            call disp%show("if (info /= 0) error stop 'Cholesky factorization failed.'")
                            if (info /= 0) error stop 'Cholesky factorization failed.'
            call disp%show("detSqrtLog")
            call disp%show( detSqrtLog )
            call disp%show("mat")
            call disp%show( mat , format = format )
            call disp%show("getMatChol(mat(:,2:ndim+1), uppDia)")
            call disp%show( getMatChol(mat(:,2:ndim+1), uppDia) , format = format )
            call disp%show("getMatMulTraceLog(getMatChol(mat(:,2:ndim+1), uppDia)) ! for comparison.")
            call disp%show( getMatMulTraceLog(getMatChol(mat(:,2:ndim+1), uppDia)) )
            call disp%show("detSqrtLog - getMatMulTraceLog(getMatChol(mat(:,2:ndim+1), uppDia)) ! must be one.")
            call disp%show( detSqrtLog - getMatMulTraceLog(getMatChol(mat(:,2:ndim+1), uppDia)) )
            call disp%skip()
        end do
    end block

    block
        use pm_kind, only: TKC => CKS
        integer(IK), parameter :: ndim = 3
        complex(TKC), allocatable :: tmp(:,:)
        complex(TKC), parameter :: mat(*,*) = reshape(  [ (9.0,  0.0), (3.0, 3.0), (3.0, -3.0) &
                                                        , (3.0, -3.0),(18.0, 0.0), (8.0, -6.0) &
                                                        , (3.0,  3.0), (8.0, 6.0),(43.0,  0.0) &
                                                        ], shape = [ndim, ndim], order = [2, 1])
        real(TKC) :: detSqrtLog
        format = getFormat(mold = [(0._TKC, 0._TKC)], ed = SK_"e")
        call disp%skip()
        call disp%show("mat")
        call disp%show( mat , format = format )
        call disp%show("call setResized(tmp, [ndim, ndim + 1])")
                        call setResized(tmp, [ndim, ndim + 1])
        call disp%show("tmp(:,1:ndim) = getMatCopy(rdpack, mat, rdpack, lowDia, init = (0._TKC, 0._TKC)) ! reset the upper.")
                        tmp(:,1:ndim) = getMatCopy(rdpack, mat, rdpack, lowDia, init = (0._TKC, 0._TKC))
        call disp%show("tmp(:,1:ndim)")
        call disp%show( tmp(:,1:ndim) , format = format )
        call disp%show("call setMatDetSqrtLog(tmp(:,1:ndim), lowDia, detSqrtLog, info, tmp(:,2:ndim+1), transHerm)")
                        call setMatDetSqrtLog(tmp(:,1:ndim), lowDia, detSqrtLog, info, tmp(:,2:ndim+1), transHerm)
        call disp%show("if (info /= 0) error stop 'Cholesky factorization failed.'")
                        if (info /= 0) error stop 'Cholesky factorization failed.'
        call disp%show("detSqrtLog")
        call disp%show( detSqrtLog )
        call disp%show("tmp")
        call disp%show( tmp , format = format )
        call disp%show("getMatMulTraceLog(getMatChol(tmp(:,1:ndim), lowDia)) ! for comparison.")
        call disp%show( getMatMulTraceLog(getMatChol(tmp(:,1:ndim), lowDia)) )
        call disp%show("detSqrtLog - getMatMulTraceLog(getMatChol(tmp(:,1:ndim), lowDia)) ! must be one")
        call disp%show( detSqrtLog - getMatMulTraceLog(getMatChol(tmp(:,1:ndim), lowDia)) )
        call disp%skip()
        call disp%show("tmp(:,2:ndim+1) = getMatCopy(rdpack, mat(:,1:ndim), rdpack, uppDia, init = (0._TKC, 0._TKC)) ! reset the upper.")
                        tmp(:,2:ndim+1) = getMatCopy(rdpack, mat(:,1:ndim), rdpack, uppDia, init = (0._TKC, 0._TKC))
        call disp%show("tmp(:,2:ndim+1)")
        call disp%show( tmp(:,2:ndim+1) , format = format )
        call disp%show("call setMatDetSqrtLog(tmp(:,2:ndim+1), uppDia, detSqrtLog, info, tmp(:,1:ndim), transHerm)")
                        call setMatDetSqrtLog(tmp(:,2:ndim+1), uppDia, detSqrtLog, info, tmp(:,1:ndim), transHerm)
        call disp%show("if (info /= 0) error stop 'Cholesky factorization failed.'")
                        if (info /= 0) error stop 'Cholesky factorization failed.'
        call disp%show("detSqrtLog")
        call disp%show( detSqrtLog )
        call disp%show("tmp")
        call disp%show( tmp , format = format )
        call disp%show("getMatMulTraceLog(getMatChol(tmp(:,2:ndim+1), uppDia)) ! for comparison.")
        call disp%show( getMatMulTraceLog(getMatChol(tmp(:,2:ndim+1), uppDia)) )
        call disp%show("detSqrtLog - getMatMulTraceLog(getMatChol(tmp(:,2:ndim+1), uppDia)) ! must be one.")
        call disp%show( detSqrtLog - getMatMulTraceLog(getMatChol(tmp(:,2:ndim+1), uppDia)) )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => CKS
        integer(IK), parameter :: ndim = 3
        complex(TKC), allocatable :: tmp(:,:)
        complex(TKC), parameter :: mat(*,*) = reshape(  [  (25.0, 0.0), (-5.0, -5.0), (10.0, 5.0) &
                                                        ,  (-5.0, 5.0),  (51.0, 0.0), (4.0, -6.0) &
                                                        , (10.0, -5.0),   (4.0, 6.0), (71.0, 0.0) &
                                                        ], shape = [ndim, ndim], order = [2, 1])
        real(TKC) :: detSqrtLog
        format = getFormat(mold = [(0._TKC, 0._TKC)], ed = SK_"e")
        call disp%skip()
        call disp%show("mat")
        call disp%show( mat , format = format )
        call disp%show("call setResized(tmp, [ndim, ndim + 1])")
                        call setResized(tmp, [ndim, ndim + 1])
        call disp%show("tmp(:,1:ndim) = getMatCopy(rdpack, mat, rdpack, lowDia, init = (0._TKC, 0._TKC)) ! reset the upper.")
                        tmp(:,1:ndim) = getMatCopy(rdpack, mat, rdpack, lowDia, init = (0._TKC, 0._TKC))
        call disp%show("tmp(:,1:ndim)")
        call disp%show( tmp(:,1:ndim) , format = format )
        call disp%show("call setMatDetSqrtLog(tmp(:,1:ndim), lowDia, detSqrtLog, info, tmp(:,2:ndim+1), transHerm)")
                        call setMatDetSqrtLog(tmp(:,1:ndim), lowDia, detSqrtLog, info, tmp(:,2:ndim+1), transHerm)
        call disp%show("if (info /= 0) error stop 'Cholesky factorization failed.'")
                        if (info /= 0) error stop 'Cholesky factorization failed.'
        call disp%show("detSqrtLog")
        call disp%show( detSqrtLog )
        call disp%show("tmp")
        call disp%show( tmp , format = format )
        call disp%show("getMatMulTraceLog(getMatChol(tmp(:,1:ndim), lowDia)) ! for comparison.")
        call disp%show( getMatMulTraceLog(getMatChol(tmp(:,1:ndim), lowDia)) )
        call disp%show("detSqrtLog - getMatMulTraceLog(getMatChol(tmp(:,1:ndim), lowDia)) ! must be one")
        call disp%show( detSqrtLog - getMatMulTraceLog(getMatChol(tmp(:,1:ndim), lowDia)) )
        call disp%skip()
        call disp%show("tmp(:,2:ndim+1) = getMatCopy(rdpack, mat(:,1:ndim), rdpack, uppDia, init = (0._TKC, 0._TKC)) ! reset the upper.")
                        tmp(:,2:ndim+1) = getMatCopy(rdpack, mat(:,1:ndim), rdpack, uppDia, init = (0._TKC, 0._TKC))
        call disp%show("tmp(:,2:ndim+1)")
        call disp%show( tmp(:,2:ndim+1) , format = format )
        call disp%show("call setMatDetSqrtLog(tmp(:,2:ndim+1), uppDia, detSqrtLog, info, tmp(:,1:ndim), transHerm)")
                        call setMatDetSqrtLog(tmp(:,2:ndim+1), uppDia, detSqrtLog, info, tmp(:,1:ndim), transHerm)
        call disp%show("if (info /= 0) error stop 'Cholesky factorization failed.'")
                        if (info /= 0) error stop 'Cholesky factorization failed.'
        call disp%show("detSqrtLog")
        call disp%show( detSqrtLog )
        call disp%show("tmp")
        call disp%show( tmp , format = format )
        call disp%show("getMatMulTraceLog(getMatChol(tmp(:,2:ndim+1), uppDia)) ! for comparison.")
        call disp%show( getMatMulTraceLog(getMatChol(tmp(:,2:ndim+1), uppDia)) )
        call disp%show("detSqrtLog - getMatMulTraceLog(getMatChol(tmp(:,2:ndim+1), uppDia)) ! must be one.")
        call disp%show( detSqrtLog - getMatMulTraceLog(getMatChol(tmp(:,2:ndim+1), uppDia)) )
        call disp%skip()
    end block

end program example