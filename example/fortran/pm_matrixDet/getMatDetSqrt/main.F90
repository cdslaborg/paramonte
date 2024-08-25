program example

    use pm_kind, only: SK, IK, LK
    use pm_matrixCopy, only: rdpack
    use pm_matrixCopy, only: getMatCopy
    use pm_matrixDet, only: getMatDetSqrt
    use pm_matrixTrace, only: getMatMulTrace
    use pm_matrixChol, only: uppDia, lowDia
    use pm_matrixChol, only: getMatChol
    use pm_matrixInv, only: getMatInv
    use pm_distCov, only: getCovRand
    use pm_arrayResize, only: setResized
    use pm_distUnif, only: getUnifRand
    use pm_io, only: display_type

    implicit none

    integer(IK) :: ndim, itry, ntry = 10
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the sqrt of the determinant of the positive definite matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKG => RKD
        real(TKG), allocatable :: mat(:,:)
        real(TKG) :: detSqrt
        do itry = 1, ntry
            call disp%skip()
            call disp%show("ndim = getUnifRand(1, 5)")
                            ndim = getUnifRand(1, 5)
            call disp%show("ndim")
            call disp%show( ndim )
            call disp%show("mat = getCovRand(mold = 1._TKG, ndim = ndim)")
                            mat = getCovRand(mold = 1._TKG, ndim = ndim)
            call disp%show("mat")
            call disp%show( mat )
            call disp%show("detSqrt = getMatDetSqrt(mat)")
                            detSqrt = getMatDetSqrt(mat)
            call disp%show("detSqrt")
            call disp%show( detSqrt )
            call disp%show("getMatMulTrace(getMatChol(mat, uppDia)) ! for comparison.")
            call disp%show( getMatMulTrace(getMatChol(mat, uppDia)) )
            call disp%show("detSqrt * getMatDetSqrt(getMatInv(mat)) ! must be one.")
            call disp%show( detSqrt * getMatDetSqrt(getMatInv(mat)) )
            call disp%skip()
            call disp%show("mat = getCovRand(mold = 1._TKG, ndim = ndim)")
                            mat = getCovRand(mold = 1._TKG, ndim = ndim)
            call disp%show("mat = getMatCopy(rdpack, mat, rdpack, lowDia, init = 0._TKG) ! reset the upper.")
                            mat = getMatCopy(rdpack, mat, rdpack, lowDia, init = 0._TKG)
            call disp%show("mat")
            call disp%show( mat )
            call disp%show("detSqrt = getMatDetSqrt(mat, lowDia)")
                            detSqrt = getMatDetSqrt(mat, lowDia)
            call disp%show("detSqrt")
            call disp%show( detSqrt )
            call disp%show("getMatMulTrace(getMatChol(mat, lowDia)) ! for comparison.")
            call disp%show( getMatMulTrace(getMatChol(mat, lowDia)) )
            call disp%show("detSqrt / getMatMulTrace(getMatChol(mat, lowDia)) ! must be one.")
            call disp%show( detSqrt / getMatMulTrace(getMatChol(mat, lowDia)) )
            call disp%skip()
            call disp%show("mat = getCovRand(mold = 1._TKG, ndim = ndim)")
                            mat = getCovRand(mold = 1._TKG, ndim = ndim)
            call disp%show("mat = getMatCopy(rdpack, mat, rdpack, uppDia, init = 0._TKG) ! reset the lower.")
                            mat = getMatCopy(rdpack, mat, rdpack, uppDia, init = 0._TKG)
            call disp%show("mat")
            call disp%show( mat )
            call disp%show("detSqrt = getMatDetSqrt(mat, uppDia)")
                            detSqrt = getMatDetSqrt(mat, uppDia)
            call disp%show("detSqrt")
            call disp%show( detSqrt )
            call disp%show("getMatMulTrace(getMatChol(mat, uppDia)) ! for comparison.")
            call disp%show( getMatMulTrace(getMatChol(mat, uppDia)) )
            call disp%show("detSqrt / getMatMulTrace(getMatChol(mat, uppDia)) ! must be one.")
            call disp%show( detSqrt / getMatMulTrace(getMatChol(mat, uppDia)) )
            call disp%skip()
        end do
    end block

    block
        use pm_kind, only: TKG => CKS
        complex(TKG), allocatable :: tmp(:,:)
        complex(TKG), parameter :: mat(*,*) = reshape(  [ (9.0,  0.0), (3.0, 3.0), (3.0, -3.0) &
                                                        , (3.0, -3.0),(18.0, 0.0), (8.0, -6.0) &
                                                        , (3.0,  3.0), (8.0, 6.0),(43.0,  0.0) &
                                                        ], shape = [3, 3], order = [2, 1])
        real(TKG) :: detSqrt
        call disp%skip()
        call disp%show("mat")
        call disp%show( mat )
        call disp%show("detSqrt = getMatDetSqrt(mat)")
                        detSqrt = getMatDetSqrt(mat)
        call disp%show("detSqrt")
        call disp%show( detSqrt )
        call disp%show("getMatMulTrace(getMatChol(mat, uppDia)) ! for comparison.")
        call disp%show( getMatMulTrace(getMatChol(mat, uppDia)) )
        call disp%show("detSqrt * getMatDetSqrt(getMatInv(mat)) ! must be one.")
        call disp%show( detSqrt * getMatDetSqrt(getMatInv(mat)) )
        call disp%skip()
        call disp%show("tmp = getMatCopy(rdpack, mat, rdpack, lowDia, init = (0._TKG, 0._TKG)) ! reset the upper.")
                        tmp = getMatCopy(rdpack, mat, rdpack, lowDia, init = (0._TKG, 0._TKG))
        call disp%show("tmp")
        call disp%show( tmp )
        call disp%show("detSqrt = getMatDetSqrt(tmp, subset = lowDia)")
                        detSqrt = getMatDetSqrt(tmp, subset = lowDia)
        call disp%show("detSqrt")
        call disp%show( detSqrt )
        call disp%show("getMatMulTrace(getMatChol(tmp, lowDia)) ! for comparison.")
        call disp%show( getMatMulTrace(getMatChol(tmp, lowDia)) )
        call disp%show("detSqrt / getMatMulTrace(getMatChol(tmp, lowDia)) ! must be one")
        call disp%show( detSqrt / getMatMulTrace(getMatChol(tmp, lowDia)) )
        call disp%skip()
        call disp%show("tmp = getMatCopy(rdpack, mat, rdpack, uppDia, init = (0._TKG, 0._TKG)) ! reset the upper.")
                        tmp = getMatCopy(rdpack, mat, rdpack, uppDia, init = (0._TKG, 0._TKG))
        call disp%show("tmp")
        call disp%show( tmp )
        call disp%show("detSqrt = getMatDetSqrt(tmp, subset = uppDia)")
                        detSqrt = getMatDetSqrt(tmp, subset = uppDia)
        call disp%show("detSqrt")
        call disp%show( detSqrt )
        call disp%show("getMatMulTrace(getMatChol(tmp, uppDia)) ! for comparison.")
        call disp%show( getMatMulTrace(getMatChol(tmp, uppDia)) )
        call disp%show("detSqrt / getMatMulTrace(getMatChol(tmp, uppDia)) ! must be one.")
        call disp%show( detSqrt / getMatMulTrace(getMatChol(tmp, uppDia)) )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKG => CKS
        complex(TKG), allocatable :: tmp(:,:)
        complex(TKG), parameter :: mat(*,*) = reshape(  [  (25.0, 0.0), (-5.0, -5.0), (10.0, 5.0) &
                                                        ,  (-5.0, 5.0),  (51.0, 0.0), (4.0, -6.0) &
                                                        , (10.0, -5.0),   (4.0, 6.0), (71.0, 0.0) &
                                                        ], shape = [3, 3], order = [2, 1])
        real(TKG) :: detSqrt
        call disp%skip()
        call disp%show("mat")
        call disp%show( mat )
        call disp%show("detSqrt = getMatDetSqrt(mat)")
                        detSqrt = getMatDetSqrt(mat)
        call disp%show("detSqrt")
        call disp%show( detSqrt )
        call disp%show("getMatMulTrace(getMatChol(mat, uppDia)) ! for comparison.")
        call disp%show( getMatMulTrace(getMatChol(mat, uppDia)) )
        call disp%show("detSqrt * getMatDetSqrt(getMatInv(mat)) ! must be one.")
        call disp%show( detSqrt * getMatDetSqrt(getMatInv(mat)) )
        call disp%skip()
        call disp%show("tmp = getMatCopy(rdpack, mat, rdpack, lowDia, init = (0._TKG, 0._TKG)) ! reset the upper.")
                        tmp = getMatCopy(rdpack, mat, rdpack, lowDia, init = (0._TKG, 0._TKG))
        call disp%show("tmp")
        call disp%show( tmp )
        call disp%show("detSqrt = getMatDetSqrt(tmp, subset = lowDia)")
                        detSqrt = getMatDetSqrt(tmp, subset = lowDia)
        call disp%show("detSqrt")
        call disp%show( detSqrt )
        call disp%show("getMatMulTrace(getMatChol(tmp, lowDia)) ! for comparison.")
        call disp%show( getMatMulTrace(getMatChol(tmp, lowDia)) )
        call disp%show("detSqrt / getMatMulTrace(getMatChol(tmp, lowDia)) ! must be one.")
        call disp%show( detSqrt / getMatMulTrace(getMatChol(tmp, lowDia)) )
        call disp%skip()
        call disp%show("tmp = getMatCopy(rdpack, mat, rdpack, uppDia, init = (0._TKG, 0._TKG)) ! reset the upper.")
                        tmp = getMatCopy(rdpack, mat, rdpack, uppDia, init = (0._TKG, 0._TKG))
        call disp%show("tmp")
        call disp%show( tmp )
        call disp%show("detSqrt = getMatDetSqrt(tmp, subset = uppDia)")
                        detSqrt = getMatDetSqrt(tmp, subset = uppDia)
        call disp%show("detSqrt")
        call disp%show( detSqrt )
        call disp%show("getMatMulTrace(getMatChol(tmp, uppDia)) ! for comparison.")
        call disp%show( getMatMulTrace(getMatChol(tmp, uppDia)) )
        call disp%show("detSqrt / getMatMulTrace(getMatChol(tmp, uppDia)) ! must be one.")
        call disp%show( detSqrt / getMatMulTrace(getMatChol(tmp, uppDia)) )
        call disp%skip()
    end block

end program example