program example

    use pm_kind, only: SK, IK, LK
    use pm_matrixDet, only: getMatDet
    use pm_matrixInv, only: getMatInv
    use pm_distUnif, only: getUnifRand
    use pm_io, only: display_type

    implicit none

    integer(IK) :: ndim, itry, ntry = 10
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the determinant of the square matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKH
        real(TKC), allocatable :: mat(:,:)
        real(TKC) :: det
        do itry = 1, ntry
            call disp%skip()
            call disp%show("ndim = getUnifRand(1, 9)")
                            ndim = getUnifRand(1, 9)
            call disp%show("ndim")
            call disp%show( ndim )
            call disp%show("mat = getUnifRand(1._TKC, 2._TKC, ndim, ndim)")
                            mat = getUnifRand(1._TKC, 2._TKC, ndim, ndim)
            call disp%show("mat")
            call disp%show( mat )
            call disp%show("det = getMatDet(mat)")
                            det = getMatDet(mat)
            call disp%show("det")
            call disp%show( det )
            call disp%show("det * getMatDet(getMatInv(mat)) ! must be one.")
            call disp%show( det * getMatDet(getMatInv(mat)) )
            call disp%skip()
        end do
    end block

    block
        use pm_kind, only: TKC => CKH
        complex(TKC), allocatable :: mat(:,:)
        real(TKC) :: det
        do itry = 1, ntry
            call disp%skip()
            call disp%show("ndim = getUnifRand(1, 9)")
                            ndim = getUnifRand(1, 9)
            call disp%show("mat = getUnifRand((1._TKC, 1._TKC), (2._TKC, 2._TKC), ndim, ndim)")
                            mat = getUnifRand((1._TKC, 1._TKC), (2._TKC, 2._TKC), ndim, ndim)
            call disp%show("mat")
            call disp%show( mat )
            call disp%show("det = getMatDet(mat)")
                            det = getMatDet(mat)
            call disp%show("det")
            call disp%show( det )
            call disp%show("det * getMatDet(getMatInv(mat)) ! must be one.")
            call disp%show( det * getMatDet(getMatInv(mat)) )
            call disp%skip()
        end do
    end block

end program example