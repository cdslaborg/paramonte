program example

    use pm_kind, only: SK, IK, LK
    use pm_matrixDet, only: getMatDet
    use pm_matrixDet, only: setMatDet
    use pm_matrixInv, only: getMatInv
    use pm_distUnif, only: getUnifRand
    use pm_io, only: display_type

    implicit none

    integer(IK) :: info, ndim, itry, ntry = 10
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the determinant of the square matrix.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: TKC => RKH
        real(TKC), allocatable :: mat(:,:), lup(:,:)
        real(TKC) :: det
        do itry = 1, ntry
            call disp%skip()
            call disp%show("ndim = getUnifRand(1, 9)")
                            ndim = getUnifRand(1, 9)
            call disp%show("mat = getUnifRand(1._TKC, 2._TKC, ndim, ndim)")
                            mat = getUnifRand(1._TKC, 2._TKC, ndim, ndim)
            call disp%show("lup = mat")
                            lup = mat
            call disp%show("lup ! matrix.")
            call disp%show( lup )
            call disp%show("call setMatDet(lup, det, info)")
                            call setMatDet(lup, det, info)
            call disp%show("lup ! lup factorization.")
            call disp%show( lup )
            call disp%show("if (info /= 0) error stop 'determinant computation failed.'")
                            if (info /= 0) error stop 'determinant computation failed.'
            call disp%show("det * getMatDet(getMatInv(mat)) ! must be one.")
            call disp%show( det * getMatDet(getMatInv(mat)) )
            call disp%skip()
        end do
    end block

    block
        use pm_kind, only: TKC => CKH
        complex(TKC), allocatable :: mat(:,:), lup(:,:)
        complex(TKC) :: det
        do itry = 1, ntry
            call disp%skip()
            call disp%show("ndim = getUnifRand(1, 9)")
                            ndim = getUnifRand(1, 9)
            call disp%show("mat = getUnifRand((1._TKC, 1._TKC), (2._TKC, 2._TKC), ndim, ndim)")
                            mat = getUnifRand((1._TKC, 1._TKC), (2._TKC, 2._TKC), ndim, ndim)
            call disp%show("lup = mat")
                            lup = mat
            call disp%show("lup ! matrix.")
            call disp%show( lup )
            call disp%show("call setMatDet(lup, det, info)")
                            call setMatDet(lup, det, info)
            call disp%show("if (info /= 0) error stop 'determinant computation failed.'")
                            if (info /= 0) error stop 'determinant computation failed.'
            call disp%show("lup ! lup factorization.")
            call disp%show( lup )
            call disp%show("det * getMatDet(getMatInv(mat)) ! must be one.")
            call disp%show( det * getMatDet(getMatInv(mat)) )
            call disp%skip()
        end do
    end block

end program example