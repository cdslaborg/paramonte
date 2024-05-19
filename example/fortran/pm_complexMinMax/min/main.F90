program example

    use pm_kind, only: SK, IK
    use pm_complexMinMax, only: min
    use pm_distUnif, only: getUnifRand
    use pm_io, only: display_type

    implicit none

    integer(IK) :: lb, ub, lenarr
    integer(IK) :: itry, ntry = 5
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: CKG => CKS
        complex(CKG) :: a1, a2
        complex(CKG), allocatable :: array1(:)
        call disp%skip()
        do itry = 1, ntry
            call disp%show("lb = -9; ub = +9")
                            lb = -9; ub = +9
            call disp%show("a1 = cmplx(getUnifRand(lb, ub), getUnifRand(lb, ub), CKG)")
                            a1 = cmplx(getUnifRand(lb, ub), getUnifRand(lb, ub), CKG)
            call disp%show("a2 = cmplx(getUnifRand(lb, ub), getUnifRand(lb, ub), CKG)")
                            a2 = cmplx(getUnifRand(lb, ub), getUnifRand(lb, ub), CKG)
            call disp%show("[a1, a2]")
            call disp%show( [a1, a2] )
            call disp%show("min(a1, a2)")
            call disp%show( min(a1, a2) )
            call disp%show("lenarr = getUnifRand(2, 5)")
                            lenarr = getUnifRand(2, 5)
            call disp%show("array1 = getUnifRand((-9._CKG, -9._CKG), (+9._CKG, +9._CKG), lenarr)")
                            array1 = getUnifRand((-9._CKG, -9._CKG), (+9._CKG, +9._CKG), lenarr)
            call disp%show("min(array1, a2)")
            call disp%show( min(array1, a2) )
            call disp%show("min(a2, array1)")
            call disp%show( min(a2, array1) )
            call disp%skip()
        end do
    end block

    block
        use pm_kind, only: CKG => CKD
        complex(CKG), allocatable :: array1(:)
        complex(CKG) :: a1, a2
        call disp%skip()
        do itry = 1, ntry
            call disp%show("lb = -9; ub = +9")
                            lb = -9; ub = +9
            call disp%show("a1 = cmplx(getUnifRand(lb, ub), getUnifRand(lb, ub), CKG)")
                            a1 = cmplx(getUnifRand(lb, ub), getUnifRand(lb, ub), CKG)
            call disp%show("a2 = cmplx(getUnifRand(lb, ub), getUnifRand(lb, ub), CKG)")
                            a2 = cmplx(getUnifRand(lb, ub), getUnifRand(lb, ub), CKG)
            call disp%show("[a1, a2]")
            call disp%show( [a1, a2] )
            call disp%show("min(a1, a2)")
            call disp%show( min(a1, a2) )
            call disp%show("lenarr = getUnifRand(2, 5)")
                            lenarr = getUnifRand(2, 5)
            call disp%show("array1 = getUnifRand((-9._CKG, -9._CKG), (+9._CKG, +9._CKG), lenarr)")
                            array1 = getUnifRand((-9._CKG, -9._CKG), (+9._CKG, +9._CKG), lenarr)
            call disp%show("min(array1, a2)")
            call disp%show( min(array1, a2) )
            call disp%show("min(a2, array1)")
            call disp%show( min(a2, array1) )
            call disp%skip()
        end do
    end block

    block
        use pm_kind, only: CKG => CKH
        complex(CKG) :: a1, a2
        complex(CKG), allocatable :: array1(:)
        call disp%skip()
        do itry = 1, ntry
            call disp%show("lb = -9; ub = +9")
                            lb = -9; ub = +9
            call disp%show("a1 = cmplx(getUnifRand(lb, ub), getUnifRand(lb, ub), CKG)")
                            a1 = cmplx(getUnifRand(lb, ub), getUnifRand(lb, ub), CKG)
            call disp%show("a2 = cmplx(getUnifRand(lb, ub), getUnifRand(lb, ub), CKG)")
                            a2 = cmplx(getUnifRand(lb, ub), getUnifRand(lb, ub), CKG)
            call disp%show("[a1, a2]")
            call disp%show( [a1, a2] )
            call disp%show("min(a1, a2)")
            call disp%show( min(a1, a2) )
            call disp%show("lenarr = getUnifRand(2, 5)")
                            lenarr = getUnifRand(2, 5)
            call disp%show("array1 = getUnifRand((-9._CKG, -9._CKG), (+9._CKG, +9._CKG), lenarr)")
                            array1 = getUnifRand((-9._CKG, -9._CKG), (+9._CKG, +9._CKG), lenarr)
            call disp%show("min(array1, a2)")
            call disp%show( min(array1, a2) )
            call disp%show("min(a2, array1)")
            call disp%show( min(a2, array1) )
            call disp%skip()
        end do
    end block

end program example