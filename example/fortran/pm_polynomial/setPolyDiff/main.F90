program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRand
    use pm_polynomial, only: getPolyStr
    use pm_polynomial, only: setPolyDiff
    use pm_arrayResize, only: setResized

    implicit none

    integer(IK) :: order
    integer(IK) :: degree
    type(display_type) :: disp
    integer(IK) :: itry, ntry = 20
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKC => RKS ! all processor real and complex kinds are supported.
        real(TKC), allocatable :: coef(:), diff(:)
        do itry = 1, ntry
            call disp%show("degree = getUnifRand(0, 9_IK)")
                            degree = getUnifRand(0, 9_IK)
            call disp%show("degree")
            call disp%show( degree )
            call disp%show("coef = getUnifRand(-9, 9, degree)")
                            coef = getUnifRand(-9, 9, degree)
            call disp%show("coef")
            call disp%show( coef )
            call disp%show("getPolyStr(coef)")
            call disp%show( getPolyStr(coef) )
            call disp%show("order = getUnifRand(0, size(coef) + 1)")
                            order = getUnifRand(0, size(coef) + 1)
            call disp%show("order")
            call disp%show( order )
            call disp%skip()
            call disp%show("call setResized(diff, max(0_IK, size(coef, 1, IK) - order))")
                            call setResized(diff, max(0_IK, size(coef, 1, IK) - order))
            call disp%show("call setPolyDiff(diff, coef, order)")
                            call setPolyDiff(diff, coef, order)
            call disp%show("diff ! derivative coefficients.")
            call disp%show( diff )
            call disp%show("getPolyStr(diff)")
            call disp%show( getPolyStr(diff) )
            call disp%skip()
        end do
    end block

    block
        use pm_kind, only: TKC => RKS ! all processor real and complex kinds are supported.
        complex(TKC), allocatable :: coef(:), diff(:)
        do itry = 1, ntry
            call disp%show("degree = getUnifRand(0, 9_IK)")
                            degree = getUnifRand(0, 9_IK)
            call disp%show("degree")
            call disp%show( degree )
            call disp%show("coef = cmplx(getUnifRand(-9, 9, degree), getUnifRand(-9, 9, degree), TKC)")
                            coef = cmplx(getUnifRand(-9, 9, degree), getUnifRand(-9, 9, degree), TKC)
            call disp%show("coef")
            call disp%show( coef )
            call disp%show("getPolyStr(coef)")
            call disp%show( getPolyStr(coef) )
            call disp%show("order = getUnifRand(0, size(coef) + 1)")
                            order = getUnifRand(0, size(coef) + 1)
            call disp%show("order")
            call disp%show( order )
            call disp%skip()
            call disp%show("call setResized(diff, max(0_IK, size(coef, 1, IK) - order))")
                            call setResized(diff, max(0_IK, size(coef, 1, IK) - order))
            call disp%show("call setPolyDiff(diff, coef, order)")
                            call setPolyDiff(diff, coef, order)
            call disp%show("diff ! derivative coefficients.")
            call disp%show( diff )
            call disp%show("getPolyStr(diff)")
            call disp%show( getPolyStr(diff) )
            call disp%skip()
        end do
    end block

end program example