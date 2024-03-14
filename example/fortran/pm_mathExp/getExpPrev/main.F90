program example

    use pm_kind, only: SK, IK, LK
    use pm_mathExp, only: getExpPrev
    use pm_io, only: display_type

    implicit none

    integer(IK), allocatable :: expPrev(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    allocate(expPrev(1))

    call disp%skip()
    call disp%show("expPrev(1) = getExpPrev(.5)")
                    expPrev(1) = getExpPrev(.5)
    call disp%show("expPrev(1)")
    call disp%show( expPrev(1) )
    call disp%show("2.**expPrev(1)")
    call disp%show( 2.**expPrev(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("expPrev(1) = getExpPrev(1.5)")
                    expPrev(1) = getExpPrev(1.5)
    call disp%show("expPrev(1)")
    call disp%show( expPrev(1) )
    call disp%show("2.**expPrev(1)")
    call disp%show( 2.**expPrev(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("expPrev(1) = getExpPrev(1.)")
                    expPrev(1) = getExpPrev(1.)
    call disp%show("expPrev(1)")
    call disp%show( expPrev(1) )
    call disp%show("2.**expPrev(1)")
    call disp%show( 2.**expPrev(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("expPrev(1) = getExpPrev(1.1)")
                    expPrev(1) = getExpPrev(1.1)
    call disp%show("expPrev(1)")
    call disp%show( expPrev(1) )
    call disp%show("2.**expPrev(1)")
    call disp%show( 2.**expPrev(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("expPrev = getExpPrev(abs([1, -2, 3, -4, 5, 9, 519]))")
                    expPrev = getExpPrev(abs([1, -2, 3, -4, 5, 9, 519]))
    call disp%show("expPrev")
    call disp%show( expPrev )
    call disp%show("2.**expPrev")
    call disp%show( 2.**expPrev )
    call disp%skip()

    call disp%skip()
    call disp%show("expPrev = getExpPrev(abs([1, -2, 3, -4, 5, 9, 519]), base = 3)")
                    expPrev = getExpPrev(abs([1, -2, 3, -4, 5, 9, 519]), base = 3)
    call disp%show("expPrev")
    call disp%show( expPrev )
    call disp%show("3.**expPrev")
    call disp%show( 3.**expPrev )
    call disp%skip()

    call disp%skip()
    call disp%show("expPrev = getExpPrev(abs([real :: 1, -2, 3, -4, 5, 9, 519]), base = exp(1.))")
                    expPrev = getExpPrev(abs([real :: 1, -2, 3, -4, 5, 9, 519]), base = exp(1.))
    call disp%show("expPrev")
    call disp%show( expPrev )
    call disp%show("exp(1.)**expPrev")
    call disp%show( exp(1.)**expPrev )
    call disp%skip()

end program example