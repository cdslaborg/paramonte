program example

    use pm_kind, only: SK, IK, LK
    use pm_mathExp, only: getExpNext
    use pm_io, only: display_type

    implicit none

    integer(IK), allocatable :: expNext(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    allocate(expNext(1))

    call disp%skip()
    call disp%show("expNext(1) = getExpNext(.5)")
                    expNext(1) = getExpNext(.5)
    call disp%show("expNext(1)")
    call disp%show( expNext(1) )
    call disp%show("2.**expNext(1)")
    call disp%show( 2.**expNext(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("expNext(1) = getExpNext(1.5)")
                    expNext(1) = getExpNext(1.5)
    call disp%show("expNext(1)")
    call disp%show( expNext(1) )
    call disp%show("2.**expNext(1)")
    call disp%show( 2.**expNext(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("expNext(1) = getExpNext(1.)")
                    expNext(1) = getExpNext(1.)
    call disp%show("expNext(1)")
    call disp%show( expNext(1) )
    call disp%show("2.**expNext(1)")
    call disp%show( 2.**expNext(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("expNext(1) = getExpNext(1.1)")
                    expNext(1) = getExpNext(1.1)
    call disp%show("expNext(1)")
    call disp%show( expNext(1) )
    call disp%show("2.**expNext(1)")
    call disp%show( 2.**expNext(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("expNext = getExpNext(abs([1, -2, 3, -4, 5, 9, 519]))")
                    expNext = getExpNext(abs([1, -2, 3, -4, 5, 9, 519]))
    call disp%show("expNext")
    call disp%show( expNext )
    call disp%show("2.**expNext")
    call disp%show( 2.**expNext )
    call disp%skip()

    call disp%skip()
    call disp%show("expNext = getExpNext(abs([1, -2, 3, -4, 5, 9, 519]), base = 3)")
                    expNext = getExpNext(abs([1, -2, 3, -4, 5, 9, 519]), base = 3)
    call disp%show("expNext")
    call disp%show( expNext )
    call disp%show("3.**expNext")
    call disp%show( 3.**expNext )
    call disp%skip()

    call disp%skip()
    call disp%show("expNext = getExpNext(abs([real :: 1, -2, 3, -4, 5, 9, 519]), base = exp(1.))")
                    expNext = getExpNext(abs([real :: 1, -2, 3, -4, 5, 9, 519]), base = exp(1.))
    call disp%show("expNext")
    call disp%show( expNext )
    call disp%show("exp(1.)**expNext")
    call disp%show( exp(1.)**expNext )
    call disp%skip()

end program example