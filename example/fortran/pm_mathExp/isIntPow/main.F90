program example

    use pm_kind, only: SK, IK, LK
    use pm_mathExp, only: isIntPow
    use pm_io, only: display_type

    implicit none

    integer(IK), allocatable :: expNext(:)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    allocate(expNext(1))

    call disp%skip()
    call disp%show("isIntPow(1)")
    call disp%show( isIntPow(1) )
    call disp%show("isIntPow(abs([1, -2, 3, -4, 5, 9, 519]))")
    call disp%show( isIntPow(abs([1, -2, 3, -4, 5, 9, 519])) )
    call disp%show("isIntPow(abs([1, -2, 3, -4, 5, 9, 519]), base = 3)")
    call disp%show( isIntPow(abs([1, -2, 3, -4, 5, 9, 519]), base = 3) )
    call disp%skip()

end program example