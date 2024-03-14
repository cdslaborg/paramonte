program example

    use pm_kind, only: LK, SK
    use pm_complexAbs, only: abs
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("abs((-1., -2.))")
    call disp%show( abs((-1., -2.)) )
    call disp%skip()

    call disp%skip()
    call disp%show("abs([(-1., -2.), (-3., +2.), (+5., -2.), (+7., +1.)])")
    call disp%show( abs([(-1., -2.), (-3., +2.), (+5., -2.), (+7., +1.)]) )
    call disp%skip()

end program example