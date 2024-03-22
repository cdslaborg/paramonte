program example

    use pm_kind, only: SK, IK
    use pm_kind, only: RKC => RKS ! all processor real and complex kinds are supported.
    use pm_io, only: display_type
    use pm_polynomial, only: setPolyDiv
    use pm_polynomial, only: getPolyStr

    implicit none

    integer(IK) :: lenQuo
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("getPolyStr([real(RKC) :: ])")
    call disp%show( getPolyStr([real(RKC) :: ]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getPolyStr([real(RKC) :: -1.5])")
    call disp%show( getPolyStr([real(RKC) :: -1.5]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getPolyStr([real(RKC) :: 0., -2., 0., +4.])")
    call disp%show( getPolyStr([real(RKC) :: 0., -2., 0., +4.]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getPolyStr([real(RKC) :: -1., 2.2, -3.33, 4.444])")
    call disp%show( getPolyStr([real(RKC) :: -1., 2.2, -3.33, 4.444]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getPolyStr([complex(RKC) :: ])")
    call disp%show( getPolyStr([complex(RKC) :: ]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getPolyStr([complex(RKC) :: (-1.5, +1.5)])")
    call disp%show( getPolyStr([complex(RKC) :: (-1.5, +1.5)]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getPolyStr([complex(RKC) :: 0., -2., 0., +4.])")
    call disp%show( getPolyStr([complex(RKC) :: 0., -2., 0., +4.]) )
    call disp%skip()

    call disp%skip()
    call disp%show("getPolyStr([complex(RKC) :: -1., 2.2, -3.33, 4.444])")
    call disp%show( getPolyStr([complex(RKC) :: -1., 2.2, -3.33, 4.444]) )
    call disp%skip()

end program example