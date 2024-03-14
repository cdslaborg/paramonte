program example

    use pm_kind, only: SK, IK, LK, RKH
    use pm_polynomial, only: getPolyVal
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    block
        use pm_kind, only: TKC => RKH ! all processor kinds are supported.
        call disp%skip()
        call disp%show("getPolyVal([real(TKC) :: 1, 3, 0, 2], x = 2._TKC) ! 23.")
        call disp%show( getPolyVal([real(TKC) :: 1, 3, 0, 2], x = 2._TKC) )
        call disp%skip()
        call disp%show("getPolyVal([real(TKC) :: -19.0, 7.0, -4.0, 6.0], x = 3._TKC) ! 128.")
        call disp%show( getPolyVal([real(TKC) :: -19.0, 7.0, -4.0, 6.0], x = 3._TKC) )
        call disp%skip()
        call disp%show("getPolyVal([real(TKC) :: -19.0, 7.0, -4.0, 6.0], x = (3._TKC, 0._TKC)) ! (128., 0.)")
        call disp%show( getPolyVal([real(TKC) :: -19.0, 7.0, -4.0, 6.0], x = (3._TKC, 0._TKC)) )
        call disp%skip()
    end block

    block
        use pm_kind, only: TKC => RKS ! all processor kinds are supported.
        call disp%skip()
        call disp%show("getPolyVal([complex(TKC) :: (+2., +4.), (+24., +12.), (-28., +36.), (-15., +0.), (+24., -42.), (-7., +0.)], x = (2._TKC, -2._TKC))")
        call disp%show( getPolyVal([complex(TKC) :: (+2., +4.), (+24., +12.), (-28., +36.), (-15., +0.), (+24., -42.), (-7., +0.)], x = (2._TKC, -2._TKC)) )
        call disp%skip()
    end block

end program example