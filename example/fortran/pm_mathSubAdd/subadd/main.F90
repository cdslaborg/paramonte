program example

    use pm_kind, only: SK, IK, RKH ! all intrinsic types and kinds are supported.
    use pm_io, only: display_type
    use pm_mathSubAdd, only: operator(.subadd.)

    implicit none

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show(".subadd. 5")
    call disp%show( .subadd. 5 )
    call disp%skip()

    call disp%skip()
    call disp%show("5 .subadd. 5")
    call disp%show( 5 .subadd. 5 )
    call disp%skip()

    call disp%skip()
    call disp%show(".subadd. 5.")
    call disp%show( .subadd. 5. )
    call disp%skip()

    call disp%skip()
    call disp%show("5.d0 .subadd. 5.d0")
    call disp%show( 5.d0 .subadd. 5.d0 )
    call disp%skip()

    call disp%skip()
    call disp%show(".subadd. (5., -5.)")
    call disp%show( .subadd. (5., -5.) )
    call disp%skip()

    call disp%skip()
    call disp%show("(5._RKH, -5._RKH) .subadd. (5._RKH, -5._RKH)")
    call disp%show( (5._RKH, -5._RKH) .subadd. (5._RKH, -5._RKH) )
    call disp%skip()

    call disp%skip()
    call disp%show(".subadd. huge(0)")
    call disp%show( .subadd. huge(0) )
    call disp%skip()

    call disp%skip()
    call disp%show(".subadd. huge(0.)")
    call disp%show( .subadd. huge(0.) )
    call disp%skip()

    call disp%skip()
    call disp%show(".subadd. cmplx(huge(0.), huge(0.))")
    call disp%show( .subadd. cmplx(huge(0.), huge(0.)) )
    call disp%skip()

end program example