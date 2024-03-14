program example

    use pm_kind, only: SK, IK, RKH ! all intrinsic types and kinds are supported.
    use pm_io, only: display_type
    use pm_mathDivMul, only: operator(.divmul.)

    implicit none

    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show(".divmul. 5")
    call disp%show( .divmul. 5 )
    call disp%skip()

    call disp%skip()
    call disp%show("5 .divmul. 5")
    call disp%show( 5 .divmul. 5 )
    call disp%skip()

    call disp%skip()
    call disp%show(".divmul. 5.")
    call disp%show( .divmul. 5. )
    call disp%skip()

    call disp%skip()
    call disp%show("5.d0 .divmul. 5.d0")
    call disp%show( 5.d0 .divmul. 5.d0 )
    call disp%skip()

    call disp%skip()
    call disp%show(".divmul. (5., -5.)")
    call disp%show( .divmul. (5., -5.) )
    call disp%skip()

    call disp%skip()
    call disp%show("(5._RKH, -5._RKH) .divmul. (5._RKH, -5._RKH)")
    call disp%show( (5._RKH, -5._RKH) .divmul. (5._RKH, -5._RKH) )
    call disp%skip()

    call disp%skip()
    call disp%show(".divmul. huge(0)/2")
    call disp%show( .divmul. huge(0)/2 )
    call disp%skip()

    call disp%skip()
    call disp%show(".divmul. huge(0.)/2")
    call disp%show( .divmul. huge(0.)/2 )
    call disp%skip()

    call disp%skip()
    call disp%show(".divmul. cmplx(huge(0.), huge(0.))/2")
    call disp%show( .divmul. cmplx(huge(0.), huge(0.))/2 )
    call disp%skip()

end program example