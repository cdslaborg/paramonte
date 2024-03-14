program example

    use pm_kind, only: SK, IK
    use pm_kind, only: IK32, IK64 ! all processor types and kinds are supported.
    use pm_kind, only: RK32, RK64 ! all processor types and kinds are supported.
    use pm_kind, only: CK32, CK64 ! all processor types and kinds are supported.
    use pm_except, only: isAddOutflowPos
    use pm_io, only: display_type

    implicit none

    type(display_type) :: disp
    disp = display_type(file = SK_"main.out.F90")

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Integer addition outflow")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("isAddOutflowPos(1, 1)")
    call disp%show( isAddOutflowPos(1, 1) )
    call disp%skip

    call disp%skip
    call disp%show("isAddOutflowPos(1, -1)")
    call disp%show( isAddOutflowPos(1, -1) )
    call disp%skip

    call disp%skip
    call disp%show("isAddOutflowPos([huge(1), -huge(1)], [huge(1), -huge(1)])")
    call disp%show( isAddOutflowPos([huge(1), -huge(1)], [huge(1), -huge(1)]) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Complex addition outflow")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("isAddOutflowPos((1., -1), (1., -1))")
    call disp%show( isAddOutflowPos((1., -1), (1., -1)) )
    call disp%skip

    call disp%skip
    call disp%show("isAddOutflowPos(cmplx(1., huge(1.)), [(0., 0.), cmplx(huge(0.), 1.), cmplx(-huge(1), -huge(1))])")
    call disp%show( isAddOutflowPos(cmplx(1., huge(1.)), [(0., 0.), cmplx(huge(0.), 1.), cmplx(-huge(1), -huge(1))]) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!Real    addition outflow")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("isAddOutflowPos(1., 1.)")
    call disp%show( isAddOutflowPos(1., 1.) )
    call disp%skip

    call disp%skip
    call disp%show("isAddOutflowPos(1., -1.)")
    call disp%show( isAddOutflowPos(1., -1.) )
    call disp%skip

    call disp%skip
    call disp%show("isAddOutflowPos([huge(1.), huge(1.), huge(1.), -huge(1.)], [-1., 1., huge(1.), -huge(1.)])")
    call disp%show( isAddOutflowPos([huge(1.), huge(1.), huge(1.), -huge(1.)], [-1., 1., huge(1.), -huge(1.)]) )
    call disp%skip

end program example