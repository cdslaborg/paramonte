program example

    use pm_kind, only: IK32, IK64
    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_mathFactorial, only: getFactorial

    implicit none

    integer(IK) :: i

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the factorial of a scalar or array of integers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getFactorial(5_IK32)")
    call disp%show( getFactorial(5_IK32) )
    call disp%skip()

    call disp%skip()
    call disp%show("getFactorial(5_IK64)")
    call disp%show( getFactorial(5_IK64) )
    call disp%skip()

    call disp%skip()
    call disp%show("getFactorial(int([(i, i = 1,13,2)],IK32))")
    call disp%show( getFactorial(int([(i, i = 1,13,2)],IK32)) )
    call disp%skip()

    call disp%skip()
    call disp%show("getFactorial(int([(i, i = 1,20,2)],IK64))")
    call disp%show( getFactorial(int([(i, i = 1,20,2)],IK64)) )
    call disp%skip()

end program example