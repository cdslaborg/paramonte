program example

    use pm_kind, only: IKS, IKD
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
    call disp%show("getFactorial(5_IKS)")
    call disp%show( getFactorial(5_IKS) )
    call disp%skip()

    call disp%skip()
    call disp%show("getFactorial(5_IKD)")
    call disp%show( getFactorial(5_IKD) )
    call disp%skip()

    call disp%skip()
    call disp%show("getFactorial(int([(i, i = 1,13,2)],IKS))")
    call disp%show( getFactorial(int([(i, i = 1,13,2)],IKS)) )
    call disp%skip()

    call disp%skip()
    call disp%show("getFactorial(int([(i, i = 1,20,2)],IKD))")
    call disp%show( getFactorial(int([(i, i = 1,20,2)],IKD)) )
    call disp%skip()

end program example