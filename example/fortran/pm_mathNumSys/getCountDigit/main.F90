program example

    use pm_kind, only: SK, IK
    use pm_kind, only: IKS, IKD ! all integer kinds are supported.
    use pm_io, only: display_type
    use pm_mathNumSys, only: getCountDigit

    implicit none

    integer(IKD), allocatable :: IntList(:)

    integer(IK) :: i
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the number of digits in a given integer.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("getCountDigit(5_IKS)")
    call disp%show( getCountDigit(5_IKS) )
    call disp%skip()

    call disp%skip()
    call disp%show("getCountDigit(5_IKD)")
    call disp%show( getCountDigit(5_IKD) )
    call disp%skip()

    call disp%skip()
    call disp%show("[-huge(0_IKS), huge(0_IKS)]")
    call disp%show( [-huge(0_IKS), huge(0_IKS)] )
    call disp%show("getCountDigit([-huge(0_IKS), huge(0_IKS)])")
    call disp%show( getCountDigit([-huge(0_IKS), huge(0_IKS)]) )
    call disp%skip()

    call disp%skip()
    call disp%show("[-huge(0_IKD), huge(0_IKD)]")
    call disp%show( [-huge(0_IKD), huge(0_IKD)] )
    call disp%show("getCountDigit([-huge(0_IKD), huge(0_IKD)])")
    call disp%show( getCountDigit([-huge(0_IKD), huge(0_IKD)]) )
    call disp%skip()

    call disp%skip()
    call disp%show("IntList = [(2_IKD**i, i = 0, 63, 16)]")
                    IntList = [(2_IKD**i, i = 0, 63, 16)]
    call disp%show("IntList")
    call disp%show( IntList )
    call disp%show("getCountDigit(IntList)")
    call disp%show( getCountDigit(IntList) )
    call disp%skip()

end program example