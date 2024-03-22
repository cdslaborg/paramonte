
program example

    use pm_io, only: display_type
    use pm_kind, only: IK, RKS, RKD, RKH
    use pm_distNorm, only: getNormFisher

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    ! Inverse variance.

    call disp%skip()
    call disp%show("getNormFisher(varInv = 1. / 3.3)")
    call disp%show( getNormFisher(varInv = 1. / 3.3) )
    call disp%skip()

    call disp%skip()
    call disp%show("getNormFisher(varInv = 1._RKS / 3.3_RKS)")
    call disp%show( getNormFisher(varInv = 1._RKS / 3.3_RKS) )
    call disp%skip()

    call disp%skip()
    call disp%show("getNormFisher(varInv = 1._RKD / 3.3_RKD)")
    call disp%show( getNormFisher(varInv = 1._RKD / 3.3_RKD) )
    call disp%skip()

    call disp%skip()
    call disp%show("getNormFisher(varInv = 1._RKH / 3.3_RKH)")
    call disp%show( getNormFisher(varInv = 1._RKH / 3.3_RKH) )
    call disp%skip()

end program example