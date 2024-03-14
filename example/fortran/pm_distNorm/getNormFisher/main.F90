
program example

    use pm_io, only: display_type
    use pm_kind, only: IK, RK32, RK64, RK128
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
    call disp%show("getNormFisher(varInv = 1._RK32 / 3.3_RK32)")
    call disp%show( getNormFisher(varInv = 1._RK32 / 3.3_RK32) )
    call disp%skip()

    call disp%skip()
    call disp%show("getNormFisher(varInv = 1._RK64 / 3.3_RK64)")
    call disp%show( getNormFisher(varInv = 1._RK64 / 3.3_RK64) )
    call disp%skip()

    call disp%skip()
    call disp%show("getNormFisher(varInv = 1._RK128 / 3.3_RK128)")
    call disp%show( getNormFisher(varInv = 1._RK128 / 3.3_RK128) )
    call disp%skip()

end program example