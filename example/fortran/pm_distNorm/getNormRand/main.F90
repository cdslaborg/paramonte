program example

    use pm_kind, only: SK
    use pm_kind, only: IK, RK ! all real kinds are supported.
    use pm_distNorm, only: getNormRand
    use pm_arraySpace, only: setLinSpace
    use pm_arraySpace, only: setLogSpace
    use pm_io, only: display_type

    implicit none

    integer(IK), parameter  :: NP = 1000_IK
    real(RK), dimension(NP) :: mean, std, rand

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call setLinSpace(mean, x1 = -5._RK, x2 = +5._RK)
    call setLogSpace(std, logx1 = log(0.1_RK), logx2 = log(10._RK))

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate random numbers from the Normal distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Normal random number with a particular mean and default unity standard deviation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mean(1)")
    call disp%show( mean(1) )
    call disp%show("rand(1:3) = getNormRand(mean(1))")
                    rand(1:3) = getNormRand(mean(1))
    call disp%show("rand(1:3)")
    call disp%show( rand(1:3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Normal random number with given mean and standard deviation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mean(1)")
    call disp%show( mean(1) )
    call disp%show("std(1)")
    call disp%show( std(1) )
    call disp%show("rand(1:2) = getNormRand(mean(1), std(1))")
                    rand(1:2) = getNormRand(mean(1), std(1))
    call disp%show("rand(1)")
    call disp%show( rand(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Normal random numbers with a fixed set of mean and standard deviation.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mean(1:NP:NP/3)")
    call disp%show( mean(1:NP:NP/3) )
    call disp%show("std(1:NP:NP/3)")
    call disp%show( std(1:NP:NP/3) )
    call disp%show("rand(1:NP:NP/3) = getNormRand(mean(1), std(1))")
                    rand(1:NP:NP/3) = getNormRand(mean(1), std(1))
    call disp%show("rand(1:NP:NP/3)")
    call disp%show( rand(1:NP:NP/3) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Normal random number with a range of means and standard deviations.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("mean(1:NP:NP/3)")
    call disp%show( mean(1:NP:NP/3) )
    call disp%show("std(1:NP:NP/3)")
    call disp%show( std(1:NP:NP/3) )
    call disp%show("rand(1:NP:NP/3) = getNormRand(mean(1:NP:NP/3), std(1:NP:NP/3))")
                    rand(1:NP:NP/3) = getNormRand(mean(1:NP:NP/3), std(1:NP:NP/3))
    call disp%show("rand(1:NP:NP/3)")
    call disp%show( rand(1:NP:NP/3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example rand array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: fileUnit, i
        integer(IK), parameter :: NP = 5000_IK
        open(newunit = fileUnit, file = "getNormRand.RK.txt")
        write(fileUnit,"(3(g0,:,' '))") ( getNormRand(+2._RK, std = 3.0_RK) &
                                        , getNormRand(+0._RK, std = 1.0_RK) &
                                        , getNormRand(-5._RK) &
                                        , i = 1, NP &
                                        )
        close(fileUnit)
    end block

end program example