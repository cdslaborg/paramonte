program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_distExp, only: getExpRand

    implicit none

    integer :: i
    integer(IK) , parameter :: NP = 5_IK
    real :: rand(NP,NP)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")


    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a scalar real random number.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("rand(1,1) = getExpRand(sigma = 2.) ! random real number with the Exponential mean 2..")
                    rand(1,1) = getExpRand(sigma = 2.)
    call disp%show("rand(1,1)")
    call disp%show( rand(1,1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a vector of real random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("rand(:,1) = getExpRand(sigma = [(10., i = 1,NP)]) ! random real vector with the default Exponential mean (1,1).")
                    rand(:,1) = getExpRand(sigma = [(10., i = 1,NP)])
    call disp%show("rand(:,1)")
    call disp%show( rand(:,1) )
    call disp%skip()

    block
        integer :: fileUnit, i
        open(newunit = fileUnit, file = "getExpRand.RK.txt")
        do i = 1, 1000
            write(fileUnit,"(2(g0,:,', '))") getExpRand(sigma = [3., 2.], mu = [0., 2.])
        end do
        close(fileUnit)
    end block

end program example