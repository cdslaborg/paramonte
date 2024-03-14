program example

    use pm_kind, only: SK, IK
    use pm_io, only: display_type
    use pm_distExp, only: setExpRand
    use pm_distUnif, only: setUnifRand

    implicit none

    integer(IK) , parameter :: NP = 5_IK
    real :: rand(NP,NP)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")


    call disp%show("call setUnifRand(rand(1,1))")
                    call setUnifRand(rand(1,1))
    call disp%show("call setExpRand(rand(1,1))")
                    call setExpRand(rand(1,1))
    call disp%show("rand(1,1)")
    call disp%show( rand(1,1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a scalar real random number.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setUnifRand(rand(1,1))")
                    call setUnifRand(rand(1,1))
    call disp%show("call setExpRand(rand(1,1), sigma = 2.) ! random real number with the Exponential mean 2..")
                    call setExpRand(rand(1,1), sigma = 2.)
    call disp%show("rand(1,1)")
    call disp%show( rand(1,1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a vector of real random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setUnifRand(rand(:,1))")
                    call setUnifRand(rand(:,1))
    call disp%show("call setExpRand(rand(:,1)) ! random real vector with the default Exponential mean (1,1).")
                    call setExpRand(rand(:,1))
    call disp%show("rand(:,1)")
    call disp%show( rand(:,1) )
    call disp%skip()

    call disp%show("call setUnifRand(rand(:,1))")
                    call setUnifRand(rand(:,1))
    call disp%show("call setExpRand(rand(:,1), sigma = 10.) ! random real vector with mean 10..")
                    call setExpRand(rand(:,1), sigma = 10.)
    call disp%show("rand(:,1)")
    call disp%show( rand(:,1) )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate a matrix of real random numbers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%show("call setUnifRand(rand)")
                    call setUnifRand(rand)
    call disp%show("call setExpRand(rand) ! random matrix with the default Exponential mean 1..")
                    call setExpRand(rand)
    call disp%show("rand")
    call disp%show( rand )
    call disp%skip()

    call disp%show("call setUnifRand(rand)")
                    call setUnifRand(rand)
    call disp%show("call setExpRand(rand, sigma = 100.) ! random real vector with mean 100..")
                    call setExpRand(rand, sigma = 100.)
    call disp%show("rand")
    call disp%show( rand )
    call disp%skip()

    block
        integer :: fileUnit, i
        real :: rand(2)
        open(newunit = fileUnit, file = "setExpRand.RK.txt")
        do i = 1, 1000
            call setUnifRand(rand)
            call setExpRand(rand, sigma = [3., 2.], mu = [0., 2.])
            write(fileUnit,"(2(g0,:,', '))") rand
        end do
        close(fileUnit)
    end block

end program example