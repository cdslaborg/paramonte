program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_distUnifPar, only: setUnifParRand

    implicit none

    real :: rand(5)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("call random_number(rand(1:2))")
                    call random_number(rand(1:2))
    call disp%show("call setUnifParRand(rand(1:2), ub = 1.)")
                    call setUnifParRand(rand(1:2), ub = 1.)
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("call random_number(rand(1:2))")
                    call random_number(rand(1:2))
    call disp%show("call setUnifParRand(rand(1:2), lb = -100., ub = 1.)")
                    call setUnifParRand(rand(1:2), lb = -100., ub = 1.)
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("call random_number(rand(1:5))")
                    call random_number(rand(1:5))
    call disp%show("call setUnifParRand(rand(1:5), ub = [1., 2., 3., 4., 5.]) ! random vector of dimension 5.")
                    call setUnifParRand(rand(1:5), ub = [1., 2., 3., 4., 5.])
    call disp%show("rand(1:5)")
    call disp%show( rand(1:5) )
    call disp%skip()

    call disp%skip()
    call disp%show("call random_number(rand(1:5))")
                    call random_number(rand(1:5))
    call disp%show("call setUnifParRand(rand(1:5), lb = -[1., 2., 3., 4., 5.], ub = [1., 2., 3., 4., 5.]) ! random vector of dimension 5.")
                    call setUnifParRand(rand(1:5), lb = -[1., 2., 3., 4., 5.], ub = [1., 2., 3., 4., 5.])
    call disp%show("rand(1:5)")
    call disp%show( rand(1:5) )
    call disp%skip()

    call disp%skip()
    call disp%show("call random_number(rand(1:2))")
                    call random_number(rand(1:2))
    call disp%show("call setUnifParRand(rand(1:2), ub = reshape([1., 1., -1., -1.], shape = [2, 2]))")
                    call setUnifParRand(rand(1:2), ub = reshape([1., 1., -1., -1.], shape = [2, 2]))
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("call random_number(rand(1:2))")
                    call random_number(rand(1:2))
    call disp%show("call setUnifParRand(rand(1:2), lb = [10., 10.], ub = reshape([1., 1., -1., -1.], shape = [2, 2]))")
                    call setUnifParRand(rand(1:2), lb = [10., 10.], ub = reshape([1., 1., -1., -1.], shape = [2, 2]))
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        real :: rand(2, 4)
        integer(IK) :: fileUnit, i
        open(newunit = fileUnit, file = "setUnifParRand.RK.txt")
        do i = 1, 1000
            call random_number(rand)
            call setUnifParRand(rand(:,1), lb = -[2., 2.], ub = -[0., 1.])
            call setUnifParRand(rand(:,2), ub = reshape([1., 1., -1., +1.], shape = [2, 2]))
            call setUnifParRand(rand(:,3), lb = [2.0, -2.], ub = reshape([1., 1., +2., -2.], shape = [2, 2]))
            call setUnifParRand(rand(:,4), lb = [1.5, -.5], ub = reshape([2., 1., 1., 2.], shape = [2, 2]))
            write(fileUnit, "(*(g0,:,','))") rand
        end do
        close(fileUnit)
    end block

end program example