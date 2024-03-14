program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_distUnifPar, only: getUnifParRand

    implicit none

    real :: rand(5)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("rand(1:2) = getUnifParRand(ub = 1., ndim = 2)")
                    rand(1:2) = getUnifParRand(ub = 1., ndim = 2)
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("rand(1:2) = getUnifParRand(lb = -100., ub = 1., ndim = 2)")
                    rand(1:2) = getUnifParRand(lb = -100., ub = 1., ndim = 2)
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("rand(1:5) = getUnifParRand(ub = [1., 2., 3., 4., 5.]) ! random vector of dimension 5.")
                    rand(1:5) = getUnifParRand(ub = [1., 2., 3., 4., 5.])
    call disp%show("rand(1:5)")
    call disp%show( rand(1:5) )
    call disp%skip()

    call disp%skip()
    call disp%show("rand(1:5) = getUnifParRand(lb = -[1., 2., 3., 4., 5.], ub = [1., 2., 3., 4., 5.]) ! random vector of dimension 5.")
                    rand(1:5) = getUnifParRand(lb = -[1., 2., 3., 4., 5.], ub = [1., 2., 3., 4., 5.])
    call disp%show("rand(1:5)")
    call disp%show( rand(1:5) )
    call disp%skip()

    call disp%skip()
    call disp%show("rand(1:2) = getUnifParRand(ub = reshape([1., 1., -1., -1.], shape = [2, 2]))")
                    rand(1:2) = getUnifParRand(ub = reshape([1., 1., -1., -1.], shape = [2, 2]))
    call disp%show("rand(1:2)")
    call disp%show( rand(1:2) )
    call disp%skip()

    call disp%skip()
    call disp%show("rand(1:2) = getUnifParRand(lb = [10., 10.], ub = reshape([1., 1., -1., -1.], shape = [2, 2]))")
                    rand(1:2) = getUnifParRand(lb = [10., 10.], ub = reshape([1., 1., -1., -1.], shape = [2, 2]))
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
        open(newunit = fileUnit, file = "getUnifParRand.RK.txt")
        do i = 1, 1000
            rand(:,1) = getUnifParRand(lb = -[2., 2.], ub = -[0., 1.])
            rand(:,2) = getUnifParRand(ub = reshape([1., 1., -1., +1.], shape = [2, 2]))
            rand(:,3) = getUnifParRand(lb = [2.0, -2.], ub = reshape([1., 1., +2., -2.], shape = [2, 2]))
            rand(:,4) = getUnifParRand(lb = [1.5, -.5], ub = reshape([2., 1., 1., 2.], shape = [2, 2]))
            write(fileUnit, "(*(g0,:,','))") rand
        end do
        close(fileUnit)
    end block

end program example