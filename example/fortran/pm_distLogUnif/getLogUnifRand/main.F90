program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_distLogUnif, only: getLogUnifRand

    implicit none

    real :: rand(3)

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Generate random value(s) from the LogUniform distribution.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("rand(1) = getLogUnifRand(minx = 2., maxx = 5.)")
                    rand(1) = getLogUnifRand(minx = 2., maxx = 5.)
    call disp%show("rand(1)")
    call disp%show( rand(1) )
    call disp%skip()

    call disp%skip()
    call disp%show("rand(1:3) = getLogUnifRand(minx = [.3, .4, .25], maxx = 5.)")
                    rand(1:3) = getLogUnifRand(minx = [.3, .4, .25], maxx = 5.)
    call disp%show("rand(1:3)")
    call disp%show( rand(1:3) )
    call disp%skip()

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example array for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        use pm_arraySpace, only: setLinSpace
        real :: MinX(2), MaxX(2), rand(2)
        integer(IK) :: fileUnit, i
        MinX = [3., 2.0]
        MaxX = [7., 10.]
        open(newunit = fileUnit, file = "getLogUnifRand.RK.txt")
        do i = 1, 2000
            rand = getLogUnifRand(MinX, MaxX)
            write(fileUnit, "(*(g0,:,', '))") rand
        end do
        close(fileUnit)
    end block

end program example