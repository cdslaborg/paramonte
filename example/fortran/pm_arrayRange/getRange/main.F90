program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_kind, only: IKL, IKS, IKD, IKH
    use pm_arrayRange, only: getRange

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("!%%%%%%%%%%")
    call disp%show("!real range")
    call disp%show("!%%%%%%%%%%")
    call disp%skip

    block
        real, allocatable :: range(:) ! all non-default kinds are also supported.

        call disp%skip
        call disp%show("range = getRange(0., 10 * nearest(0., 1.))")
                        range = getRange(0., 10 * nearest(0., 1.))
        call disp%show("range")
        call disp%show( range , deliml = """" )
        call disp%skip

        call disp%skip
        call disp%show("range = getRange(1., 1. + 10 * epsilon(1.))")
                        range = getRange(1., 1. + 10 * epsilon(1.))
        call disp%show("range")
        call disp%show( range , deliml = """" )
        call disp%skip

        call disp%skip
        call disp%show("range = getRange(1., 1. - 10 * epsilon(1.))")
                        range = getRange(1., 1. - 10 * epsilon(1.))
        call disp%show("range")
        call disp%show( range , deliml = """" )
        call disp%skip

        call disp%skip
        call disp%show("range = getRange(0., 9., 2.)")
                        range = getRange(0., 9., 2.)
        call disp%show("range")
        call disp%show( range , deliml = """" )
        call disp%skip

        call disp%skip
        call disp%show("range = getRange(9., 0., -2.)")
                        range = getRange(9., 0., -2.)
        call disp%show("range")
        call disp%show( range , deliml = """" )
        call disp%skip

    end block

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%show("!character range")
    call disp%show("!%%%%%%%%%%%%%%%")
    call disp%skip

    block
        character(:), allocatable :: range ! all non-default kinds are also supported.

        call disp%skip
        call disp%show("range = getRange('A', 'Z')")
                        range = getRange('A', 'Z')
        call disp%show("range")
        call disp%show( range , deliml = """" )
        call disp%skip

        call disp%skip
        call disp%show("range = getRange('A', 'Z', 2)")
                        range = getRange('A', 'Z', 2)
        call disp%show("range")
        call disp%show( range , deliml = """" )
        call disp%skip

        call disp%skip
        call disp%show("range = getRange('A', 'Z', 3)")
                        range = getRange('A', 'Z', 3)
        call disp%show("range")
        call disp%show( range , deliml = """" )
        call disp%skip

        call disp%skip
        call disp%show("range = getRange('A', 'z', 3)")
                        range = getRange('A', 'z', 3)
        call disp%show("range")
        call disp%show( range , deliml = """" )
        call disp%skip

        call disp%skip
        call disp%show("range = getRange('Z', 'A', -2)")
                        range = getRange('Z', 'A', -2)
        call disp%show("range")
        call disp%show( range , deliml = """" )
        call disp%skip
    end block

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%")
    call disp%show("!8-bit integer")
    call disp%show("!%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getRange(0_IKL, 10_IKL)")
    call disp%show( getRange(0_IKL, 10_IKL) )
    call disp%skip
    call disp%show("getRange(0_IKL, 10_IKL, 2_IKL)")
    call disp%show( getRange(0_IKL, 10_IKL, 2_IKL) )
    call disp%skip
    call disp%show("getRange(0_IKL, 10_IKL, 3_IKL)")
    call disp%show( getRange(0_IKL, 10_IKL, 3_IKL) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getRange(10_IKL, 0_IKL, -2_IKL)")
    call disp%show( getRange(10_IKL, 0_IKL, -2_IKL) )
    call disp%skip
    call disp%show("getRange(10_IKL, 0_IKL, -3_IKL)")
    call disp%show( getRange(10_IKL, 0_IKL, -3_IKL) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!16-bit integer")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getRange(0_IKS, 10_IKS)")
    call disp%show( getRange(0_IKS, 10_IKS) )
    call disp%skip
    call disp%show("getRange(0_IKS, 10_IKS, 2_IKS)")
    call disp%show( getRange(0_IKS, 10_IKS, 2_IKS) )
    call disp%skip
    call disp%show("getRange(0_IKS, 10_IKS, 3_IKS)")
    call disp%show( getRange(0_IKS, 10_IKS, 3_IKS) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getRange(10_IKS, 0_IKS, -2_IKS)")
    call disp%show( getRange(10_IKS, 0_IKS, -2_IKS) )
    call disp%skip
    call disp%show("getRange(10_IKS, 0_IKS, -3_IKS)")
    call disp%show( getRange(10_IKS, 0_IKS, -3_IKS) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!32-bit integer")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getRange(0_IKD, 10_IKD)")
    call disp%show( getRange(0_IKD, 10_IKD) )
    call disp%skip
    call disp%show("getRange(0_IKD, 10_IKD, 2_IKD)")
    call disp%show( getRange(0_IKD, 10_IKD, 2_IKD) )
    call disp%skip
    call disp%show("getRange(0_IKD, 10_IKD, 3_IKD)")
    call disp%show( getRange(0_IKD, 10_IKD, 3_IKD) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getRange(10_IKD, 0_IKD, -2_IKD)")
    call disp%show( getRange(10_IKD, 0_IKD, -2_IKD) )
    call disp%skip
    call disp%show("getRange(10_IKD, 0_IKD, -3_IKD)")
    call disp%show( getRange(10_IKD, 0_IKD, -3_IKD) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!64-bit integer")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getRange(0_IKH, 10_IKH)")
    call disp%show( getRange(0_IKH, 10_IKH) )
    call disp%skip
    call disp%show("getRange(0_IKH, 10_IKH, 2_IKH)")
    call disp%show( getRange(0_IKH, 10_IKH, 2_IKH) )
    call disp%skip
    call disp%show("getRange(0_IKH, 10_IKH, 3_IKH)")
    call disp%show( getRange(0_IKH, 10_IKH, 3_IKH) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getRange(10_IKH, 0_IKH, -2_IKH)")
    call disp%show( getRange(10_IKH, 0_IKH, -2_IKH) )
    call disp%skip
    call disp%show("getRange(10_IKH, 0_IKH, -3_IKH)")
    call disp%show( getRange(10_IKH, 0_IKH, -3_IKH) )
    call disp%skip

end program example