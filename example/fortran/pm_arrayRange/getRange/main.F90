program example

    use pm_kind, only: SK, IK, LK
    use pm_io, only: display_type
    use pm_kind, only: IK8, IK16, IK32, IK64
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
    call disp%show("getRange(0_IK8, 10_IK8)")
    call disp%show( getRange(0_IK8, 10_IK8) )
    call disp%skip
    call disp%show("getRange(0_IK8, 10_IK8, 2_IK8)")
    call disp%show( getRange(0_IK8, 10_IK8, 2_IK8) )
    call disp%skip
    call disp%show("getRange(0_IK8, 10_IK8, 3_IK8)")
    call disp%show( getRange(0_IK8, 10_IK8, 3_IK8) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getRange(10_IK8, 0_IK8, -2_IK8)")
    call disp%show( getRange(10_IK8, 0_IK8, -2_IK8) )
    call disp%skip
    call disp%show("getRange(10_IK8, 0_IK8, -3_IK8)")
    call disp%show( getRange(10_IK8, 0_IK8, -3_IK8) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!16-bit integer")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getRange(0_IK16, 10_IK16)")
    call disp%show( getRange(0_IK16, 10_IK16) )
    call disp%skip
    call disp%show("getRange(0_IK16, 10_IK16, 2_IK16)")
    call disp%show( getRange(0_IK16, 10_IK16, 2_IK16) )
    call disp%skip
    call disp%show("getRange(0_IK16, 10_IK16, 3_IK16)")
    call disp%show( getRange(0_IK16, 10_IK16, 3_IK16) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getRange(10_IK16, 0_IK16, -2_IK16)")
    call disp%show( getRange(10_IK16, 0_IK16, -2_IK16) )
    call disp%skip
    call disp%show("getRange(10_IK16, 0_IK16, -3_IK16)")
    call disp%show( getRange(10_IK16, 0_IK16, -3_IK16) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!32-bit integer")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getRange(0_IK32, 10_IK32)")
    call disp%show( getRange(0_IK32, 10_IK32) )
    call disp%skip
    call disp%show("getRange(0_IK32, 10_IK32, 2_IK32)")
    call disp%show( getRange(0_IK32, 10_IK32, 2_IK32) )
    call disp%skip
    call disp%show("getRange(0_IK32, 10_IK32, 3_IK32)")
    call disp%show( getRange(0_IK32, 10_IK32, 3_IK32) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getRange(10_IK32, 0_IK32, -2_IK32)")
    call disp%show( getRange(10_IK32, 0_IK32, -2_IK32) )
    call disp%skip
    call disp%show("getRange(10_IK32, 0_IK32, -3_IK32)")
    call disp%show( getRange(10_IK32, 0_IK32, -3_IK32) )
    call disp%skip

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%show("!64-bit integer")
    call disp%show("!%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip
    call disp%show("getRange(0_IK64, 10_IK64)")
    call disp%show( getRange(0_IK64, 10_IK64) )
    call disp%skip
    call disp%show("getRange(0_IK64, 10_IK64, 2_IK64)")
    call disp%show( getRange(0_IK64, 10_IK64, 2_IK64) )
    call disp%skip
    call disp%show("getRange(0_IK64, 10_IK64, 3_IK64)")
    call disp%show( getRange(0_IK64, 10_IK64, 3_IK64) )
    call disp%skip

    ! Generate sequence in reverse.

    call disp%skip
    call disp%show("getRange(10_IK64, 0_IK64, -2_IK64)")
    call disp%show( getRange(10_IK64, 0_IK64, -2_IK64) )
    call disp%skip
    call disp%show("getRange(10_IK64, 0_IK64, -3_IK64)")
    call disp%show( getRange(10_IK64, 0_IK64, -3_IK64) )
    call disp%skip

end program example