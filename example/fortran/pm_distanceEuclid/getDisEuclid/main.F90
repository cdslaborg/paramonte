program example

    use pm_kind, only: SK, IK, LK, RKH
    use pm_kind, only: RKG => RKS ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distanceEuclid, only: getDisEuclid, euclid, euclidu, euclidsq

    implicit none

    type(display_type) :: disp
    complex(RKG)                :: z
    real(RKG)   , allocatable   :: point(:)
    real(RKG)   , allocatable   :: ref(:)
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("point = [real(RKG) :: 1, 3, 2]")
                    point = [real(RKG) :: 1, 3, 2]
    call disp%show("point")
    call disp%show( point )
    call disp%show("[getDisEuclid(point), getDisEuclid(point(1), point(2), point(3))]") ! norm2(point), 
    call disp%show( [getDisEuclid(point), getDisEuclid(point(1), point(2), point(3))] ) ! norm2(point), 
    call disp%skip()

    call disp%skip()
    call disp%show("point = [real(RKG) :: 1, 1, 1] * huge(0._RKG) / 10")
                    point = [real(RKG) :: 1, 1, 1] * huge(0._RKG) / 10
    call disp%show("point")
    call disp%show( point )
    call disp%show("[getDisEuclid(point), getDisEuclid(point(1), point(2), point(3))]   !, sqrt(dot_product(point, point)) norm2(point), ! Note that GNU gfortran `norm2()` respects robustness, while Intel ifort does not.")
    call disp%show( [getDisEuclid(point), getDisEuclid(point(1), point(2), point(3))] ) !, sqrt(dot_product(point, point))
    call disp%skip()

    call disp%skip()
    call disp%show("z = (1._RKG, -1._RKG) * huge(0._RKG) / 10")
                    z = (1._RKG, -1._RKG) * huge(0._RKG) / 10
    call disp%show("point")
    call disp%show( point )
    call disp%show("[hypot(z%re, z%im), abs(z), getDisEuclid([z%re, z%im])] ! norm2([z%re, z%im]), Note that GNU gfortran `norm2()` respects robustness, while Intel ifort does not.")
    call disp%show( [hypot(z%re, z%im), abs(z), getDisEuclid([z%re, z%im])] )
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the distance of a point with respect to a reference in arbitrary dimensions without undue overflow/underflow.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    call disp%skip()
    call disp%show("point = [real(RKG) :: 1, 0, 5, 2]")
                    point = [real(RKG) :: 1, 0, 5, 2]
    call disp%show("ref = [real(RKG) :: 2, 1, 3, 0]")
                    ref = [real(RKG) :: 2, 1, 3, 0]
    call disp%show("point")
    call disp%show( point )
    call disp%show("[getDisEuclid(point - ref), getDisEuclid(point, ref)]") ! norm2(point - ref), 
    call disp%show( [getDisEuclid(point - ref), getDisEuclid(point, ref)] ) ! norm2(point - ref), 
    call disp%skip()

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the asymmetric matrix of (squared) distances of a set of points from a set of reference points with or without undue overflow/underflow.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(RKG), allocatable :: point(:,:), ref(:,:), distance(:,:)
        call disp%skip()
        call disp%show("point = reshape([real(RKG) :: 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1], shape = [3, 5])")
                        point = reshape([real(RKG) :: 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1], shape = [3, 5])
        call disp%show("point")
        call disp%show( point )
        call disp%show("ref = reshape([real(RKG) :: -1, -1, -1, 1, 1, 1], shape = [3, 2]) ! two reference points.")
                        ref = reshape([real(RKG) :: -1, -1, -1, 1, 1, 1], shape = [3, 2])
        call disp%show("ref")
        call disp%show( ref )
        call disp%show("distance = getDisEuclid(point, ref)")
                        distance = getDisEuclid(point, ref)
        call disp%show("distance")
        call disp%show( distance )
        call disp%show("distance = getDisEuclid(point, ref, euclid)")
                        distance = getDisEuclid(point, ref, euclid)
        call disp%show("distance")
        call disp%show( distance )
        call disp%show("distance = getDisEuclid(point, ref, euclidu)")
                        distance = getDisEuclid(point, ref, euclidu)
        call disp%show("distance")
        call disp%show( distance )
        call disp%show("distance = getDisEuclid(point, ref, euclidsq)")
                        distance = getDisEuclid(point, ref, euclidsq)
        call disp%show("distance")
        call disp%show( distance )
        call disp%skip()
    end block

end program example