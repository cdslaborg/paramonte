program example

    use pm_kind, only: SK, IK, LK, RKH
    use pm_kind, only: RKC => RKS ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distanceEuclid, only: getDisMatEuclid, rdpack, uppLow, uppLowDia, euclid, euclidu, euclidsq
    use pm_arrayResize, only: setResized
    use pm_distUnif, only: getUnifRand

    implicit none

    integer(IK) :: ndim, npnt, itry, ntry = 5
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the distance matrix of a set of points.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(RKC), allocatable :: distance(:,:), point(:,:)
        do itry = 1, ntry
            call disp%skip()
            call disp%show("ndim = getUnifRand(1, 3); npnt = getUnifRand(1, 7)")
                            ndim = getUnifRand(1, 3); npnt = getUnifRand(1, 7)
            call disp%show("[ndim, npnt]")
            call disp%show( [ndim, npnt] )
            call disp%show("point = getUnifRand(1, 10, ndim, npnt)")
                            point = getUnifRand(1, 10, ndim, npnt)
            call disp%show("point")
            call disp%show( point )
            call disp%show("distance = getDisMatEuclid(point)")
                            distance = getDisMatEuclid(point)
            call disp%show("distance")
            call disp%show( distance )
            call disp%show("distance = getDisMatEuclid(rdpack, uppLow, point) ! drop the zero-valued diagonal elements of the distance matrix.")
                            distance = getDisMatEuclid(rdpack, uppLow, point) ! drop the zero-valued diagonal elements of the distance matrix.
            call disp%show("distance")
            call disp%show( distance )
            call disp%show("distance = getDisMatEuclid(point, euclid)")
                            distance = getDisMatEuclid(point, euclid)
            call disp%show("distance")
            call disp%show( distance )
            call disp%show("distance = getDisMatEuclid(rdpack, uppLow, point, euclid) ! drop the zero-valued diagonal elements of the distance matrix.")
                            distance = getDisMatEuclid(rdpack, uppLow, point, euclid) ! drop the zero-valued diagonal elements of the distance matrix.
            call disp%show("distance")
            call disp%show( distance )
            call disp%show("distance = getDisMatEuclid(point, euclidu)")
                            distance = getDisMatEuclid(point, euclidu)
            call disp%show("distance")
            call disp%show( distance )
            call disp%show("distance = getDisMatEuclid(rdpack, uppLow, point, euclidsq) ! drop the zero-valued diagonal elements of the distance matrix.")
                            distance = getDisMatEuclid(rdpack, uppLow, point, euclidsq) ! drop the zero-valued diagonal elements of the distance matrix.
            call disp%show("distance")
            call disp%show( distance )
            call disp%show("distance = getDisMatEuclid(point, euclidsq)")
                            distance = getDisMatEuclid(point, euclidsq)
            call disp%show("distance")
            call disp%show( distance )
            call disp%show("distance = getDisMatEuclid(rdpack, uppLow, point, euclidsq) ! drop the zero-valued diagonal elements of the distance matrix.")
                            distance = getDisMatEuclid(rdpack, uppLow, point, euclidsq) ! drop the zero-valued diagonal elements of the distance matrix.
            call disp%show("distance")
            call disp%show( distance )
            call disp%skip()
        end do
    end block

end program example