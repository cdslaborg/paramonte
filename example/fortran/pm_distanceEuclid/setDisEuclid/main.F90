program example

    use pm_kind, only: SK, IK, LK, RKH
    use pm_kind, only: TKG => RKS ! all processor kinds are supported.
    use pm_io, only: display_type
    use pm_distanceEuclid, only: setDisEuclid, euclid, euclidu, euclidv, euclidsq
    use pm_ellipsoid, only: getVolUnitBall
    use pm_arrayResize, only: setResized
    use pm_distUnif, only: getUnifRand

    implicit none

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the distance of a point with respect to a reference in arbitrary dimensions without undue overflow/underflow.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKG) :: distance
        real(TKG), allocatable :: point(:), ref(:)
        call disp%skip()
        call disp%show("point = getUnifRand(1, 10, 3_IK)")
                        point = getUnifRand(1, 10, 3_IK)
        call disp%show("point")
        call disp%show( point )
        call disp%show("ref = getUnifRand(1, 10, 3_IK)")
                        ref = getUnifRand(1, 10, 3_IK)
        call disp%show("ref")
        call disp%show( ref )
        call disp%show("call setDisEuclid(distance, point, ref, euclidsq)")
                        call setDisEuclid(distance, point, ref, euclidsq)
        call disp%show("[norm2(point - ref), sqrt(distance), distance]")
        call disp%show( [norm2(point - ref), sqrt(distance), distance] )
        call disp%skip()
    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the asymmetric matrix of (squared) distances of a set of points from a set of reference points with or without undue overflow/underflow.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        real(TKG), allocatable :: point(:,:), ref(:,:), distance(:,:)
        integer(IK) :: ndim, npnt, nref
        call disp%skip()
        call disp%show("ndim = getUnifRand(1, 3); npnt = getUnifRand(1, 4); nref = getUnifRand(1, 3)")
                        ndim = getUnifRand(1, 3); npnt = getUnifRand(1, 4); nref = getUnifRand(1, 3)
        call disp%show("[ndim, npnt, nref]")
        call disp%show( [ndim, npnt, nref] )
        call disp%show("point = getUnifRand(1, 10, ndim, npnt)")
                        point = getUnifRand(1, 10, ndim, npnt)
        call disp%show("point")
        call disp%show( point )
        call disp%show("ref = getUnifRand(1, 10, ndim, nref)")
                        ref = getUnifRand(1, 10, ndim, nref)
        call disp%show("ref")
        call disp%show( ref )
        call disp%show("call setResized(distance, [npnt, nref])")
                        call setResized(distance, [npnt, nref])
        call disp%show("call setDisEuclid(distance, point, ref, euclidsq)")
                        call setDisEuclid(distance, point, ref, euclidsq)
        call disp%show("distance")
        call disp%show( distance )
        call disp%show("call setDisEuclid(distance(1,:), point(:,1), ref, euclidsq)")
                        call setDisEuclid(distance(1,:), point(:,1), ref, euclidsq)
        call disp%show("distance(1,:)")
        call disp%show( distance(1,:) )
        call disp%show("call setDisEuclid(distance(:,1), point, ref(:,1), euclidsq)")
                        call setDisEuclid(distance(:,1), point, ref(:,1), euclidsq)
        call disp%show("distance(:,1)")
        call disp%show( distance(:,1) )
        call disp%show("call setDisEuclid(distance(1,1), point(:,1), ref(:,1), euclidsq)")
                        call setDisEuclid(distance(1,1), point(:,1), ref(:,1), euclidsq)
        call disp%show("distance(1,1)")
        call disp%show( distance(1,1) )
        call disp%show("call setDisEuclid(distance(1,1), point(:,1), euclidsq) ! `ref` is the origin.")
                        call setDisEuclid(distance(1,1), point(:,1), euclidsq)
        call disp%show("distance(1,1)")
        call disp%show( distance(1,1) )
        call disp%show("call setDisEuclid(distance(:,1), point, euclidsq) ! `ref` is the origin.")
                        call setDisEuclid(distance(:,1), point, euclidsq)
        call disp%show("distance(:,1)")
        call disp%show( distance(:,1) )
        call disp%show("call setDisEuclid(distance, point(1,:), ref(1,:), euclid) ! 1D point and ref (faster than below).")
                        call setDisEuclid(distance, point(1,:), ref(1,:), euclid)
        call disp%show("distance")
        call disp%show( distance )
        call disp%show("call setDisEuclid(distance, point(1:1,:), ref(1:1,:), euclid) ! For comparison with the above.")
                        call setDisEuclid(distance, point(1:1,:), ref(1:1,:), euclid)
        call disp%show("distance")
        call disp%show( distance )
        call disp%show("call setDisEuclid(distance, point(1,:), ref(1,:), euclidsq) ! 1D point and ref (faster than below).")
                        call setDisEuclid(distance, point(1,:), ref(1,:), euclidsq)
        call disp%show("distance")
        call disp%show( distance )
        call disp%show("call setDisEuclid(distance, point(1:1,:), ref(1:1,:), euclidsq) ! For comparison with the above.")
                        call setDisEuclid(distance, point(1:1,:), ref(1:1,:), euclidsq)
        call disp%show("distance")
        call disp%show( distance )
        call disp%show("call setDisEuclid(distance, point(1,:), ref(1,:), euclidv) ! volume.")
                        call setDisEuclid(distance, point(1,:), ref(1,:), euclidv)
        call disp%show("distance")
        call disp%show( distance )
        call disp%show("call setDisEuclid(distance, point(1,:), ref(1,:), euclid) ! reference volume.")
                        call setDisEuclid(distance, point(1,:), ref(1,:), euclid)
        call disp%show("getVolUnitBall(1._TKG) * distance")
        call disp%show( getVolUnitBall(1._TKG) * distance )
        call disp%skip()
        call disp%show("call setDisEuclid(distance, point, ref, euclidv) ! volume")
                        call setDisEuclid(distance, point, ref, euclidv) ! volume
        call disp%show("distance")
        call disp%show( distance )
        call disp%show("call setDisEuclid(distance, point, ref, euclid) ! reference volume")
                        call setDisEuclid(distance, point, ref, euclid) ! reference volume
        call disp%show("getVolUnitBall(real(ndim, TKG)) * distance**ndim")
        call disp%show( getVolUnitBall(real(ndim, TKG)) * distance**ndim )
    end block

end program example