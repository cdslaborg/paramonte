program example

    use pm_kind, only: SK, IK, LK, RKH
    use pm_io, only: display_type
    use pm_distanceEuclid, only: getDisMatEuclid, rdpack, uppLow, uppLowDia, euclid, euclidu, euclidsq
    use pm_arrayResize, only: setResized
    use pm_distUnif, only: getUnifRand
    use pm_knn, only: setKnnSorted

    implicit none

    integer(IK) :: nref, npnt, itry, ntry = 5
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Sort a precomputed distance matrix according to the nearest neighbors.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: RKC => RKS ! all processor kinds are supported.
        integer(IK), allocatable :: rank(:,:)
        real(RKC), allocatable :: dismat(:,:), point(:,:)
        do itry = 1, ntry
            call disp%skip()
            call disp%show("nref = getUnifRand(1, 3); npnt = getUnifRand(5, 7)")
                            nref = getUnifRand(1, 3); npnt = getUnifRand(5, 7)
            call disp%show("[nref, npnt]")
            call disp%show( [nref, npnt] )
            call disp%show("point = getUnifRand(1._RKC, 2._RKC, nref, npnt)")
                            point = getUnifRand(1._RKC, 2._RKC, nref, npnt)
            call disp%show("point")
            call disp%show( point )
            call disp%show("dismat = getDisMatEuclid(point)")
                            dismat = getDisMatEuclid(point)
            call disp%show("dismat")
            call disp%show( dismat )
            call disp%show("call setKnnSorted(dismat, k = 1)")
                            call setKnnSorted(dismat, k = 1)
            call disp%show("dismat")
            call disp%show( dismat )
            call disp%show("call setKnnSorted(dismat, k = 2)")
                            call setKnnSorted(dismat, k = 2)
            call disp%show("dismat")
            call disp%show( dismat )
            call disp%show("dismat = getDisMatEuclid(point)")
                            dismat = getDisMatEuclid(point)
            call disp%show("dismat")
            call disp%show( dismat )
            call disp%show("call setResized(rank, shape(dismat))")
                            call setResized(rank, shape(dismat))
            call disp%show("call setKnnSorted(dismat, rank)")
                            call setKnnSorted(dismat, rank)
            call disp%show("dismat")
            call disp%show( dismat )
            call disp%show("rank")
            call disp%show( rank )
            call disp%skip()
            call disp%show("dismat = getDisMatEuclid(rdpack, uppLow, point) ! drop the zero-valued diagonal elements of the dismat matrix.")
                            dismat = getDisMatEuclid(rdpack, uppLow, point) ! drop the zero-valued diagonal elements of the dismat matrix.
            call disp%show("dismat")
            call disp%show( dismat )
            call disp%show("call setKnnSorted(dismat, k = 1)")
                            call setKnnSorted(dismat, k = 1)
            call disp%show("dismat")
            call disp%show( dismat )
            call disp%show("dismat = getDisMatEuclid(rdpack, uppLow, point) ! drop the zero-valued diagonal elements of the dismat matrix.")
                            dismat = getDisMatEuclid(rdpack, uppLow, point) ! drop the zero-valued diagonal elements of the dismat matrix.
            call disp%show("dismat")
            call disp%show( dismat )
            call disp%show("call setResized(rank, shape(dismat))")
                            call setResized(rank, shape(dismat))
            call disp%show("call setKnnSorted(dismat, rank)")
                            call setKnnSorted(dismat, rank)
            call disp%show("dismat")
            call disp%show( dismat )
            call disp%show("rank")
            call disp%show( rank )
            call disp%skip()
        end do
    end block

end program example