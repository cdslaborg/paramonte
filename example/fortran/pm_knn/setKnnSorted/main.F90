program example

    use pm_kind, only: SK, IK, LK, RKH
    use pm_io, only: display_type
    use pm_distanceEuclid, only: getDisMatEuclid, rdpack, uppLow, uppLowDia, euclid, euclidu, euclidsq
    use pm_arrayResize, only: setResized
    use pm_distUnif, only: getUnifRand
    use pm_knn, only: setKnnSorted

    implicit none

    integer(IK) :: npnt, nref, itry, ntry = 5
    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Sort a precomputed distance matrix according to the nearest neighbors.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block
        use pm_kind, only: RKG => RKS ! all processor kinds are supported.
        integer(IK), allocatable :: rank(:,:)
        real(RKG), allocatable :: dismat(:,:), point(:,:)
        do itry = 1, ntry
            call disp%skip()
            call disp%show("npnt = getUnifRand(1, 3); nref = getUnifRand(5, 7)")
                            npnt = getUnifRand(1, 3); nref = getUnifRand(5, 7)
            call disp%show("[npnt, nref]")
            call disp%show( [npnt, nref] )
            call disp%show("point = getUnifRand(1._RKG, 2._RKG, npnt, nref)")
                            point = getUnifRand(1._RKG, 2._RKG, npnt, nref)
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