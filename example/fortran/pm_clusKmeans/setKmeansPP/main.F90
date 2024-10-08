program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RKS ! all other real kinds are also supported.
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRand
    use pm_arrayResize, only: setResized
    use pm_clusKmeans, only: setKmeansPP, rngf

    implicit none

    integer(IK) :: ndim, nsam, ncls, itry
    real(RKG)   , allocatable  :: sample(:,:), center(:,:), disq(:), csdisq(:), potential(:)
    integer(IK) , allocatable  :: membership(:), size(:)
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute cluster centers based on an input sample and cluster memberships and member-center distances.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    do itry = 1, 10
    call disp%skip()
    call disp%show("ndim = getUnifRand(1, 5); ncls = getUnifRand(1, 5); nsam = getUnifRand(ncls, 2 * ncls);")
                    ndim = getUnifRand(1, 5); ncls = getUnifRand(1, 5); nsam = getUnifRand(ncls, 2 * ncls);
    call disp%show("[ndim, nsam, ncls]")
    call disp%show( [ndim, nsam, ncls] )
    call disp%show("sample = getUnifRand(0., 5., ndim, nsam) ! Create a random sample.")
                    sample = getUnifRand(0., 5., ndim, nsam) ! Create a random sample.
    call disp%show("sample")
    call disp%show( sample )
    call disp%show("call setResized(disq, nsam)")
                    call setResized(disq, nsam)
    call disp%show("call setResized(csdisq, nsam + 1_IK)")
                    call setResized(csdisq, nsam + 1_IK)
    call disp%show("call setResized(membership, nsam)")
                    call setResized(membership, nsam)
    call disp%show("call setResized(center, [ndim, ncls])")
                    call setResized(center, [ndim, ncls])
    call disp%show("call setResized(potential, ncls)")
                    call setResized(potential, ncls)
    call disp%show("call setResized(size, ncls)")
                    call setResized(size, ncls)
    call disp%skip()

    call disp%show("call setKmeansPP(rngf, membership, disq, csdisq, sample, center, size, potential) ! compute the new clusters and memberships.")
                    call setKmeansPP(rngf, membership, disq, csdisq, sample, center, size, potential) ! compute the new clusters and memberships.
    call disp%show("disq")
    call disp%show( disq )
    call disp%show("csdisq")
    call disp%show( csdisq )
    call disp%show("membership")
    call disp%show( membership )
    call disp%show("potential")
    call disp%show( potential )
    call disp%show("center")
    call disp%show( center )
    call disp%show("size")
    call disp%show( size )
    call disp%skip()
    end do

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output an example for visualization.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    block
        integer(IK) :: funit, i
        ndim = 2
        ncls = 5
        nsam = 5000
        center = getUnifRand(0., 1., ndim, ncls)
        sample = getUnifRand(0., 1., ndim, nsam)
        call setResized(csdisq, nsam + 1_IK)
        call setResized(disq, nsam)
        call setResized(membership, nsam)
        call setResized(center, [ndim, ncls])
        call setResized(potential, ncls)
        call setResized(size, ncls)
        call setKmeansPP(rngf, membership, disq, csdisq, sample, center, size, potential)
        open(newunit = funit, file = "setKmeansPP.center.txt")
            do i = 1, ncls
                write(funit, "(*(g0,:,','))") i, center(:,i)
            end do
        close(funit)
        open(newunit = funit, file = "setKmeansPP.sample.txt")
            do i = 1, nsam
                write(funit, "(*(g0,:,','))") membership(i), sample(:,i)
            end do
        close(funit)
    end block

end program example