program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RKS ! all other real kinds are also supported.
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRand
    use pm_arrayResize, only: setResized
    use pm_clusKmeans, only: setMember
    use pm_arrayRange, only: getRange

    implicit none

    logical(LK) :: changed
    integer(IK) :: ndim, nsam, ncls
    real(RKG)   , allocatable  :: sample(:,:), center(:,:), disq(:)
    integer(IK) , allocatable  :: membership(:)
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute memberships of a sample of points from arbitrary dimensional cluster centers.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip

    call disp%skip()
    call disp%show("ndim = getUnifRand(1, 5); nsam = getUnifRand(1, 5); ncls = getUnifRand(1, 5);")
                    ndim = getUnifRand(1, 5); nsam = getUnifRand(1, 5); ncls = getUnifRand(1, 5);
    call disp%show("[ndim, nsam, ncls]")
    call disp%show( [ndim, nsam, ncls] )
    call disp%show("center = getUnifRand(0., 5., ndim, ncls) ! initialize random centers.")
                    center = getUnifRand(0., 5., ndim, ncls) ! initialize random centers.
    call disp%show("center")
    call disp%show( center )
    call disp%show("sample = getUnifRand(0., 5., ndim, nsam) ! Create a random sample.")
                    sample = getUnifRand(0., 5., ndim, nsam) ! Create a random sample.
    call disp%show("sample")
    call disp%show( sample )
    call disp%show("call setResized(disq, nsam)")
                    call setResized(disq, nsam)
    call disp%show("call setResized(membership, nsam)")
                    call setResized(membership, nsam)
    call disp%skip()

    call disp%show("call setMember(membership, disq, sample, center) ! sample points memberships.")
                    call setMember(membership, disq, sample, center) ! sample points memberships.
    call disp%show("membership")
    call disp%show( membership )
    call disp%show("disq")
    call disp%show( disq )
    call disp%skip()

    call disp%show("call setMember(membership, disq, sample, center, changed) ! sample points memberships.")
                    call setMember(membership, disq, sample, center, changed) ! sample points memberships.
    call disp%show("membership")
    call disp%show( membership )
    call disp%show("disq")
    call disp%show( disq )
    call disp%show("changed")
    call disp%show( changed )
    call disp%skip()

    call disp%show("call setMember(membership(1), disq(1), sample(:,1), center) ! single point membership.")
                    call setMember(membership(1), disq(1), sample(:,1), center) ! single point membership.
    call disp%show("membership")
    call disp%show( membership )
    call disp%show("disq")
    call disp%show( disq )
    call disp%skip()

    call disp%show("call setMember(membership(1), disq(1), sample(:,1), center, changed) ! single point membership.")
                    call setMember(membership(1), disq(1), sample(:,1), center, changed) ! single point membership.
    call disp%show("membership")
    call disp%show( membership )
    call disp%show("disq")
    call disp%show( disq )
    call disp%show("changed")
    call disp%show( changed )
    call disp%skip()

    call disp%show("call setMember(membership, disq, sample(1,:), center(1,:)) ! sample points memberships in one-dimension.")
                    call setMember(membership, disq, sample(1,:), center(1,:)) ! sample points memberships in one-dimension.
    call disp%show("membership")
    call disp%show( membership )
    call disp%show("disq")
    call disp%show( disq )
    call disp%skip()

    call disp%show("call setMember(membership, disq, sample(1,:), center(1,:), changed) ! sample points memberships in one-dimension.")
                    call setMember(membership, disq, sample(1,:), center(1,:), changed) ! sample points memberships in one-dimension.
    call disp%show("membership")
    call disp%show( membership )
    call disp%show("disq")
    call disp%show( disq )
    call disp%show("changed")
    call disp%show( changed )
    call disp%skip()

    call disp%show("call setMember(membership(1), disq(1), sample(1,1), center(1,:)) ! single point membership in one-dimension.")
                    call setMember(membership(1), disq(1), sample(1,1), center(1,:)) ! single point membership in one-dimension.
    call disp%show("membership")
    call disp%show( membership )
    call disp%show("disq")
    call disp%show( disq )
    call disp%skip()

    call disp%show("call setMember(membership(1), disq(1), sample(1,1), center(1,:), changed) ! single point membership in one-dimension.")
                    call setMember(membership(1), disq(1), sample(1,1), center(1,:), changed) ! single point membership in one-dimension.
    call disp%show("membership")
    call disp%show( membership )
    call disp%show("disq")
    call disp%show( disq )
    call disp%show("changed")
    call disp%show( changed )
    call disp%skip()

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
        call setResized(disq, nsam)
        call setResized(membership, nsam)
        call setMember(membership, disq, sample, center)
        call setMember(membership, disq, sample, center, changed)
        open(newunit = funit, file = "setMember.center.txt")
            do i = 1, ncls
                write(funit, "(*(g0,:,','))") i, center(:,i)
            end do
        close(funit)
        open(newunit = funit, file = "setMember.sample.txt")
            do i = 1, nsam
                write(funit, "(*(g0,:,','))") membership(i), sample(:,i)
            end do
        close(funit)
    end block

end program example