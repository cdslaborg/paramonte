program example

    use pm_kind, only: SK, IK, LK
    use pm_kind, only: RKG => RKS ! all other real kinds are also supported.
    use pm_io, only: display_type
    use pm_distUnif, only: getUnifRand
    use pm_arrayResize, only: setResized
    use pm_clusKmeans, only: setCenter
    use pm_clusKmeans, only: setMember

    implicit none

    integer(IK) :: ndim, nsam, ncls
    real(RKG)   , allocatable  :: sample(:,:), center(:,:), disq(:), potential(:)
    integer(IK) , allocatable  :: membership(:), size(:)
    type(display_type) :: disp

    disp = display_type(file = "main.out.F90")

    call disp%skip
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute cluster centers based on an input sample and cluster memberships and member-center distances.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
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
    call disp%show("call setResized(potential, ncls)")
                    call setResized(potential, ncls)
    call disp%show("call setResized(size, ncls)")
                    call setResized(size, ncls)
    call disp%show("call setMember(membership, disq, sample, center) ! get sample points memberships.")
                    call setMember(membership, disq, sample, center) ! get sample points memberships.
    call disp%skip()

    call disp%show("call setCenter(membership, disq, sample, center, size, potential) ! now compute the new clusters.")
                    call setCenter(membership, disq, sample, center, size, potential) ! now compute the new clusters.
    call disp%show("size")
    call disp%show( size )
    call disp%show("center")
    call disp%show( center )
    call disp%show("potential")
    call disp%show( potential )
    call disp%skip()

    call disp%show("call setCenter(membership, disq, sample(1,:), center(1,:), size, potential) ! sample points memberships in one-dimension.")
                    call setCenter(membership, disq, sample(1,:), center(1,:), size, potential) ! sample points memberships in one-dimension.
    call disp%show("size")
    call disp%show( size )
    call disp%show("center(1,:)")
    call disp%show( center(1,:) )
    call disp%show("potential")
    call disp%show( potential )
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
        call setResized(potential, ncls)
        call setResized(size, ncls)
        ! output the sample and the initial centers.
        call setMember(membership, disq, sample, center)
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
        call setCenter(membership, disq, sample, center, size, potential)
        open(newunit = funit, file = "setCenter.center.txt")
            do i = 1, ncls
                write(funit, "(*(g0,:,','))") i, center(:,i)
            end do
        close(funit)
        open(newunit = funit, file = "setCenter.sample.txt")
            do i = 1, nsam
                write(funit, "(*(g0,:,','))") membership(i), sample(:,i)
            end do
        close(funit)
    end block

end program example