program example

    use pm_kind, only: SK, IK
    use pm_arrayFill, only: getFilled
    use pm_distUnif, only: getUnifRand
    use pm_arrayResize, only: setResized
    use pm_distCov, only: getCovRand
    use pm_distanceMahal, only: getMahalSq
    use pm_sampleShift, only: getShifted
    use pm_matrixInit, only: getMatInit
    use pm_matrixInit, only: uppLowDia
    use pm_io, only: display_type

    implicit none

    integer(IK) :: ndim, npnt, nsam, isam

    type(display_type) :: disp
    disp = display_type(file = "main.out.F90")

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Mahalanobis distance squared for real-valued arguments.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        use pm_kind, only: RKC => RKS ! all other real kinds are also supported.
        real(RKC), allocatable :: point(:,:), mahalSq(:,:), mean(:,:), invCov(:,:,:)

        call disp%show("ndim = 3; npnt = 5; nsam = 2")
                        ndim = 3; npnt = 5; nsam = 2
        call disp%show("point = getUnifRand(0., 1., ndim, npnt)")
                        point = getUnifRand(0., 1., ndim, npnt)
        call disp%show("point")
        call disp%show( point )
        call disp%show("call setResized(mean, [ndim, nsam])")
                        call setResized(mean, [ndim, nsam])
        call disp%show("call setResized(invCov, [ndim, ndim, nsam])")
                        call setResized(invCov, [ndim, ndim, nsam])
        call disp%show("do isam = 1, nsam")
                        do isam = 1, nsam
            call disp%show("invCov(:,:,isam) = getCovRand(mold = 1., ndim = ndim)")
                            invCov(:,:,isam) = getCovRand(mold = 1., ndim = ndim)
            call disp%show("mean(:,isam) = getFilled(isam - 1, ndim)")
                            mean(:,isam) = getFilled(isam - 1, ndim)
        call disp%show("end do")
                        end do
        call disp%show("invCov")
        call disp%show( invCov )
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("call setResized(mahalSq, [nsam, npnt])")
                        call setResized(mahalSq, [nsam, npnt])


        call disp%show("mahalSq = getMahalSq(point, invCov) ! for a set of ndim-dimensional points against nsam independent samples centered at origin.")
                        mahalSq = getMahalSq(point, invCov)
        call disp%show("mahalSq")
        call disp%show( mahalSq )
        call disp%show("mahalSq = getMahalSq(point, invCov, mean) ! for a set of ndim-dimensional points against nsam independent samples centered at mean.")
                        mahalSq = getMahalSq(point, invCov, mean)
        call disp%show("mahalSq")
        call disp%show( mahalSq )

        call disp%show("mahalSq(1:nsam, 1) = getMahalSq(point(1:ndim, 1), invCov) ! for a single ndim-dimensional point against nsam independent samples centered at origin.")
                        mahalSq(1:nsam, 1) = getMahalSq(point(1:ndim, 1), invCov)
        call disp%show("mahalSq(1:nsam, 1)")
        call disp%show( mahalSq(1:nsam, 1) )
        call disp%show("mahalSq(1:nsam, 1) = getMahalSq(point(1:ndim, 1), invCov, mean) ! for a single ndim-dimensional point against nsam independent samples centered at mean.")
                        mahalSq(1:nsam, 1) = getMahalSq(point(1:ndim, 1), invCov, mean)
        call disp%show("mahalSq(1:nsam, 1)")
        call disp%show( mahalSq(1:nsam, 1) )


        call disp%show("mahalSq(1, 1:npnt) = getMahalSq(point, invCov(:,:,1)) ! for a set of ndim-dimensional points against one sample centered at origin.")
                        mahalSq(1, 1:npnt) = getMahalSq(point, invCov(:,:,1))
        call disp%show("mahalSq(1, 1:npnt)")
        call disp%show( mahalSq(1, 1:npnt) )
        call disp%show("mahalSq(1, 1:npnt) = getMahalSq(point, invCov(:,:,1), mean(:,1)) ! for a set of ndim-dimensional points against one sample centered at mean.")
                        mahalSq(1, 1:npnt) = getMahalSq(point, invCov(:,:,1), mean(:,1))
        call disp%show("mahalSq(1, 1:npnt)")
        call disp%show( mahalSq(1, 1:npnt) )

        call disp%show("mahalSq(1, 1) = getMahalSq(point(1:ndim, 1), invCov(:,:,1)) ! for a single ndim-dimensional point one sample centered at origin.")
                        mahalSq(1, 1) = getMahalSq(point(1:ndim, 1), invCov(:,:,1))
        call disp%show("mahalSq(1, 1)")
        call disp%show( mahalSq(1, 1) )
        call disp%show("mahalSq(1, 1) = getMahalSq(point(1:ndim, 1), invCov(:,:,1), mean(:,1)) ! for a single ndim-dimensional point one sample centered at mean.")
                        mahalSq(1, 1) = getMahalSq(point(1:ndim, 1), invCov(:,:,1), mean(:,1))
        call disp%show("mahalSq(1, 1)")
        call disp%show( mahalSq(1, 1) )


        call disp%show("mahalSq(1, 1:npnt) = getMahalSq(point(1, 1:npnt), invCov(1,1,1)) ! for a set of single-dimensional points against one sample centered at origin.")
                        mahalSq(1, 1:npnt) = getMahalSq(point(1, 1:npnt), invCov(1,1,1))
        call disp%show("mahalSq(1, 1:npnt)")
        call disp%show( mahalSq(1, 1:npnt) )
        call disp%show("mahalSq(1, 1:npnt) = getMahalSq(point(1, 1:npnt), invCov(1,1,1), mean(1,1)) ! for a set of single-dimensional points against one sample centered at mean.")
                        mahalSq(1, 1:npnt) = getMahalSq(point(1, 1:npnt), invCov(1,1,1), mean(1,1))
        call disp%show("mahalSq(1, 1:npnt)")
        call disp%show( mahalSq(1, 1:npnt) )

        call disp%show("mahalSq(1, 1) = getMahalSq(point(1, 1), invCov(1,1,1)) ! for a single single-dimensional point one sample centered at origin.")
                        mahalSq(1, 1) = getMahalSq(point(1, 1), invCov(1,1,1))
        call disp%show("mahalSq(1, 1)")
        call disp%show( mahalSq(1, 1) )
        call disp%show("mahalSq(1, 1) = getMahalSq(point(1, 1), invCov(1,1,1), mean(1,1)) ! for a single single-dimensional point one sample centered at mean.")
                        mahalSq(1, 1) = getMahalSq(point(1, 1), invCov(1,1,1), mean(1,1))
        call disp%show("mahalSq(1, 1)")
        call disp%show( mahalSq(1, 1) )

    end block

    call disp%skip()
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("! Compute the Mahalanobis distance squared for complex-valued arguments.")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%show("!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    call disp%skip()

    block

        use pm_kind, only: CKC => CKS ! all other complex kinds are also supported.
        complex(CKC), allocatable :: point(:,:), mahalSq(:,:), mean(:,:), invCov(:,:,:)

        call disp%show("ndim = 3; npnt = 5; nsam = 2")
                        ndim = 3; npnt = 5; nsam = 2
        call disp%show("point = getUnifRand(0., 1., ndim, npnt)")
                        point = getUnifRand(0., 1., ndim, npnt)
        call disp%show("point")
        call disp%show( point )
        call disp%show("call setResized(mean, [ndim, nsam])")
                        call setResized(mean, [ndim, nsam])
        call disp%show("call setResized(invCov, [ndim, ndim, nsam])")
                        call setResized(invCov, [ndim, ndim, nsam])
        call disp%show("do isam = 1, nsam")
                        do isam = 1, nsam
            call disp%show("invCov(:,:,isam) = getCovRand(mold = (1., 1.), ndim = ndim)")
                            invCov(:,:,isam) = getCovRand(mold = (1., 1.), ndim = ndim)
            call disp%show("mean(:,isam) = getFilled(isam - 1, ndim)")
                            mean(:,isam) = getFilled(isam - 1, ndim)
        call disp%show("end do")
                        end do
        call disp%show("invCov")
        call disp%show( invCov )
        call disp%show("mean")
        call disp%show( mean )
        call disp%show("call setResized(mahalSq, [nsam, npnt])")
                        call setResized(mahalSq, [nsam, npnt])


        call disp%show("mahalSq = getMahalSq(point, invCov) ! for a set of ndim-dimensional points against nsam independent samples centered at origin.")
                        mahalSq = getMahalSq(point, invCov)
        call disp%show("mahalSq")
        call disp%show( mahalSq )
        call disp%show("mahalSq = getMahalSq(point, invCov, mean) ! for a set of ndim-dimensional points against nsam independent samples centered at mean.")
                        mahalSq = getMahalSq(point, invCov, mean)
        call disp%show("mahalSq")
        call disp%show( mahalSq )

        call disp%show("mahalSq(1:nsam, 1) = getMahalSq(point(1:ndim, 1), invCov) ! for a single ndim-dimensional point against nsam independent samples centered at origin.")
                        mahalSq(1:nsam, 1) = getMahalSq(point(1:ndim, 1), invCov)
        call disp%show("mahalSq(1:nsam, 1)")
        call disp%show( mahalSq(1:nsam, 1) )
        call disp%show("mahalSq(1:nsam, 1) = getMahalSq(point(1:ndim, 1), invCov, mean) ! for a single ndim-dimensional point against nsam independent samples centered at mean.")
                        mahalSq(1:nsam, 1) = getMahalSq(point(1:ndim, 1), invCov, mean)
        call disp%show("mahalSq(1:nsam, 1)")
        call disp%show( mahalSq(1:nsam, 1) )


        call disp%show("mahalSq(1, 1:npnt) = getMahalSq(point, invCov(:,:,1)) ! for a set of ndim-dimensional points against one sample centered at origin.")
                        mahalSq(1, 1:npnt) = getMahalSq(point, invCov(:,:,1))
        call disp%show("mahalSq(1, 1:npnt)")
        call disp%show( mahalSq(1, 1:npnt) )
        call disp%show("mahalSq(1, 1:npnt) = getMahalSq(point, invCov(:,:,1), mean(:,1)) ! for a set of ndim-dimensional points against one sample centered at mean.")
                        mahalSq(1, 1:npnt) = getMahalSq(point, invCov(:,:,1), mean(:,1))
        call disp%show("mahalSq(1, 1:npnt)")
        call disp%show( mahalSq(1, 1:npnt) )

        call disp%show("mahalSq(1, 1) = getMahalSq(point(1:ndim, 1), invCov(:,:,1)) ! for a single ndim-dimensional point one sample centered at origin.")
                        mahalSq(1, 1) = getMahalSq(point(1:ndim, 1), invCov(:,:,1))
        call disp%show("mahalSq(1, 1)")
        call disp%show( mahalSq(1, 1) )
        call disp%show("mahalSq(1, 1) = getMahalSq(point(1:ndim, 1), invCov(:,:,1), mean(:,1)) ! for a single ndim-dimensional point one sample centered at mean.")
                        mahalSq(1, 1) = getMahalSq(point(1:ndim, 1), invCov(:,:,1), mean(:,1))
        call disp%show("mahalSq(1, 1)")
        call disp%show( mahalSq(1, 1) )


        call disp%show("mahalSq(1, 1:npnt) = getMahalSq(point(1, 1:npnt), invCov(1,1,1)) ! for a set of single-dimensional points against one sample centered at origin.")
                        mahalSq(1, 1:npnt) = getMahalSq(point(1, 1:npnt), invCov(1,1,1))
        call disp%show("mahalSq(1, 1:npnt)")
        call disp%show( mahalSq(1, 1:npnt) )
        call disp%show("mahalSq(1, 1:npnt) = getMahalSq(point(1, 1:npnt), invCov(1,1,1), mean(1,1)) ! for a set of single-dimensional points against one sample centered at mean.")
                        mahalSq(1, 1:npnt) = getMahalSq(point(1, 1:npnt), invCov(1,1,1), mean(1,1))
        call disp%show("mahalSq(1, 1:npnt)")
        call disp%show( mahalSq(1, 1:npnt) )

        call disp%show("mahalSq(1, 1) = getMahalSq(point(1, 1), invCov(1,1,1)) ! for a single single-dimensional point one sample centered at origin.")
                        mahalSq(1, 1) = getMahalSq(point(1, 1), invCov(1,1,1))
        call disp%show("mahalSq(1, 1)")
        call disp%show( mahalSq(1, 1) )
        call disp%show("mahalSq(1, 1) = getMahalSq(point(1, 1), invCov(1,1,1), mean(1,1)) ! for a single single-dimensional point one sample centered at mean.")
                        mahalSq(1, 1) = getMahalSq(point(1, 1), invCov(1,1,1), mean(1,1))
        call disp%show("mahalSq(1, 1)")
        call disp%show( mahalSq(1, 1) )

    end block

end program example