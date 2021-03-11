
        ! Compute the cumulative sum of the cluster Size array.

        Kmeans%Prop%CumSumSize(0) = nps - 1
        do ic = 1, nc ! loopDetermineClusterBoundaries
            KmeansMemberCounter(ic) = Kmeans%Prop%CumSumSize(ic-1)
            Kmeans%Prop%CumSumSize(ic) = KmeansMemberCounter(ic) + Kmeans%Size(ic)
        end do ! loopDetermineClusterBoundaries

        do ip = nps, npe
            jp = ip - nps + 1
            KmeansMemberCounter(Kmeans%Membership(jp)) = KmeansMemberCounter(Kmeans%Membership(jp)) + 1_IK
            StanPoint(1:nd,KmeansMemberCounter(Kmeans%Membership(jp))) = Point(1:nd,ip)
            Kmeans%Prop%Index(KmeansMemberCounter(Kmeans%Membership(jp))) = PointIndex(ip)
        end do
        PointIndex(nps:npe) = Kmeans%Prop%Index
        Point(1:nd,nps:npe) = StanPoint

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the properties of non-singular clusters.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        maxLogVolNormed = NEGINF_RK
        pointLogVolNormedDefault = -1._RK ! indicator: no singular cluster detected.
        Kmeans%Prop%logSumVolNormed = NEGINF_RK
        do ic = 1, nc ! loopComputeClusterProperties

            if (Kmeans%Size(ic) > nd) then ! blockMinimumClusterSize

                ipstart = Kmeans%Prop%CumSumSize(ic-1) + 1_IK
                ipend = Kmeans%Prop%CumSumSize(ic)

                ! Correct the cluster memberships.

                Kmeans%Membership(ipstart-nps+1:ipend-nps+1) = ic

                ! Normalize points.

                do concurrent(ip = nps:npe)
                    StanPoint(1:nd,ip) = Point(1:nd,ip) - Kmeans%Center(1:nd,ic)
                end do

                ! Compute the upper covariance matrix of the cluster covariance matrices.

                Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic) = 0._RK
                do ip = ipstart, ipend
                    do j = 1, nd
                        Kmeans%Prop%ChoLowCovUpp(1:j,j,ic) = Kmeans%Prop%ChoLowCovUpp(1:j,j,ic) + StanPoint(1:j,ip) * StanPoint(j,ip)
                    end do
                end do
                Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic) = Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic) / real(Kmeans%Size(ic)-1_IK, kind = RK)

                ! Compute the Cholesky Factor of the cluster covariance matrices.

                call getCholeskyFactor(nd, Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic), Kmeans%Prop%ChoDia(1:nd,ic))
                if (Kmeans%Prop%ChoDia(1,ic) < 0._RK) then
                    ! LCOV_EXCL_START
                    !Kmeans%Err%msg = PROCEDURE_NAME//"Cholesky factorization failed."
                    !Kmeans%Err%occurred = .true.
                    PartitionSize(1) = npc
                    neopt = 1
                    return
                    ! LCOV_EXCL_STOP
                end if

                ! Compute the inverse of the cluster covariance matrices.

                Kmeans%Prop%InvCovMat(1:nd,1:nd,ic) = getInvMatFromCholFac( nd, Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic), Kmeans%Prop%ChoDia(1:nd,ic) )

                ! Compute the MahalSq of as many points as needed.

                do concurrent(ip = 1:nps-1)
                    NormedPoint(1:nd) = Point(1:nd,ip) - Kmeans%Center(1:nd,ic)
                    Kmeans%Prop%MahalSq(ip,ic) = dot_product( NormedPoint(1:nd) , matmul(Kmeans%Prop%InvCovMat(1:nd,1:nd,ic), NormedPoint(1:nd)) )
                end do

                do concurrent(ip = nps:npe)
                    Kmeans%Prop%MahalSq(ip,ic) = dot_product( StanPoint(1:nd,ip) , matmul(Kmeans%Prop%InvCovMat(1:nd,1:nd,ic), StanPoint(1:nd,ip)) )
                end do

                do concurrent(ip = npe+1:np)
                    NormedPoint(1:nd) = Point(1:nd,ip) - Kmeans%Center(1:nd,ic)
                    Kmeans%Prop%MahalSq(ip,ic) = dot_product( NormedPoint(1:nd) , matmul(Kmeans%Prop%InvCovMat(1:nd,1:nd,ic), NormedPoint(1:nd)) )
                end do

                ! Compute the scaleFcator of the bounding region and scale the properties to form the bounded regions.

                Kmeans%Prop%ScaleFactorSq(ic) = maxval(Kmeans%Prop%MahalSq(ipstart:ipend,ic))

                scaleFactorSqInverse = 1._RK / Kmeans%Prop%ScaleFactorSq(ic)
                scaleFactor = sqrt(Kmeans%Prop%ScaleFactorSq(ic))
                Kmeans%Prop%ChoDia(1:nd,ic) = Kmeans%Prop%ChoDia(1:nd,ic) * scaleFactor
                Kmeans%Prop%InvCovMat(1:nd,1:nd,ic) = Kmeans%Prop%InvCovMat(1:nd,1:nd,ic) * scaleFactorSqInverse
                do j = 1, nd
                    do i = j + 1, nd
                        Kmeans%Prop%ChoLowCovUpp(i,j,ic) = Kmeans%Prop%ChoLowCovUpp(i,j,ic) * scaleFactor
                    end do
                end do

                ! Rescale the MahalSq to the bounding regions

                !Kmeans%Prop%MahalSq(1:np,ic) = Kmeans%Prop%MahalSq(1:np,ic) * scaleFactorSqInverse

                ! Compute the effective fraction of points inside the bounded region.

                if (inclusionFraction > 0._RK) then
                    Kmeans%Prop%EffectiveSize(ic)   = Kmeans%Size(ic) & ! LCOV_EXCL_LINE
                                                    + nint(inclusionFraction * & ! LCOV_EXCL_LINE
                                                    ( count(Kmeans%Prop%MahalSq(1:Kmeans%Prop%CumSumSize(ic-1),ic)<=Kmeans%Prop%ScaleFactorSq(ic)) & ! LCOV_EXCL_LINE
                                                    + count(Kmeans%Prop%MahalSq(Kmeans%Prop%CumSumSize(ic)+1:np,ic)<=Kmeans%Prop%ScaleFactorSq(ic)) & ! LCOV_EXCL_LINE
                                                    ) )
                else
                    Kmeans%Prop%EffectiveSize(ic) = Kmeans%Size(ic)
                end if

                Kmeans%Prop%LogVolNormed(ic) = sum( log(Kmeans%Prop%ChoDia(1:nd,ic)) )
                maxLogVolNormed = max( maxLogVolNormed, Kmeans%Prop%LogVolNormed(ic) )

            else ! blockMinimumClusterSize

                pointLogVolNormedDefault = 1._RK ! This is an indicator
                Kmeans%Prop%LogVolNormed(ic) = NEGINF_RK ! This initialization is essential in dependent procedures.
                Kmeans%Prop%EffectiveSize(ic) = 0_IK

            end if ! blockMinimumClusterSize

        end do ! loopComputeClusterProperties

        if (maxLogVolNormed /= NEGINF_RK) Kmeans%Prop%logSumVolNormed = getLogSumExp(nc, Kmeans%Prop%LogVolNormed, maxLogVolNormed)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the properties of singular clusters.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (pointLogVolNormedDefault > 0._RK) then ! blockSingularClusterDetected

            if (pointLogVolNormed > NEGINF_RK) then
                pointLogVolNormedDefault = pointLogVolNormed
            elseif (maxLogVolNormed /= NEGINF_RK) then
                pointLogVolNormedDefault = 0._RK
                do ic = 1, nc ! add all points.
                    if (Kmeans%Size(ic) > nd) pointLogVolNormedDefault = pointLogVolNormedDefault + real(Kmeans%Size(ic),RK)
                end do
                pointLogVolNormedDefault = Kmeans%Prop%logSumVolNormed - log(pointLogVolNormedDefault)
            else ! There is not any non-singular cluster for which to compute the properties.
                ! LCOV_EXCL_START
                PartitionSize(1) = npc
                neopt = 1
                return
                ! LCOV_EXCL_STOP
            end if

            if (maxLogVolNormed /= NEGINF_RK) then
                Kmeans%Prop%logSumVolNormed = exp(Kmeans%Prop%logSumVolNormed - maxLogVolNormed)
            else
                Kmeans%Prop%logSumVolNormed = 0._RK
                maxLogVolNormed = pointLogVolNormedDefault
            end if

            ! Assume a hyper-spherical ellipsoid as the bounding region of the points.

            do ic = 1, nc ! loopComputeSingularClusterProperties

                if (Kmeans%Size(ic) <= nd .and. Kmeans%Size(ic) > 0_IK) then

                    ipstart = Kmeans%Prop%CumSumSize(ic-1) + 1_IK
                    ipend = Kmeans%Prop%CumSumSize(ic)

                    ! Correct the cluster memberships.

                    Kmeans%Membership(ipstart-nps+1:ipend-nps+1) = ic

                    ! Compute the scale factor and other properties.

                    do concurrent(ip = 1:np)
                        Kmeans%Prop%MahalSq(ip,ic) = sum( ( Point(1:nd,ip) - Kmeans%Center(1:nd,ic) )**2 )
                    end do

                    Kmeans%Prop%ScaleFactorSq(ic) = max(exp(2*pointLogVolNormedDefault/nd), maxval(Kmeans%Prop%MahalSq(ipstart:ipend,ic)))
                    scaleFactorSqInverse = 1._RK / Kmeans%Prop%ScaleFactorSq(ic)
                    scaleFactor = sqrt(Kmeans%Prop%ScaleFactorSq(ic))

                    !Kmeans%Prop%MahalSq(ip,ic) = Kmeans%Prop%MahalSq(ip,ic) * scaleFactorSqInverse

                    Kmeans%Prop%ChoDia(1:nd,ic) = scaleFactor
                    Kmeans%Prop%LogVolNormed(ic) = nd * log(scaleFactor)
                    Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic) = getEye(nd, nd, Kmeans%Prop%ScaleFactorSq(ic))
                    Kmeans%Prop%InvCovMat(1:nd,1:nd,ic) = getEye(nd, nd, scaleFactorSqInverse)

                    ! Compute the effective fraction of points inside the bounded region.

                    if (inclusionFraction > 0._RK) then
                        Kmeans%Prop%EffectiveSize(ic)   = Kmeans%Size(ic) & ! LCOV_EXCL_LINE
                                                        + nint(inclusionFraction * & ! LCOV_EXCL_LINE
                                                        ( count(Kmeans%Prop%MahalSq(1:Kmeans%Prop%CumSumSize(ic-1),ic)<=Kmeans%Prop%ScaleFactorSq(ic)) & ! LCOV_EXCL_LINE
                                                        + count(Kmeans%Prop%MahalSq(Kmeans%Prop%CumSumSize(ic)+1:np,ic)<=Kmeans%Prop%ScaleFactorSq(ic)) & ! LCOV_EXCL_LINE
                                                        ), kind = IK)
                    else
                        Kmeans%Prop%EffectiveSize(ic) = Kmeans%Size(ic)
                    end if

                    Kmeans%Prop%logSumVolNormed = Kmeans%Prop%logSumVolNormed + exp(Kmeans%Prop%LogVolNormed(ic) - maxLogVolNormed)

                end if

            end do ! loopComputeSingularClusterProperties

            Kmeans%Prop%logSumVolNormed = maxLogVolNormed + log(Kmeans%Prop%logSumVolNormed)

        end if ! blockSingularClusterDetected
