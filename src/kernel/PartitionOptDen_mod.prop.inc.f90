
                    ! Compute the cumulative sum of the cluster Size array.

                    ParLev(is)%CumSumSize(0) = ParLev(is)%nps - 1
                    do ic = 1, ParLev(is)%nc ! loopDetermineClusterBoundaries
                        KmeansMemberCounter(ic) = ParLev(is)%CumSumSize(ic-1)
                        ParLev(is)%CumSumSize(ic) = KmeansMemberCounter(ic) + ParLev(is)%Size(ic)
                    end do ! loopDetermineClusterBoundaries

                    do ip = ParLev(is)%nps, ParLev(is)%npe
                        KmeansMemberCounter(ParLev(is)%Membership(ip)) = KmeansMemberCounter(ParLev(is)%Membership(ip)) + 1
                        NormedPoint(1:nd,KmeansMemberCounter(ParLev(is)%Membership(ip))) = Point(1:nd,ip)
                        PointIndex(KmeansMemberCounter(ParLev(is)%Membership(ip))) = Partition%PointIndex(ip)
                    end do
                    Point(1:nd,ParLev(is)%nps:ParLev(is)%npe) = NormedPoint
                    Partition%PointIndex(ParLev(is)%nps:ParLev(is)%npe) = PointIndex

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ! Compute the properties of non-singular clusters.
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    maxLogVolNormed = NEGINF_RK
                    do ic = 1, ParLev(is)%nc ! loopComputeClusterProperties

                        if (ParLev(is)%Size(ic) > nd) then ! blockMinimumClusterSize

                            ipstart = ParLev(is)%CumSumSize(ic-1) + 1
                            ipend = ParLev(is)%CumSumSize(ic)

                            ! Correct the cluster memberships.

                            ParLev(is)%Membership(ipstart:ipend) = ic

                            ! Normalize points.

                            do concurrent(ip = 1:Partition%np) ! ParLev(is)%nps:ParLev(is)%npe)
                                NormedPoint(1:nd,ip) = Point(1:nd,ip) - ParLev(is)%Center(1:nd,ic)
                            end do

                            ! Compute the upper covariance matrix of the cluster covariance matrices.

                            ParLev(is)%ChoLowCovUpp(1:nd,1:nd,ic) = 0._RK
                            do ip = ipstart, ipend
                                do j = 1, nd
                                    ParLev(is)%ChoLowCovUpp(1:j,j,ic) = ParLev(is)%ChoLowCovUpp(1:j,j,ic) + NormedPoint(1:j,ip) * NormedPoint(j,ip)
                                end do
                            end do
                            ParLev(is)%ChoLowCovUpp(1:nd,1:nd,ic) = ParLev(is)%ChoLowCovUpp(1:nd,1:nd,ic) / real(ParLev(is)%Size(ic)-1,RK)

                            ! Compute the Cholesky Factor of the cluster covariance matrices.

                            call getCholeskyFactor(nd, ParLev(is)%ChoLowCovUpp(1:nd,1:nd,ic), ParLev(is)%ChoDia(1:nd,ic))
                            if (ParLev(is)%ChoDia(1,ic) < 0._RK) then
                                ! LCOV_EXCL_START
                                !ParLev(is)%Err%msg = PROCEDURE_NAME//"Cholesky factorization failed."
                                !ParLev(is)%Err%occurred = .true.
                                error stop ! xxx this must be fixed.
                                Partition%Size(ic) = ParLev(is)%np
                                neopt = 1
                                return
                                ! LCOV_EXCL_STOP
                            end if

                            ! Compute the inverse of the cluster covariance matrices.

                            ParLev(is)%InvCovMat(1:nd,1:nd,ic) = getInvMatFromCholFac( nd, ParLev(is)%ChoLowCovUpp(1:nd,1:nd,ic), ParLev(is)%ChoDia(1:nd,ic) )

                            ! Compute the MahalSq of as many points as needed and the scaleFcator of the bounding region.

                            do concurrent(ip = 1:ipstart-1)
                                ParLev(is)%MahalSq(ip,ic) = dot_product( NormedPoint(1:nd) , matmul(ParLev(is)%InvCovMat(1:nd,1:nd,ic), NormedPoint(1:nd)) )
                            end do

                            ParLev(is)%ScaleFactorSq(ic) = NEGINF_RK
                            do ip = ipstart, ipend
                                ParLev(is)%MahalSq(ip,ic) = dot_product( NormedPoint(1:nd,ip) , matmul(ParLev(is)%InvCovMat(1:nd,1:nd,ic), NormedPoint(1:nd,ip)) )
                                if (ParLev(is)%ScaleFactorSq(ic) < ParLev(is)%MahalSq(ip,ic)) ParLev(is)%ScaleFactorSq(ic) = ParLev(is)%MahalSq(ip,ic)
                            end do

                            do concurrent(ip = ipend+1:Partition%np)
                                ParLev(is)%MahalSq(ip,ic) = dot_product( NormedPoint(1:nd) , matmul(ParLev(is)%InvCovMat(1:nd,1:nd,ic), NormedPoint(1:nd)) )
                            end do

                            ! Scale the properties to form the bounded regions.

                            scaleFactorSqInverse = 1._RK / ParLev(is)%ScaleFactorSq(ic)
                            scaleFactor = sqrt(ParLev(is)%ScaleFactorSq(ic))
                            ParLev(is)%ChoDia(1:nd,ic) = ParLev(is)%ChoDia(1:nd,ic) * scaleFactor
                            ParLev(is)%InvCovMat(1:nd,1:nd,ic) = ParLev(is)%InvCovMat(1:nd,1:nd,ic) * scaleFactorSqInverse
                            do j = 1, nd
                                do i = j + 1, nd
                                    ParLev(is)%ChoLowCovUpp(i,j,ic) = ParLev(is)%ChoLowCovUpp(i,j,ic) * scaleFactor
                                end do
                            end do

                            ! Rescale the MahalSq to the bounding regions

                            !ParLev(is)%MahalSq(1:ParLev(is)%np,ic) = ParLev(is)%MahalSq(1:ParLev(is)%np,ic) * scaleFactorSqInverse

                            ! Compute the effective fraction of points inside the bounded region.

                            if (Partition%inclusionFraction > 0._RK) then
                                ParLev(is)%EffectiveSize(ic)    = ParLev(is)%Size(ic) & ! LCOV_EXCL_LINE
                                                                + nint(Partition%inclusionFraction * & ! LCOV_EXCL_LINE
                                                                ( count(ParLev(is)%MahalSq(1:ParLev(is)%CumSumSize(ic-1),ic)<=ParLev(is)%ScaleFactorSq(ic)) & ! LCOV_EXCL_LINE
                                                                + count(ParLev(is)%MahalSq(ParLev(is)%CumSumSize(ic)+1:Partition%np,ic)<=ParLev(is)%ScaleFactorSq(ic)) & ! LCOV_EXCL_LINE
                                                                ) )
                            else
                                ParLev(is)%EffectiveSize(ic) = ParLev(is)%Size(ic)
                            end if

                            ParLev(is)%LogVolNormed(ic) = sum( log(ParLev(is)%ChoDia(1:nd,ic)) )
                            maxLogVolNormed = max( maxLogVolNormed, ParLev(is)%LogVolNormed(ic) )

                        else ! blockMinimumClusterSize

                            ParLev(is)%LogVolNormed(ic) = NEGINF_RK ! This initialization is essential in dependent procedures.
                            ParLev(is)%EffectiveSize(ic) = 0

                        end if ! blockMinimumClusterSize

                    end do ! loopComputeClusterProperties

                    if (maxLogVolNormed > NEGINF_RK) then ! blockSingularClusterDetected
                    
                        ParLev(is)%logSumVolNormed = getLogSumExp(ParLev(is)%nc, ParLev(is)%LogVolNormed, maxLogVolNormed)

                    else ! blockSingularClusterDetected

                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        ! Compute the properties of singular clusters.
                        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        ParLev(is)%logSumVolNormed = 0._RK ! exp(ParLev(is)%logSumVolNormed - maxLogVolNormed)

                        ! Assume a hyper-spherical ellipsoid as the bounding region of the points.

                        do ic = 1, ParLev(is)%nc ! loopComputeSingularClusterProperties

                            if (ParLev(is)%Size(ic) <= nd .and. ParLev(is)%Size(ic) > 0) then

                                ipstart = ParLev(is)%CumSumSize(ic-1) + 1
                                ipend = ParLev(is)%CumSumSize(ic)

                                ! Correct the cluster memberships.

                                ParLev(is)%Membership(ipstart:ipend) = ic

                                ! Compute the scale factor and other properties.

                                do concurrent(ip = ipstart:ipend) ! 1:Partition%np)
                                    ParLev(is)%MahalSq(ip,ic) = sum( ( Point(1:nd,ip) - ParLev(is)%Center(1:nd,ic) )**2 )
                                end do

                                ParLev(is)%ScaleFactorSq(ic) = max(exp(2*Partition%pointLogVolNormed/nd), maxval(ParLev(is)%MahalSq(ipstart:ipend,ic)))
                                scaleFactorSqInverse = 1._RK / ParLev(is)%ScaleFactorSq(ic)
                                scaleFactor = sqrt(ParLev(is)%ScaleFactorSq(ic))

                                !ParLev(is)%MahalSq(ip,ic) = ParLev(is)%MahalSq(ip,ic) * scaleFactorSqInverse

                                ParLev(is)%ChoDia(1:nd,ic) = scaleFactor
                                ParLev(is)%LogVolNormed(ic) = nd * log(scaleFactor)
                                ParLev(is)%InvCovMat(1:nd,1:nd,ic) = getEye(nd, nd, scaleFactorSqInverse)
                                ParLev(is)%ChoLowCovUpp(1:nd,1:nd,ic) = getEye(nd, nd, ParLev(is)%ScaleFactorSq(ic))

                                ! Compute the effective fraction of points inside the bounded region.

                                !if (Partition%inclusionFraction > 0._RK) then
                                !    ParLev(is)%EffectiveSize(ic)   = ParLev(is)%Size(ic) & ! LCOV_EXCL_LINE
                                !                                    + nint(Partition%inclusionFraction * & ! LCOV_EXCL_LINE
                                !                                    ( count(ParLev(is)%MahalSq(1:ParLev(is)%CumSumSize(ic-1),ic)<=ParLev(is)%ScaleFactorSq(ic)) & ! LCOV_EXCL_LINE
                                !                                    + count(ParLev(is)%MahalSq(ParLev(is)%CumSumSize(ic)+1:Partition%np,ic)<=ParLev(is)%ScaleFactorSq(ic)) & ! LCOV_EXCL_LINE
                                !                                    ), kind = IK)
                                !else
                                !    ParLev(is)%EffectiveSize(ic) = ParLev(is)%Size(ic)
                                !end if

                                ParLev(is)%logSumVolNormed = ParLev(is)%logSumVolNormed + exp(ParLev(is)%LogVolNormed(ic) - maxLogVolNormed)

                            end if

                        end do ! loopComputeSingularClusterProperties

                        ParLev(is)%logSumVolNormed = maxLogVolNormed + log(ParLev(is)%logSumVolNormed)

                    end if ! blockSingularClusterDetected

#if (defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED)
                    block
                        use Statistics_mod, only: getSamCholFac, getMahalSq
                        real(RK), parameter :: TOLERANCE = 1.e-2_RK
                        real(RK) :: mahalSqScalar
                        real(RK) :: NormedPoint(nd)
                        real(RK) :: dum, diff, scaleFac
                        real(RK) :: ChoLow(nd,nd), ChoDia(nd)
                        real(RK) :: ChoLowDiff(nd,nd), ChoDiaDiff(nd)
                        loopComputeScaleFactor: do ic = 1, ParLev(is)%nc
                            call getSamCholFac  ( nd = nd & ! LCOV_EXCL_LINE
                                                , np = ParLev(is)%Size(ic) & ! LCOV_EXCL_LINE
                                                , Mean = ParLev(is)%Center(1:nd,ic) & ! LCOV_EXCL_LINE
                                                , Point = Point(1:nd,ParLev(is)%CumSumSize(ic-1)+1:ParLev(is)%CumSumSize(ic)) & ! LCOV_EXCL_LINE
                                                , CholeskyLower = ChoLow & ! LCOV_EXCL_LINE
                                                , CholeskyDiago = ChoDia & ! LCOV_EXCL_LINE
                                                )
                            scaleFac = sqrt(ParLev(is)%ScaleFactorSq(ic))
                            ChoDia = ChoDia * scaleFac
                            do j = 1, nd
                                do i = j+1, nd
                                    ChoLow(i,j) = ChoLow(i,j) * scaleFac
                                end do
                            end do
                            ChoLowDiff = 2 * abs(ParLev(is)%ChoLowCovUpp(1:nd,1:nd,ic)-ChoLow) / abs(ParLev(is)%ChoLowCovUpp(1:nd,1:nd,ic)+ChoLow)
                            if ( any(ChoLowDiff > TOLERANCE) ) then
                                ! LCOV_EXCL_START
                                write(*,*)
                                write(*,*) PROCEDURE_NAME
                                write(*,*) "ParLev(is)%ChoLowCovUpp(1:nd,1:nd,ic)-ChoLow) > TOLERANCE"
                                write(*,*) "TOLERANCE"
                                write(*,*) TOLERANCE
                                write(*,*) "ChoLowDiff"
                                write(*,*) ChoLowDiff
                                write(*,*) "ChoLow"
                                write(*,*) ChoLow
                                write(*,*) "Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic)"
                                write(*,*) ParLev(is)%ChoLowCovUpp(1:nd,1:nd,ic)
                                write(*,*)
                                error stop
                                ! LCOV_EXCL_STOP
                            end if
                            ChoDiaDiff = 2 * abs(ParLev(is)%ChoDia(1:nd,ic)-ChoDia) / abs(ParLev(is)%ChoDia(1:nd,ic)+ChoDia)
                            if ( any(ChoDiaDiff > TOLERANCE) ) then
                                ! LCOV_EXCL_START
                                write(*,*)
                                write(*,*) "ParLev(is)%ChoDia(1:nd,ic)-ChoDia) > TOLERANCE"
                                write(*,*) "TOLERANCE"
                                write(*,*) TOLERANCE
                                write(*,*) "ChoDiaDiff"
                                write(*,*) ChoDiaDiff
                                write(*,*) "ChoDia"
                                write(*,*) ChoDia
                                write(*,*) "ParLev(is)%ChoDia(1:nd,1:nd,ic)"
                                write(*,*) ParLev(is)%ChoDia(1:nd,ic)
                                write(*,*)
                                !error stop
                                ! LCOV_EXCL_STOP
                            end if
                            do ip = ParLev(is)%nps, ParLev(is)%npe
                                NormedPoint(1:nd) = Point(1:nd,ip) - ParLev(is)%Center(1:nd,ic)
                                mahalSqScalar = dot_product( NormedPoint , matmul(ParLev(is)%InvCovMat(1:nd,1:nd,ic), NormedPoint) )
                                if (mahalSqScalar<0._RK) then
                                ! LCOV_EXCL_START
                                    mahalSqScalar = -1._RK
                                    write(*,*) PROCEDURE_NAME
                                    write(*,*) "mahalSqScalar<0._RK", ip, mahalSqScalar
                                    error stop
                                end if
                                ! LCOV_EXCL_STOP
                                dum = getMahalSq( nd, ParLev(is)%Center(1:nd,ic), ParLev(is)%InvCovMat(1:nd,1:nd,ic), Point(1:nd,ip) )
                                diff = 2 * abs(dum-mahalSqScalar) / abs(dum+mahalSqScalar)
                                if (diff > TOLERANCE) then
                                    ! LCOV_EXCL_START
                                    write(*,*)
                                    write(*,*) "dum/=mahalSqScalar", ip, mahalSqScalar, dum
                                    write(*,*) "TOLERANCE"
                                    write(*,*) TOLERANCE
                                    write(*,*) "diff"
                                    write(*,*) diff
                                    write(*,*)
                                    !error stop
                                    ! LCOV_EXCL_STOP
                                end if
                            end do
                            !MahalSqScalar(1:np,ic) = MahalSqScalar(1:np,ic) / Kmeans%Prop%ScaleFactorSq(ic)
                        end do loopComputeScaleFactor
                    end block
#endif