!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!   MIT License
!!!!
!!!!   ParaMonte: plain powerful parallel Monte Carlo library.
!!!!
!!!!   Copyright (C) 2012-present, The Computational Data Science Lab
!!!!
!!!!   This file is part of the ParaMonte library.
!!!!
!!!!   Permission is hereby granted, free of charge, to any person obtaining a
!!!!   copy of this software and associated documentation files (the "Software"),
!!!!   to deal in the Software without restriction, including without limitation
!!!!   the rights to use, copy, modify, merge, publish, distribute, sublicense,
!!!!   and/or sell copies of the Software, and to permit persons to whom the
!!!!   Software is furnished to do so, subject to the following conditions:
!!!!
!!!!   The above copyright notice and this permission notice shall be
!!!!   included in all copies or substantial portions of the Software.
!!!!
!!!!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!!!!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!!!!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
!!!!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
!!!!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
!!!!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
!!!!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!!!!
!!!!   ACKNOWLEDGMENT
!!!!
!!!!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!!!!   As per the ParaMonte library license agreement terms, if you use any parts of
!!!!   this library for any purposes, kindly acknowledge the use of ParaMonte in your
!!!!   work (education/research/industry/development/...) by citing the ParaMonte
!!!!   library as described on this page:
!!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \brief
!> This file is part of PartitionOptDen_mod.f90 source file.
!> \author
!> Amir Shahmoradi
!> Wednesday 3:37 am, March 10, 2021, Dallas, TX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Compute the cumulative sum of the cluster Size array.
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            Partition%Try(it)%CumSumSize(0) = Partition%Try(it)%nps - 1
            if (Partition%Try(it)%nc > 1) then
                do ic = 1, Partition%Try(it)%nc ! loopDetermineClusterBoundaries
                    KmeansMemberCounter(ic) = Partition%Try(it)%CumSumSize(ic-1)
                    Partition%Try(it)%CumSumSize(ic) = KmeansMemberCounter(ic) + Partition%Try(it)%Size(ic)
                end do ! loopDetermineClusterBoundaries
                do ip = Partition%Try(it)%nps, Partition%Try(it)%npe
                    KmeansMemberCounter(Membership(ip)) = KmeansMemberCounter(Membership(ip)) + 1
                    NormedPoint(1:nd,KmeansMemberCounter(Membership(ip))) = Point(1:nd,ip)
                    PointIndex(KmeansMemberCounter(Membership(ip))) = Partition%PointIndex(ip)
                end do
                Point(1:nd,Partition%Try(it)%nps:Partition%Try(it)%npe) = NormedPoint(1:nd,Partition%Try(it)%nps:Partition%Try(it)%npe)
                Partition%PointIndex(Partition%Try(it)%nps:Partition%Try(it)%npe) = PointIndex(Partition%Try(it)%nps:Partition%Try(it)%npe)
            else
                Partition%Try(it)%CumSumSize(1) = Partition%Try(it)%Size(1)
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Compute the properties of partitions.
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            maxLogVolNormed = NEGINF_RK
            do ic = 1, Partition%Try(it)%nc ! loopComputeClusterProperties

                ipstart = Partition%Try(it)%CumSumSize(ic-1) + 1
                ipend = Partition%Try(it)%CumSumSize(ic)

                ! Correct the cluster memberships.

                if (kvolumeEnabled) Membership(ipstart:ipend) = ic ! needed only when kvolume is activated

                ! Normalize points.

                do concurrent(ip = 1:Partition%np) ! Partition%Try(it)%nps:Partition%Try(it)%npe)
                    NormedPoint(1:nd,ip) = Point(1:nd,ip) - Partition%Try(it)%Center(1:nd,ic)
                end do

                if (Partition%Try(it)%Size(ic) > nd) then ! blockMinimumClusterSize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ! Compute the properties of non-singular clusters.
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    ! Compute the upper covariance matrix of the cluster covariance matrices.

                    Partition%Try(it)%ChoLowCovUpp(1:nd,1:nd,ic) = 0._RK
                    do ip = ipstart, ipend
                        do j = 1, nd
                            Partition%Try(it)%ChoLowCovUpp(1:j,j,ic) = Partition%Try(it)%ChoLowCovUpp(1:j,j,ic) + NormedPoint(1:j,ip) * NormedPoint(j,ip)
                        end do
                    end do
                    Partition%Try(it)%ChoLowCovUpp(1:nd,1:nd,ic) = Partition%Try(it)%ChoLowCovUpp(1:nd,1:nd,ic) / real(Partition%Try(it)%Size(ic)-1,RK)

                    ! Compute the Cholesky Factor of the cluster covariance matrices.

                    call getCholeskyFactor(nd, Partition%Try(it)%ChoLowCovUpp(1:nd,1:nd,ic), Partition%Try(it)%ChoDia(1:nd,ic))
                    if (Partition%Try(it)%ChoDia(1,ic) < 0._RK) then
                        ! LCOV_EXCL_START
                        !Partition%Try(it)%Err%msg = PROCEDURE_NAME//"Cholesky factorization failed."
                        !Partition%Try(it)%Err%occurred = .true.
                        error stop ! xxx this must be fixed.
                        Partition%Size(ic) = Partition%Try(it)%np
                        neopt = 1
                        return
                        ! LCOV_EXCL_STOP
                    end if

                    ! Compute the inverse of the cluster covariance matrices.

                    Partition%Try(it)%InvCovMat(1:nd,1:nd,ic) = getInvMatFromCholFac( nd, Partition%Try(it)%ChoLowCovUpp(1:nd,1:nd,ic), Partition%Try(it)%ChoDia(1:nd,ic) )

                    ! Compute the MahalSq of as many points as needed and the scaleFcator of the bounding region.

                    do concurrent(ip = 1:ipstart-1)
                        Partition%Try(it)%MahalSq(ip,ic) = dot_product( NormedPoint(1:nd) , matmul(Partition%Try(it)%InvCovMat(1:nd,1:nd,ic), NormedPoint(1:nd)) )
                    end do

                    Partition%Try(it)%ScaleFactorSq(ic) = NEGINF_RK
                    do ip = ipstart, ipend
                        Partition%Try(it)%MahalSq(ip,ic) = dot_product( NormedPoint(1:nd,ip) , matmul(Partition%Try(it)%InvCovMat(1:nd,1:nd,ic), NormedPoint(1:nd,ip)) )
                        if (Partition%Try(it)%ScaleFactorSq(ic) < Partition%Try(it)%MahalSq(ip,ic)) Partition%Try(it)%ScaleFactorSq(ic) = Partition%Try(it)%MahalSq(ip,ic)
                    end do

                    do concurrent(ip = ipend+1:Partition%np)
                        Partition%Try(it)%MahalSq(ip,ic) = dot_product( NormedPoint(1:nd) , matmul(Partition%Try(it)%InvCovMat(1:nd,1:nd,ic), NormedPoint(1:nd)) )
                    end do

                    ! Scale the properties to form the bounded regions.

                    scaleFactor = sqrt(Partition%Try(it)%ScaleFactorSq(ic))
                    scaleFactorSqInverse = 1._RK / Partition%Try(it)%ScaleFactorSq(ic)
                    !Partition%Try(it)%MahalSq(1:Partition%Try(it)%np,ic) = Partition%Try(it)%MahalSq(1:Partition%Try(it)%np,ic) * scaleFactorSqInverse

                    Partition%Try(it)%ChoDia(1:nd,ic) = Partition%Try(it)%ChoDia(1:nd,ic) * scaleFactor
                    Partition%Try(it)%InvCovMat(1:nd,1:nd,ic) = Partition%Try(it)%InvCovMat(1:nd,1:nd,ic) * scaleFactorSqInverse
                    do j = 1, nd
                        Partition%Try(it)%ChoLowCovUpp(j+1:nd,j,ic) = Partition%Try(it)%ChoLowCovUpp(j+1:nd,j,ic) * scaleFactor
                    end do

                    Partition%Try(it)%LogVolNormed(ic) = sum( log(Partition%Try(it)%ChoDia(1:nd,ic)) )
                    maxLogVolNormed = max(Partition%Try(it)%LogVolNormed(ic), maxLogVolNormed)

                else ! blockMinimumClusterSize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if (defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED)
                if (Partition%Try(it)%Size(ic) < 0) then
                        write(*,*) "Partition%Try(it)%Size(ic) < 0, ic, Partition%Try(it)%Size(ic):", ic, Partition%Try(it)%Size(ic)
                        error stop
                end if
#endif
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ! Compute the properties of singular clusters. Assume a hyper-spherical ellipsoid as the bounding region.
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    ! Compute the scale factor and other properties.

                    do concurrent(ip = 1:ipstart-1)
                        Partition%Try(it)%MahalSq(ip,ic) = sum( ( Point(1:nd,ip) - Partition%Try(it)%Center(1:nd,ic) )**2 )
                    end do

                    Partition%Try(it)%ScaleFactorSq(ic) = NEGINF_RK
                    do ip = ipstart, ipend ! 1:Partition%np
                        Partition%Try(it)%MahalSq(ip,ic) = sum( ( Point(1:nd,ip) - Partition%Try(it)%Center(1:nd,ic) )**2 )
                        if (Partition%Try(it)%ScaleFactorSq(ic) < Partition%Try(it)%MahalSq(ip,ic)) Partition%Try(it)%ScaleFactorSq(ic) = Partition%Try(it)%MahalSq(ip,ic)
                    end do

                    do concurrent(ip = ipend+1:Partition%np)
                        Partition%Try(it)%MahalSq(ip,ic) = sum( ( Point(1:nd,ip) - Partition%Try(it)%Center(1:nd,ic) )**2 )
                    end do

                    !Partition%Try(it)%ScaleFactorSq(ic) = max(exp(Partition%Try(it)%Size(ic)*ndHalfInverse*Partition%pointLogVolNormed), maxval(Partition%Try(it)%MahalSq(ipstart:ipend,ic)))
                    scaleFactorSqInverse = 1._RK / Partition%Try(it)%ScaleFactorSq(ic)
                    scaleFactor = sqrt(Partition%Try(it)%ScaleFactorSq(ic))
                    !Partition%Try(it)%MahalSq(ip,ic) = Partition%Try(it)%MahalSq(ip,ic) * scaleFactorSqInverse

                    Partition%Try(it)%ChoDia(1:nd,ic) = scaleFactor
                    Partition%Try(it)%InvCovMat(1:nd,1:nd,ic) = getEye(nd, nd, scaleFactorSqInverse)
                    Partition%Try(it)%ChoLowCovUpp(1:nd,1:nd,ic) = getEye(nd, nd, Partition%Try(it)%ScaleFactorSq(ic))

                    Partition%Try(it)%LogVolNormed(ic) = nd * log(scaleFactor)
                    maxLogVolNormed = max(Partition%Try(it)%LogVolNormed(ic), maxLogVolNormed)

                end if ! blockMinimumClusterSize %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end do ! loopComputeClusterProperties

            !if (Partition%Try(it)%nc > 1) then
            !    Partition%Try(it)%logSumVolNormed = getLogSumExp(Partition%Try(it)%nc, Partition%Try(it)%LogVolNormed, maxLogVolNormed)
            !else
            !    Partition%Try(it)%logSumVolNormed = Partition%Try(it)%LogVolNormed(1)
            !end if

#if (defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED)
            block
                use Statistics_mod, only: getSamCholFac, getMahalSq
                real(RK), parameter :: TOLERANCE = 1.e-2_RK
                real(RK) :: mahalSqScalar
                real(RK) :: NormedPoint(nd)
                real(RK) :: dum, diff, scaleFac
                real(RK) :: ChoLow(nd,nd), ChoDia(nd)
                real(RK) :: ChoLowDiff(nd,nd), ChoDiaDiff(nd)
                loopComputeScaleFactor: do ic = 1, Partition%Try(it)%nc
                    call getSamCholFac  ( nd = nd & ! LCOV_EXCL_LINE
                                        , np = Partition%Try(it)%Size(ic) & ! LCOV_EXCL_LINE
                                        , Mean = Partition%Try(it)%Center(1:nd,ic) & ! LCOV_EXCL_LINE
                                        , Point = Point(1:nd,Partition%Try(it)%CumSumSize(ic-1)+1:Partition%Try(it)%CumSumSize(ic)) & ! LCOV_EXCL_LINE
                                        , CholeskyLower = ChoLow & ! LCOV_EXCL_LINE
                                        , CholeskyDiago = ChoDia & ! LCOV_EXCL_LINE
                                        )
                    scaleFac = sqrt(Partition%Try(it)%ScaleFactorSq(ic))
                    ChoDia = ChoDia * scaleFac
                    do j = 1, nd
                        do i = j+1, nd
                            ChoLow(i,j) = ChoLow(i,j) * scaleFac
                        end do
                    end do
                    ChoLowDiff = 2 * abs(Partition%Try(it)%ChoLowCovUpp(1:nd,1:nd,ic)-ChoLow) / abs(Partition%Try(it)%ChoLowCovUpp(1:nd,1:nd,ic)+ChoLow)
                    if ( any(ChoLowDiff > TOLERANCE) ) then
                        ! LCOV_EXCL_START
                        write(*,*)
                        write(*,*) PROCEDURE_NAME
                        write(*,*) "Partition%Try(it)%ChoLowCovUpp(1:nd,1:nd,ic)-ChoLow) > TOLERANCE"
                        write(*,*) "TOLERANCE"
                        write(*,*) TOLERANCE
                        write(*,*) "ChoLowDiff"
                        write(*,*) ChoLowDiff
                        write(*,*) "ChoLow"
                        write(*,*) ChoLow
                        write(*,*) "Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic)"
                        write(*,*) Partition%Try(it)%ChoLowCovUpp(1:nd,1:nd,ic)
                        write(*,*)
                        error stop
                        ! LCOV_EXCL_STOP
                    end if
                    ChoDiaDiff = 2 * abs(Partition%Try(it)%ChoDia(1:nd,ic)-ChoDia) / abs(Partition%Try(it)%ChoDia(1:nd,ic)+ChoDia)
                    if ( any(ChoDiaDiff > TOLERANCE) ) then
                        ! LCOV_EXCL_START
                        write(*,*)
                        write(*,*) "Partition%Try(it)%ChoDia(1:nd,ic)-ChoDia) > TOLERANCE"
                        write(*,*) "TOLERANCE"
                        write(*,*) TOLERANCE
                        write(*,*) "ChoDiaDiff"
                        write(*,*) ChoDiaDiff
                        write(*,*) "ChoDia"
                        write(*,*) ChoDia
                        write(*,*) "Partition%Try(it)%ChoDia(1:nd,1:nd,ic)"
                        write(*,*) Partition%Try(it)%ChoDia(1:nd,ic)
                        write(*,*)
                        !error stop
                        ! LCOV_EXCL_STOP
                    end if
                    do ip = Partition%Try(it)%nps, Partition%Try(it)%npe
                        NormedPoint(1:nd) = Point(1:nd,ip) - Partition%Try(it)%Center(1:nd,ic)
                        mahalSqScalar = dot_product( NormedPoint , matmul(Partition%Try(it)%InvCovMat(1:nd,1:nd,ic), NormedPoint) )
                        if (mahalSqScalar<0._RK) then
                        ! LCOV_EXCL_START
                            mahalSqScalar = -1._RK
                            write(*,*) PROCEDURE_NAME
                            write(*,*) "mahalSqScalar<0._RK", ip, mahalSqScalar
                            error stop
                        end if
                        ! LCOV_EXCL_STOP
                        dum = getMahalSq( nd, Partition%Try(it)%Center(1:nd,ic), Partition%Try(it)%InvCovMat(1:nd,1:nd,ic), Point(1:nd,ip) )
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

!Partition%Try(it)%Size(Partition%Try(it)%ic) & ! LCOV_EXCL_LINE
!+ nint(Partition%inclusionFraction * & ! LCOV_EXCL_LINE
!( count(Partition%Try(it)%MahalSq(1:Partition%Try(it)%CumSumSize(Partition%Try(it)%ic-1),Partition%Try(it)%ic)<=Partition%Try(it)%ScaleFactorSq(Partition%Try(it)%ic)) & ! LCOV_EXCL_LINE
!+ count(Partition%Try(it)%MahalSq(Partition%Try(it)%CumSumSize(Partition%Try(it)%ic)+1:Partition%np,Partition%Try(it)%ic)<=Partition%Try(it)%ScaleFactorSq(Partition%Try(it)%ic)) & ! LCOV_EXCL_LINE
!) )
