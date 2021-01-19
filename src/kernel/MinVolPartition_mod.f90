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
!> This module contains procedures and routines for the minimum bounding ellipsoid partitioning of a given set of points.
!> \author
!> Amir Shahmoradi

module MinVolPartition_mod

    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type

    implicit none

    character(*), parameter :: MODULE_NAME = "@MinVolPartition_mod"

    !> The `MinVolPartition_type` class. Partitions an input array `Point(nd,np)` with `nd` attributes and `np` observations (points).
    type :: MinVolPartition_type
        integer(IK)                 :: nd,np                    !< See the interface of [runMinVolKernel()](@ref runminvolkernel).
        integer(IK)                 :: nemax                    !< See the interface of [runMinVolKernel()](@ref runminvolkernel).
        integer(IK)                 :: neopt                    !< See the interface of [runMinVolKernel()](@ref runminvolkernel).
        integer(IK)                 :: convergenceFailureCount  !< See the interface of [runMinVolKernel()](@ref runminvolkernel).
        integer(IK)                 :: maxKvolumeLoopRecursion  !< See the interface of [runMinVolKernel()](@ref runminvolkernel).
        integer(IK)                 :: maxAllowedKmeansFailure  !< See the interface of [runMinVolKernel()](@ref runminvolkernel).
        integer(IK)                 :: maxAllowedMinVolFailure  !< See the interface of [runMinVolKernel()](@ref runminvolkernel).
        logical                     :: stabilizationRequested   !< See the interface of [runMinVolKernel()](@ref runminvolkernel).
        real(RK)                    :: mahalSqWeightExponent    !< See the interface of [runMinVolKernel()](@ref runminvolkernel).
        real(RK)                    :: inclusionFraction        !< See the interface of [runMinVolKernel()](@ref runminvolkernel).
        real(RK)                    :: parLogVol                !< See the interface of [runMinVolKernel()](@ref runminvolkernel).
        real(RK)                    :: tightness                !< See the interface of [runMinVolKernel()](@ref runminvolkernel).
        integer(IK) , allocatable   :: Size(:)                  !< An array of size `(nemax)` representing the sizes of the corresponding bounding ellipsoids.
        integer(IK) , allocatable   :: Membership(:)            !< An array of size `(np)` representing the bounding-ellipsoid membership IDs of the corresponding data points.
        integer(IK) , allocatable   :: PointIndex(:)            !< An array of size `(np)` representing the indices of the original data points before partitioning.
        real(RK)    , allocatable   :: LogVolNormed(:)          !< An array of size `(nemax)` representing the log-volumes of the corresponding bounding ellipsoids.
        real(RK)    , allocatable   :: ChoLowCovUpp(:,:,:)      !< An array of size `(nd,nd,nemax)` representing the Cholesky lower triangle, diagonal, and covariance matrices of the bounding ellipsoids.
        real(RK)    , allocatable   :: InvCovMat(:,:,:)         !< An array of size `(nd,nd,nemax)` representing the full symmetric inverse covariance matrices of the bounding ellipsoids.
        real(RK)    , allocatable   :: ChoDia(:,:)              !< An array of size `(nd,nemax)` representing the Cholesky lower triangle, diagonal, and covariance matrices of the bounding ellipsoids.
        real(RK)    , allocatable   :: Center(:,:)              !< An array of size `(nd,nemax)` representing the centers of the bounding ellipsoids.
        type(Err_type)              :: Err                      !< An object of class [Err_type](@ref err_mod::err_type) containing information about error occurrence.
    contains
        procedure, pass :: write => writeMinVolPartitionOOP
    end type MinVolPartition_type

    interface MinVolPartition_type
        module procedure :: constructor
    end interface MinVolPartition_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructor( Point, nd, np, nemax & ! LCOV_EXCL_LINE
                        , maxKvolumeLoopRecursion & ! LCOV_EXCL_LINE
                        , maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                        , maxAllowedMinVolFailure & ! LCOV_EXCL_LINE
                        , stabilizationRequested & ! LCOV_EXCL_LINE
                        , mahalSqWeightExponent & ! LCOV_EXCL_LINE
                        , inclusionFraction & ! LCOV_EXCL_LINE
                        , parLogVol & ! LCOV_EXCL_LINE
                        , tightness & ! LCOV_EXCL_LINE
                        , trimEnabled & ! LCOV_EXCL_LINE
                        ) result(self)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructor
#endif
        use Constants_mod, only: IK, RK !, NEGINF_RK
        implicit none
        real(RK)    , intent(inout)         :: Point(nd,np)
        integer(IK) , intent(in)            :: nd, np
        integer(IK) , intent(in), optional  :: nemax
        integer(IK) , intent(in), optional  :: maxKvolumeLoopRecursion
        integer(IK) , intent(in), optional  :: maxAllowedKmeansFailure
        integer(IK) , intent(in), optional  :: maxAllowedMinVolFailure
        logical     , intent(in), optional  :: stabilizationRequested
        real(RK)    , intent(in), optional  :: mahalSqWeightExponent
        real(RK)    , intent(in), optional  :: inclusionFraction
        real(RK)    , intent(in), optional  :: parLogVol
        real(RK)    , intent(in), optional  :: tightness
        logical     , intent(in), optional  :: trimEnabled
        logical                             :: trimEnabledDefault
        type(MinVolPartition_type)          :: self

        self%nd = nd
        self%np = np
        self%Err%occurred = .false.
        trimEnabledDefault = .false.; if (present(trimEnabled)) trimEnabledDefault = trimEnabled

        self%nemax = np / (nd + 1)
        if (present(nemax)) then
            if (nemax > self%nemax) then
                self%Err%occurred = .true.
                self%Err%msg = "The input value for `nemax` cannot be larger than `np / (nd + 1)`."
                return
            else
                self%nemax = nemax
            end if
        end if

        self%maxKvolumeLoopRecursion = 50_IK; if (present(maxKvolumeLoopRecursion)) self%maxKvolumeLoopRecursion = maxKvolumeLoopRecursion
        self%maxAllowedKmeansFailure = 10_IK; if (present(maxAllowedKmeansFailure)) self%maxAllowedKmeansFailure = maxAllowedKmeansFailure
        self%maxAllowedMinVolFailure = 10_IK; if (present(maxAllowedMinVolFailure)) self%maxAllowedMinVolFailure = maxAllowedMinVolFailure
        self%stabilizationRequested = .true.; if (present(stabilizationRequested)) self%stabilizationRequested = stabilizationRequested
        self%mahalSqWeightExponent = 0._RK; if (present(mahalSqWeightExponent)) self%mahalSqWeightExponent = mahalSqWeightExponent
        self%inclusionFraction = 0._RK; if (present(inclusionFraction)) self%inclusionFraction = inclusionFraction
        self%parLogVol = -1.e-100_RK; if (present(parLogVol)) self%parLogVol = parLogVol
        self%tightness = 1.1_RK; if (present(tightness)) self%tightness = tightness

        allocate( self%Size         (self%nemax) & ! LCOV_EXCL_LINE
                , self%Center       (nd,self%nemax) & ! LCOV_EXCL_LINE
                , self%ChoDia       (nd,self%nemax) & ! LCOV_EXCL_LINE
                , self%InvCovMat    (nd,nd,self%nemax) & ! LCOV_EXCL_LINE
                , self%ChoLowCovUpp (nd,nd,self%nemax) & ! LCOV_EXCL_LINE
                , self%LogVolNormed (self%nemax) & ! LCOV_EXCL_LINE
                , self%Membership   (np) & ! LCOV_EXCL_LINE
                , self%PointIndex   (np) & ! LCOV_EXCL_LINE
                )

        call runMinVolPartition ( Point = Point & ! LCOV_EXCL_LINE
                                , nd = nd & ! LCOV_EXCL_LINE
                                , np = np & ! LCOV_EXCL_LINE
                                , nemax = self%nemax & ! LCOV_EXCL_LINE
                                , maxKvolumeLoopRecursion = self%maxKvolumeLoopRecursion & ! LCOV_EXCL_LINE
                                , maxAllowedKmeansFailure = self%maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                , maxAllowedMinVolFailure = self%maxAllowedMinVolFailure & ! LCOV_EXCL_LINE
                                , stabilizationRequested = self%stabilizationRequested & ! LCOV_EXCL_LINE
                                , mahalSqWeightExponent = self%mahalSqWeightExponent & ! LCOV_EXCL_LINE
                                , inclusionFraction = self%inclusionFraction & ! LCOV_EXCL_LINE
                                , tightness = self%tightness & ! LCOV_EXCL_LINE
                                , parLogVol = self%parLogVol & ! LCOV_EXCL_LINE
                                , neopt = self%neopt & ! LCOV_EXCL_LINE
                                , PointIndex = self%PointIndex & ! LCOV_EXCL_LINE
                                , ClusterSize = self%Size & ! LCOV_EXCL_LINE
                                , ClusterCenter = self%Center & ! LCOV_EXCL_LINE
                                , ClusterChoDia = self%ChoDia & ! LCOV_EXCL_LINE
                                , ClusterLogVolume = self%LogVolNormed & ! LCOV_EXCL_LINE
                                , ClusterInvCovMat = self%InvCovMat & ! LCOV_EXCL_LINE
                                , ClusterMembership = self%Membership & ! LCOV_EXCL_LINE
                                , ClusterChoLowCovUpp = self%ChoLowCovUpp & ! LCOV_EXCL_LINE
                                , convergenceFailureCount = self%convergenceFailureCount & ! LCOV_EXCL_LINE
                                , Err = self%Err & ! LCOV_EXCL_LINE
                                )

        if (trimEnabledDefault) then
            self%Size           = self%Size         (1:self%neopt)
            self%Center         = self%Center       (1:nd,1:self%neopt)
            self%ChoDia         = self%ChoDia       (1:nd,1:self%neopt)
            self%InvCovMat      = self%InvCovMat    (1:nd,1:nd,1:self%neopt)
            self%ChoLowCovUpp   = self%ChoLowCovUpp (1:nd,1:nd,1:self%neopt)
            self%LogVolNormed   = self%LogVolNormed (1:self%neopt)
        end if

    end function constructor

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine runMinVolPartition   ( Point & ! LCOV_EXCL_LINE
                                    , nd,np & ! LCOV_EXCL_LINE
                                    , nemax & ! LCOV_EXCL_LINE
                                    , maxKvolumeLoopRecursion & ! LCOV_EXCL_LINE
                                    , maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                    , maxAllowedMinVolFailure & ! LCOV_EXCL_LINE
                                    , stabilizationRequested & ! LCOV_EXCL_LINE
                                    , mahalSqWeightExponent & ! LCOV_EXCL_LINE
                                    , inclusionFraction & ! LCOV_EXCL_LINE
                                    , tightness & ! LCOV_EXCL_LINE
                                    , parLogVol & ! LCOV_EXCL_LINE
                                    ! output
                                    , neopt & ! LCOV_EXCL_LINE
                                    , PointIndex & ! LCOV_EXCL_LINE
                                    , ClusterSize & ! LCOV_EXCL_LINE
                                    , ClusterCenter & ! LCOV_EXCL_LINE
                                    , ClusterChoDia & ! LCOV_EXCL_LINE
                                    , ClusterLogVolume & ! LCOV_EXCL_LINE
                                    , ClusterInvCovMat & ! LCOV_EXCL_LINE
                                    , ClusterMembership & ! LCOV_EXCL_LINE
                                    , ClusterChoLowCovUpp & ! LCOV_EXCL_LINE
                                    , convergenceFailureCount & ! LCOV_EXCL_LINE
                                    , Err & ! LCOV_EXCL_LINE
                                    )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: runMinVolPartition
#endif
        use Statistics_mod, only: getSamCholFac, getMahalSq
        use Constants_mod, only: IK, RK, NEGINF_RK
        use Matrix_mod, only: getInvMatFromCholFac
        use Matrix_mod, only: getCholeskyFactor
        use String_mod, only: num2str
        use Err_mod, only: Err_type

        implicit none

        real(RK)        , intent(inout) :: Point(nd,np)                 ! the input data (points) to be clustered. on output, it is ordered.
        integer(IK)     , intent(in)    :: nd,np                        ! the maximum allowed number of ellipsoidal clusters, number of dimensions, number of points
        integer(IK)     , intent(in)    :: nemax                        ! the maximum allowed number of ellipsoidal clusters, number of dimensions, number of points
        integer(IK)     , intent(in)    :: maxKvolumeLoopRecursion
        integer(IK)     , intent(in)    :: maxAllowedKmeansFailure
        integer(IK)     , intent(in)    :: maxAllowedMinVolFailure
        logical         , intent(in)    :: stabilizationRequested
        real(RK)        , intent(in)    :: mahalSqWeightExponent
        real(RK)        , intent(in)    :: inclusionFraction                ! the fraction of non-cluster-member points that ARE inside the cluster, to be considered in the estimation of the true volume of the cluster.
        real(RK)        , intent(in)    :: tightness                        ! the factor deremining how large the bounding ellipsoid can be, compared to the true volume of the ellipsoid.
        real(RK)        , intent(in)    :: parLogVol                        ! input log-estimate of the total volume of the points.
        integer(IK)     , intent(out)   :: neopt                            ! this variable is saved throughout recursive iterations of subroutine
        integer(IK)     , intent(out)   :: PointIndex(np)                   ! Vector containing the scrambled indices of Point elements. It represent the mapping from the input Point elements to output scrambled Point vector.
        integer(IK)     , intent(out)   :: ClusterSize(nemax)               ! cluster sizes (nemax)
        real(RK)        , intent(out)   :: ClusterLogVolume(nemax)             ! CBE: Cluster Bounding Ellipsoid - On output: Each element of the array contains the volume of the corresponding Cluster Bounding Ellipsoid.
        integer(IK)     , intent(out)   :: ClusterMembership(np)            ! the ClusterMembership id of each point (representing the cluster number to which the point belongs)
        real(RK)        , intent(out)   :: ClusterCenter(nd,nemax)          ! Centers of the clusters (nd,nemax)
        real(RK)        , intent(out)   :: ClusterChoDia(nd,nemax)          ! CBE: Cluster Bounding Ellipsoid - On output: The diagonal elements of the Cholesky factor of the covariance matrix of the input points, on output, contains the cholesky factor of all found clusters.
        real(RK)        , intent(out)   :: ClusterInvCovMat(nd,nd,nemax)    ! CBE: Cluster Bounding Ellipsoid - On output: Full symmetric Inverse Covariance Matrices of the Bounding Ellipsoids of all clusters.
        real(RK)        , intent(out)   :: ClusterChoLowCovUpp(nd,nd,nemax) ! CBE: Cluster Bounding Ellipsoid - On output: The lower triangle contains the Cholesky factor of the covariance matrix of the input points, the upper triangle maybe the CovMat itself (not needed). On output contains in the lower triangle, the Cholesky factor of all clusters found.
        integer(IK)     , intent(inout) :: convergenceFailureCount          ! number of times the algorithm failed to converge.
        type(Err_type)  , intent(out)   :: Err
        real(RK)                        :: NormedPoint(nd,np)
        real(RK)                        :: npMinusOneInverse
        real(RK)                        :: scaleFactorSq
        real(RK)                        :: parVolNormed                     ! The normalized volume of parent ellipsoid
        real(RK)                        :: scaleFactor
        real(RK)                        :: logNormFac
        real(RK)                        :: mahalSq
        integer                         :: i, j, ip

        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME//"@runMinVolPartition()"

        Err%occurred = .false.

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Initialize the first cluster's center
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ClusterCenter(1:nd,1) = 0._RK
        do ip = 1, np
          PointIndex(ip) = ip
          ClusterCenter(1:nd,1) = ClusterCenter(1:nd,1) + Point(1:nd,ip)
        end do
        ClusterCenter(1:nd,1) = ClusterCenter(1:nd,1) / real(np,kind=RK)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the covariance matrix and the Cholesky factorization of the first cluster.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        npMinusOneInverse = 1._RK / real(np-1,kind=RK)

        ! Compute the covariance matrix upper.

        ClusterChoLowCovUpp(1:nd,1:nd,1) = 0._RK
        do ip = 1, np
            do j = 1, nd
                NormedPoint(j,ip) = Point(j,ip) - ClusterCenter(j,1)
                ClusterChoLowCovUpp(1:j,j,1) = ClusterChoLowCovUpp(1:j,j,1) + NormedPoint(1:j,ip) * NormedPoint(j,ip)
            end do
        end do
        ClusterChoLowCovUpp(1:nd,1:nd,1) = ClusterChoLowCovUpp(1:nd,1:nd,1) * npMinusOneInverse

        ! Compute the Cholesky factor lower.

        call getCholeskyFactor(nd, ClusterChoLowCovUpp(1:nd,1:nd,1), ClusterChoDia(1:nd,1))
        if (ClusterChoDia(1,1)<0._RK) then
            Err%msg = PROCEDURE_NAME//": Singular covariance matrix detected for nemax, nd, np: "//num2str(nemax)//", "//num2str(nd)//", "//num2str(np)
            write(*,*) 'Max cluster number, dimension, total points:'
            write(*,*) nemax, nd, np
            error stop
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the inverse covariance matrix.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ClusterInvCovMat(1:nd,1:nd,1) = getInvMatFromCholFac(nd, ClusterChoLowCovUpp(1:nd,1:nd,1), ClusterChoDia(1:nd,1))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the Mahal distances squared and the scale factor which makes the covariance matrix a bounding ellipsoid.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        scaleFactorSq = NEGINF_RK
        do concurrent(ip = 1:np)
            mahalSq = dot_product( NormedPoint(1:nd,ip) , matmul(ClusterInvCovMat(1:nd,1:nd,1), NormedPoint(1:nd,ip)) )
            if (scaleFactorSq < mahalSq) scaleFactorSq = mahalSq
        end do
        scaleFactor = sqrt(scaleFactorSq)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the largest possible LogVolume of the first cluster based on the parent cluster and the input estimated volume.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ClusterLogVolume(1) = sum( log( ClusterChoDia(1:nd,1) ) ) + nd * log(scaleFactor)

        if (ClusterLogVolume(1) > parLogVol) then
            ClusterLogVolume(1) = 1._RK
            logNormFac = ClusterLogVolume(1)
            parVolNormed = exp(parLogVol - logNormFac)
        else    ! enlarge the bounding ellipsoid
            ClusterLogVolume(1) = 1._RK
            logNormFac = parLogVol
            parVolNormed = 1._RK
            !scaleFactor = scaleFactor * exp( (logNormFac-ClusterLogVolume(1))/nd )
            !scaleFactorSq = scaleFactor**2
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX What's going on here?
        end if

        if (nemax > 1 .and. ClusterLogVolume(1) > tightness * parVolNormed) then
            convergenceFailureCount = 0
            call runMinVolKernel( Point = Point & ! LCOV_EXCL_LINE
                                , nd = nd & ! LCOV_EXCL_LINE
                                , np = np & ! LCOV_EXCL_LINE
                                , nemax = nemax & ! LCOV_EXCL_LINE
                                , maxKvolumeLoopRecursion = maxKvolumeLoopRecursion & ! LCOV_EXCL_LINE
                                , maxAllowedKmeansFailure = maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                , maxAllowedMinVolFailure = maxAllowedMinVolFailure & ! LCOV_EXCL_LINE
                                , stabilizationRequested = stabilizationRequested & ! LCOV_EXCL_LINE
                                , mahalSqWeightExponent = mahalSqWeightExponent & ! LCOV_EXCL_LINE
                                , inclusionFraction = inclusionFraction & ! LCOV_EXCL_LINE
                                , parVolNormed = parVolNormed & ! LCOV_EXCL_LINE
                                , logNormFac = logNormFac & ! LCOV_EXCL_LINE
                                , tightness = tightness & ! LCOV_EXCL_LINE
                                , neopt = neopt & ! LCOV_EXCL_LINE
                                , PointIndex = PointIndex & ! LCOV_EXCL_LINE
                                , ClusterSize = ClusterSize & ! LCOV_EXCL_LINE
                                , ClusterCenter = ClusterCenter & ! LCOV_EXCL_LINE
                                , ClusterChoDia = ClusterChoDia & ! LCOV_EXCL_LINE
                                , ClusterLogVolume = ClusterLogVolume & ! LCOV_EXCL_LINE
                                , ClusterInvCovMat = ClusterInvCovMat & ! LCOV_EXCL_LINE
                                , ClusterMembership = ClusterMembership & ! LCOV_EXCL_LINE
                                , ClusterChoLowCovUpp = ClusterChoLowCovUpp & ! LCOV_EXCL_LINE
                                , convergenceFailureCount = convergenceFailureCount & ! LCOV_EXCL_LINE
                                )
        else ! There is only one cluster
            neopt = 1_IK
        end if

        if (neopt==1_IK) then

            ! Now rescale Cholesky factor and its diagonals, such that they describe the bounding ellipsoid of the cluster.

            ClusterSize(1) = np
            ClusterMembership(1:np) = 1_IK ! @todo : xxx this is redundant and should be removed in the future.
            do j = 1, nd - 1
                ClusterChoDia(j,1) = ClusterChoDia(j,1) * scaleFactor
                !ClusterChoLowCovUpp(1:j,j,1) = ClusterChoLowCovUpp(1:j,j,1) * scaleFactorSq
                do i = j+1, nd
                    ClusterChoLowCovUpp(i,j,1) = ClusterChoLowCovUpp(i,j,1) * scaleFactor
                end do
            end do
            !ClusterChoLowCovUpp(1:nd,nd,1) = ClusterChoLowCovUpp(1:nd,nd,1) * scaleFactorSq
            ClusterChoDia(nd,1) = ClusterChoDia(nd,1) * scaleFactor
            ClusterInvCovMat(1:nd,1:nd,1) = ClusterInvCovMat(1:nd,1:nd,1) / scaleFactorSq

        end if

    end subroutine runMinVolPartition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This code assumes that nemax > 1 on input.
    ! important input requirements:
    recursive subroutine runMinVolKernel( Point & ! LCOV_EXCL_LINE
                                        , nd & ! LCOV_EXCL_LINE
                                        , np & ! LCOV_EXCL_LINE
                                        , nemax & ! LCOV_EXCL_LINE
                                        , maxKvolumeLoopRecursion & ! LCOV_EXCL_LINE
                                        , maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                        , maxAllowedMinVolFailure & ! LCOV_EXCL_LINE
                                        , stabilizationRequested & ! LCOV_EXCL_LINE
                                        , mahalSqWeightExponent & ! LCOV_EXCL_LINE
                                        , inclusionFraction & ! LCOV_EXCL_LINE
                                        , parVolNormed & ! LCOV_EXCL_LINE
                                        , logNormFac & ! LCOV_EXCL_LINE
                                        , tightness & ! LCOV_EXCL_LINE
                                        , neopt & ! LCOV_EXCL_LINE
                                        , PointIndex & ! LCOV_EXCL_LINE
                                        , ClusterSize & ! LCOV_EXCL_LINE
                                        , ClusterCenter & ! LCOV_EXCL_LINE
                                        , ClusterChoDia & ! LCOV_EXCL_LINE
                                        , ClusterLogVolume & ! LCOV_EXCL_LINE
                                        , ClusterInvCovMat & ! LCOV_EXCL_LINE
                                        , ClusterMembership & ! LCOV_EXCL_LINE
                                        , ClusterChoLowCovUpp & ! LCOV_EXCL_LINE
                                        , convergenceFailureCount & ! LCOV_EXCL_LINE
                                        )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: runMinVolKernel
#endif
        use Statistics_mod, only: getSamCholFac, getMahalSq
        use Matrix_mod, only: getInvMatFromCholFac
        use Kmeans_mod, only: runKmeans
        use Err_mod, only: Err_type
        implicit none
        integer(IK) , save          :: ncall = 0
        integer(IK) , intent(in)    :: nemax                             ! maximum number of clusters determined from the nd and np
        integer(IK) , intent(in)    :: nd,np                             ! number of dimensions, number of points
        integer(IK) , intent(in)    :: maxKvolumeLoopRecursion
        integer(IK) , intent(in)    :: maxAllowedKmeansFailure
        integer(IK) , intent(in)    :: maxAllowedMinVolFailure
        logical     , intent(in)    :: stabilizationRequested
        real(RK)    , intent(in)    :: mahalSqWeightExponent
        real(RK)    , intent(in)    :: tightness                        ! the factor deremining how large the bounding ellipsoid can be, compared to the true volume of the ellipsoid.
        real(RK)    , intent(in)    :: inclusionFraction                ! the fraction of non-cluster-member points that ARE inside the cluster, to be considered in the estimation of the true volume of the cluster.
        real(RK)    , intent(in)    :: parVolNormed                     ! parent's sqrt(det(invCovMat)), where invCovMat is the inverse covariance matrix of the bounding ellipsoid of the input points.
        real(RK)    , intent(in)    :: logNormFac                       ! parent's sqrt(det(invCovMat)), where invCovMat is the inverse covariance matrix of the bounding ellipsoid of the input points.
        integer(IK) , intent(out)   :: neopt                            ! The predicted optimal number of clusters identified in the input data.
        integer(IK) , intent(out)   :: ClusterSize(nemax)                     ! cluster sizes (nemax)
        integer(IK) , intent(inout) :: convergenceFailureCount          ! returns the number of times the algorithm failed to converge.
        integer(IK) , intent(inout) :: ClusterMembership(np)                  ! the ClusterMembership id of each point (representing the cluster number to which the point belongs)
        integer(IK) , intent(inout) :: PointIndex(np)                   ! Vector containing the scrambled indices of Point elements. It represent the mapping from the input Point elements to output scrambled Point vector.
        real(RK)    , intent(inout) :: Point(nd,np)                     ! the input data (points) to be clustered.
        real(RK)    , intent(inout) :: ClusterCenter(nd,nemax)                ! Centers of the clusters (nd,nemax)
        real(RK)    , intent(inout) :: ClusterChoLowCovUpp(nd,nd,nemax)     ! CBE: Cluster Bounding Ellipsoid - On input: The lower triangle contains the Cholesky factor of the covariance matrix of the input points, the upper triangle maybe the CovMat itself (not needed). On output contains in the lower triangle, the Cholesky factor of all clusters found.
        real(RK)    , intent(inout) :: ClusterInvCovMat(nd,nd,nemax)        ! CBE: Cluster Bounding Ellipsoid - On input: The lower triangle contains the Cholesky factor of the covariance matrix of the input points, the upper triangle maybe the CovMat itself (not needed). On output contains in the lower triangle, the Cholesky factor of all clusters found.
        real(RK)    , intent(inout) :: ClusterChoDia(nd,nemax)              ! CBE: Cluster Bounding Ellipsoid - On input: The diagonal elements of the Cholesky factor of the covariance matrix of the input points, on output, contains the cholesky factor of all found clusters.
        real(RK)    , intent(inout) :: ClusterLogVolume(nemax)               ! CBE: Cluster Bounding Ellipsoid - On input: The volume of the Bounding Ellipsoid of the parent cluster. On output: The Bounding volume of all child clusters.

        real(RK)                    :: CKCenter(nd,2)                    ! Cluster Centers for nc=2 Kmeans algorithm
        real(RK)                    :: CKBEChoLowCovUpp(nd,nd,2)              ! Covariance Matrices of the bounding ellipsoid of the Kmeans' 2 clusters, reported in the upper triangle. The lower triangle contains the Cholesky factors
        real(RK)                    :: CKBEInvCovMat(nd,nd,2)            ! Inverse of the Covariance Matrices of the bounding ellipsoid of the Kmeans' 2 clusters.
        real(RK)                    :: CKBEChoDia(nd,2)                ! determinants of the inverse covariance matrices of the bounding ellipsoids of the Kmeans' 2 clusters.
        real(RK)                    :: CKMahalSq(np,2)                   ! Vector of Maximum Mahalanobis distances for each cluster
        real(RK)                    :: CKBEVolNormed(2)                  ! New Cluster Bounding Ellisoid Volumes
        real(RK)                    :: CKParVolNormed(2)                 ! New Cluster Volume estimates
        real(RK)                    :: ScaleFacSq(2),scaleFactor(2)         ! The factor by which determinant must be multiplied in order for the ellipsoid to become a bounding ellipsoid.
        real(RK)                    :: CKPoint(nd,np)                    ! Cluster Points in 2-means partition result
        real(RK)                    :: MahalSqWeight(2)                  ! MahalSqWeight definition of Lu et al 2007
        integer(IK)                 :: CKPointIndex(np)                  ! Vector containing the scrambled indices of Point elements. It represent the mapping from the input Point elements to output scrambled Point vector.
        integer(IK)                 :: CKSize(2)                         ! cluster sizes for nc=2 Kmeans algorithm
        integer(IK)                 :: NumPointInCK(2)                   ! cluster sizes for nc=2 Kmeans algorithm
        integer(IK)                 :: CKnemax(2)                        ! maximum allowed number of clusters in each of the two K-means clusters.
        integer(IK)                 :: CKneopt(2)                    ! optimal number of clusters in each of the child clusters
        integer(IK)                 :: CMemberCounter(2)                 ! cluster member counter
        integer(IK)                 :: zcsIteration                      ! Counts the number of times zero sized clusters are found, zcs stands for zero cluster size
        integer(IK)                 :: IpStart(2),IpEnd(2)
        integer(IK)                 :: i,j,ip,ic,icOther,iteration,ndEffective
        integer(IK)                 :: icStart,icEnd,recursionCounter,kmeansFailureCounter,minVolFailureCounter
        integer(IK)                 :: minClusterSize
        logical                     :: reclusteringNeeded

        ! stabilization variables:

        real(RK)                :: CKCenterOld(nd,2)            ! Cluster Centers for nc=2 Kmeans algorithm
        real(RK)                :: CKBEChoLowCovUppOld(nd,nd,2)      ! Covariance Matrices of the bounding ellipsoid of the Kmeans' 2 clusters, reported in the upper triangle. The lower triangle contains the Cholesky factors
        real(RK)                :: CKBEInvCovMatOld(nd,nd,2)    ! Inverse of the Covariance Matrices of the bounding ellipsoid of the Kmeans' 2 clusters.
        real(RK)                :: CKBEChoDiaOld(nd,2)        ! determinants of the inverse covariance matrices of the bounding ellipsoids of the Kmeans' 2 clusters.
        real(RK)                :: CKNormedPoint(nd,np,2)       ! Vector of Maximum Mahalanobis distances for each cluster
        real(RK)                :: CKMinDistanceSq(np)          ! Vector of Maximum Mahalanobis distances for each cluster
        real(RK)                :: CKMahalSqOld(np,2)           ! Vector of Maximum Mahalanobis distances for each cluster
        real(RK)                :: CKBEVolNormedOld(2)          ! New Cluster Bounding Ellisoid Volumes
        real(RK)                :: ScaleFacSqOld(2)             ! The factor by which determinant must be multiplied in order for the ellipsoid to become a bounding ellipsoid.
        real(RK)                :: ScaleFacOld(2)               ! The factor by which determinant must be multiplied in order for the ellipsoid to become a bounding ellipsoid.
        real(RK)                :: potential                    ! The factor by which determinant must be multiplied in order for the ellipsoid to become a bounding ellipsoid.
        type(Err_type)          :: Err

        minClusterSize = nd + 1
        ndEffective = np / nemax + 1
        ncall = ncall + 1

        kmeansFailureCounter = 0
        loopRunKmeans: do
            call runKmeans  ( nd = nd & ! LCOV_EXCL_LINE
                            , np = np & ! LCOV_EXCL_LINE
                            , nc = 2_IK & ! LCOV_EXCL_LINE
                            , Point = Point & ! LCOV_EXCL_LINE
                            , Size = CKSize & ! LCOV_EXCL_LINE
                            , Center = CKCenter & ! LCOV_EXCL_LINE
                            , Membership = ClusterMembership & ! LCOV_EXCL_LINE
                            , NormedPoint = CKNormedPoint & ! LCOV_EXCL_LINE
                            , MinDistanceSq = CKMinDistanceSq & ! LCOV_EXCL_LINE
                            , potential = potential & ! LCOV_EXCL_LINE
                            , niter = iteration & ! LCOV_EXCL_LINE
                            , nzsci = zcsIteration & ! LCOV_EXCL_LINE
                            , Err = Err & ! LCOV_EXCL_LINE
                            !, niterMax & ! LCOV_EXCL_LINE
                            !, nzsciMax & ! LCOV_EXCL_LINE
                            !, relTol & ! LCOV_EXCL_LINE
                            )
            if ( any(CKSize<minClusterSize) ) then
                kmeansFailureCounter = kmeansFailureCounter + 1
                if (kmeansFailureCounter>maxAllowedKmeansFailure) then
                    neopt = 1
                    ClusterMembership = 1
                    ClusterSize(1) = np
                    convergenceFailureCount = convergenceFailureCount + 1
                    return
                end if
                cycle loopRunKmeans
            end if
            exit loopRunKmeans
        end do loopRunKmeans

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CKBEVolNormedOld = huge(0._RK)
        recursionCounter = 0
        minVolFailureCounter = 0

        loopRecursiveRunMinVolKernel: do    ! find two approximately minimum-volume sub-cllusters recursively

            recursionCounter = recursionCounter + 1
            if(mod(recursionCounter,maxKvolumeLoopRecursion)==0) then
                write(*,*) 'recursionCounter: ', recursionCounter
                neopt = 1
                ClusterMembership = 1
                ClusterSize(1) = np
                convergenceFailureCount = convergenceFailureCount + 1
                return
            end if

            CMemberCounter(1) = 0
            CMemberCounter(2) = CKSize(1)
            do ip = 1,np
                if (ClusterMembership(ip)==1) then
                    CMemberCounter(1) = CMemberCounter(1) + 1
                    CKPoint(1:nd,CMemberCounter(1)) = Point(1:nd,ip)
                    CKPointIndex(CMemberCounter(1)) = PointIndex(ip)
                else
                    CMemberCounter(2) = CMemberCounter(2) + 1
                    CKPoint(1:nd,CMemberCounter(2)) = Point(1:nd,ip)
                    CKPointIndex(CMemberCounter(2)) = PointIndex(ip)
                end if
            end do
            Point = CKPoint           ! Point is now ordered
            PointIndex = CKPointIndex ! PointIndex is now ordered

            IpStart(1) = 1             ; IpEnd(1) = CKSize(1)
            IpStart(2) = CKSize(1) + 1 ; IpEnd(2) = np

            do ic = 1,2   ! Get scaled Mahalanobis distances

                icOther = 3 - ic

                ! get Cholesky Factor of cluster samples

                call getSamCholFac  ( nd &
                                    , CKSize(ic) &
                                    , CKCenter(1:nd,ic) &
                                    , Point(1:nd,IpStart(ic):IpEnd(ic)) &
                                    , CKBEChoLowCovUpp(1:nd,1:nd,ic) &
                                    , CKBEChoDia(1:nd,ic) &
                                    )
                if (CKBEChoDia(1,ic)<0._RK) then
                    write(*,*) "    FATAL: runMinVolKernel() @ getSamCholFac() @ getCholeskyFactor() failed."
                    write(*,*) "           Cluster id, dimension, size:"
                    write(*,*) ic, nd, CKSize(ic)
                    stop
                end if

                CKBEInvCovMat(1:nd,1:nd,ic) = getInvMatFromCholFac( nd, CKBEChoLowCovUpp(1:nd,1:nd,ic), CKBEChoDia(1:nd,ic) )
                CKMahalSq(1:np,ic) = getMahalSq( nd, np, CKCenter(1:nd,ic), CKBEInvCovMat(1:nd,1:nd,ic), Point )
                ScaleFacSq(ic) = maxval(CKMahalSq(IpStart(ic):IpEnd(ic),ic))
                scaleFactor(ic) = sqrt( ScaleFacSq(ic) )
                CKMahalSq(1:np,ic) = CKMahalSq(1:np,ic) / ScaleFacSq(ic)
                CKBEVolNormed(ic) = exp( nd*log(scaleFactor(ic)) + sum(log(CKBEChoDia(1:nd,ic))) - logNormFac ) ! Get the LogVolume of the cluster
                NumPointInCK(ic) = CKSize(ic) + nint( inclusionFraction * count(CKMahalSq(IpStart(icOther):IpEnd(icOther),ic)<1._RK), kind = IK )
                CKParVolNormed(ic) = real(NumPointInCK(ic),kind=RK) * parVolNormed / real(np,kind=RK)

                blockStabilityCheck: if (stabilizationRequested) then

                    if ( CKBEVolNormed(ic) > CKBEVolNormedOld(ic) ) then    ! This is to ensure the convergence of the algorithm

                        CKCenter(1:nd,ic) = CKCenterOld(1:nd,ic)
                        CKBEChoLowCovUpp(1:nd,1:nd,ic) = CKBEChoLowCovUppOld(1:nd,1:nd,ic)
                        CKBEChoDia(1:nd,ic) = CKBEChoDiaOld(1:nd,ic)
                        CKBEInvCovMat(1:nd,1:nd,ic) = CKBEInvCovMatOld(1:nd,1:nd,ic)
                        ScaleFacSq(ic) = ScaleFacSqOld(ic)
                        scaleFactor(ic) = ScaleFacOld(ic)
                        CKBEVolNormed(ic) = CKBEVolNormedOld(ic)

                        CMemberCounter(1) = 0
                        CMemberCounter(2) = CKSize(1)
                        do ip = 1,np
                            if (ClusterMembership(ip)==1) then
                                CMemberCounter(1) = CMemberCounter(1) + 1
                                CKMahalSq(CMemberCounter(1),ic) = CKMahalSqOld(ip,ic)
                            else
                                CMemberCounter(2) = CMemberCounter(2) + 1
                                CKMahalSq(CMemberCounter(2),ic) = CKMahalSqOld(ip,ic)
                            end if
                        end do
                        CKMahalSqOld(1:np,ic) = CKMahalSq(1:np,ic)

                    else  ! no stability correction is required. Everything just fine...

                        CKCenterOld(1:nd,ic) = CKCenter(1:nd,ic)
                        CKBEChoLowCovUppOld(1:nd,1:nd,ic) = CKBEChoLowCovUpp(1:nd,1:nd,ic)
                        CKBEChoDiaOld(1:nd,ic) = CKBEChoDia(1:nd,ic)
                        CKBEInvCovMatOld(1:nd,1:nd,ic) = CKBEInvCovMat(1:nd,1:nd,ic)
                        CKMahalSqOld(1:np,ic) = CKMahalSq(1:np,ic)
                        ScaleFacSqOld(ic) = ScaleFacSq(ic)
                        ScaleFacOld(ic) = scaleFactor(ic)
                        CKBEVolNormedOld(ic) = CKBEVolNormed(ic)

                    end if

                end if blockStabilityCheck

                !if (CKBEVolNormed(ic)<CKParVolNormed(ic)) then  ! enlarge the bounding ellipsoid
                !  !write(*,*) 'enlargement happening...'
                !  scaleFacSq(ic) = scaleFacSq(ic) * (CKParVolNormed(ic)/CKBEVolNormed(ic))**(2._RK/nd)  ! ATTN: (1._RK/nd) this should likely be (2._RK/nd)
                !  scaleFactor(ic) = sqrt(scaleFacSq(ic))
                !  CKBEVolNormed(ic) = CKParVolNormed(ic)
                !end if

#if defined DEBUG_ENABLED
                debugVolNorm: block
                    real(RK) :: dummy
                    dummy = real(CKSize(ic),kind=RK)*parVolNormed/real(np,kind=RK)
                    if ( CKBEVolNormed(ic) < dummy ) then  ! enlarge the bounding ellipsoid
                        !write(*,*) 'enlargement happening...'
                        scaleFacSq(ic) = scaleFacSq(ic) * (dummy/CKBEVolNormed(ic))**(2._RK/nd)  ! ATTN: (1._RK/nd) this should likely be (2._RK/nd)
                        scaleFactor(ic) = sqrt(scaleFacSq(ic))
                        CKBEVolNormed(ic) = dummy
                    end if
                end block debugVolNorm
#endif
              !end if blockStabilityCheck

            end do

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! Reassign points to clusters if needed
            MahalSqWeight = ( CKBEVolNormed / CKParVolNormed )**mahalSqWeightExponent
            reclusteringNeeded = .false.
            CKCenter(1:nd,1) = CKCenter(1:nd,1) * CKSize(1)  ! Now it is sum instead of mean
            CKCenter(1:nd,2) = CKCenter(1:nd,2) * CKSize(2)  ! Now it is sum instead of mean

            ic = 1
            icOther = 2
            do ip = IpStart(ic), IpEnd(ic)
                ClusterMembership(ip) = ic    ! This has not been assigned correctly in the past, here is the first attempt
                if ( MahalSqWeight(ic)*CKMahalSq(ip,ic) >= MahalSqWeight(icOther)*CKMahalSq(ip,icOther) ) then
                    ClusterMembership(ip) = icOther
                    CKSize(ic) = CKSize(ic) - 1
                    CKSize(icOther) = CKSize(icOther) + 1
                    CKCenter(1:nd,ic) = CKCenter(1:nd,ic) - Point(1:nd,ip)
                    CKCenter(1:nd,icOther) = CKCenter(1:nd,icOther) + Point(1:nd,ip)
                    reclusteringNeeded = .true.
                end if
            end do

            ic = 2
            icOther = 1
            do ip = IpStart(ic), IpEnd(ic)
                ClusterMembership(ip) = ic    ! This has not been assigned correctly in the past, here is the first attempt
                if ( MahalSqWeight(ic)*CKMahalSq(ip,ic) > MahalSqWeight(icOther)*CKMahalSq(ip,icOther) ) then
                    ClusterMembership(ip) = icOther
                    CKSize(ic) = CKSize(ic) - 1
                    CKSize(icOther) = CKSize(icOther) + 1
                    CKCenter(1:nd,ic) = CKCenter(1:nd,ic) - Point(1:nd,ip)
                    CKCenter(1:nd,icOther) = CKCenter(1:nd,icOther) + Point(1:nd,ip)
                    reclusteringNeeded = .true.
                end if
            end do

            blockReclusteringNeeded: if (reclusteringNeeded) then  ! perform reassignment

                blockIllegalClusterSize: if ( any(CKSize<minClusterSize) ) then

                    minVolFailureCounter = minVolFailureCounter + 1
                    if (minVolFailureCounter>maxAllowedMinVolFailure) then
                        neopt = 1
                        ClusterMembership = 1
                        ClusterSize(1) = np
                        convergenceFailureCount = convergenceFailureCount + 1
                        !debugMinVolFailureCounter: block
                        !  write(*,*) 'minVolFailureCounter: ', minVolFailureCounter
                        !end block debugMinVolFailureCounter
                        return
                    end if

                    kmeansFailureCounter = 0
                    loopRerunKmeans: do
                        call runKmeans  ( nd = nd & ! LCOV_EXCL_LINE
                                        , np = np & ! LCOV_EXCL_LINE
                                        , nc = 2_IK & ! LCOV_EXCL_LINE
                                        , Point = Point & ! LCOV_EXCL_LINE
                                        , Size = CKSize & ! LCOV_EXCL_LINE
                                        , Center = CKCenter & ! LCOV_EXCL_LINE
                                        , Membership = ClusterMembership & ! LCOV_EXCL_LINE
                                        , NormedPoint = CKNormedPoint & ! LCOV_EXCL_LINE
                                        , MinDistanceSq = CKMinDistanceSq & ! LCOV_EXCL_LINE
                                        , potential = potential & ! LCOV_EXCL_LINE
                                        , niter = iteration & ! LCOV_EXCL_LINE
                                        , nzsci = zcsIteration & ! LCOV_EXCL_LINE
                                        , Err = Err & ! LCOV_EXCL_LINE
                                        !, niterMax & ! LCOV_EXCL_LINE
                                        !, nzsciMax & ! LCOV_EXCL_LINE
                                        !, relTol & ! LCOV_EXCL_LINE
                                        )
                        if ( any(CKSize<minClusterSize) ) then
                        kmeansFailureCounter = kmeansFailureCounter + 1
                        if (kmeansFailureCounter>maxAllowedKmeansFailure) then
                            neopt = 1
                            ClusterMembership = 1
                            ClusterSize(1) = np
                            convergenceFailureCount = convergenceFailureCount + 1
                            !debugKmeansFailureCounter: block
                            !  write(*,*) 'kmeansFailureCounter: ', kmeansFailureCounter
                            !end block debugKmeansFailureCounter
                            return
                        end if
                        cycle loopRerunKmeans
                        end if
                        exit loopRerunKmeans
                    end do loopRerunKmeans

                    CKBEVolNormedOld = huge(0._RK)
                    cycle loopRecursiveRunMinVolKernel

                end if blockIllegalClusterSize

                CKCenter(1:nd,1) = CKCenter(1:nd,1) / CKSize(1)  ! Now it is sum instead of mean
                CKCenter(1:nd,2) = CKCenter(1:nd,2) / CKSize(2)  ! Now it is sum instead of mean

                cycle loopRecursiveRunMinVolKernel

            end if blockReclusteringNeeded

            CKCenter(1:nd,1) = CKCenter(1:nd,1) / CKSize(1)  ! Now it is sum instead of mean
            CKCenter(1:nd,2) = CKCenter(1:nd,2) / CKSize(2)  ! Now it is sum instead of mean

            exit loopRecursiveRunMinVolKernel

        end do loopRecursiveRunMinVolKernel

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        furtherClusteringCheck: if (sum(CKBEVolNormed)<ClusterLogVolume(1) .or. ClusterLogVolume(1)>tightness*parVolNormed) then ! at least two clusters is better than one. now search for more child clusters

            CKnemax(1) = max( CKSize(1)/(ndEffective+1) , 1 )
            CKnemax(2) = max( CKSize(2)/(ndEffective+1) , 1 )

            IpStart(1) = 1             ; IpEnd(1) = CKSize(1)
            IpStart(2) = CKSize(1) + 1 ; IpEnd(2) = np

            loopSearchForChildren: do ic = 1,2   ! Search for more grand-child clusters

                if (ic==1) then
                    icStart = 1
                    icEnd = CKnemax(1)
                else
                    icStart = CKneopt(1) + 1
                    icEnd = CKneopt(1) + CKnemax(2)
                end if

                ClusterLogVolume(icStart) = CKBEVolNormed(ic)
                if (CKnemax(ic)>1) then ! .and. CKBEVolNormed(ic)>1.1_RK*CKParVolNormed(ic)) then
                    call runMinVolKernel( nemax = CKnemax(ic) & ! LCOV_EXCL_LINE
                                        , nd = nd & ! LCOV_EXCL_LINE
                                        , np = CKSize(ic) & ! LCOV_EXCL_LINE
                                        , stabilizationRequested = stabilizationRequested & ! LCOV_EXCL_LINE
                                        , maxKvolumeLoopRecursion = maxKvolumeLoopRecursion & ! LCOV_EXCL_LINE
                                        , maxAllowedKmeansFailure = maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                        , maxAllowedMinVolFailure = maxAllowedMinVolFailure & ! LCOV_EXCL_LINE
                                        , mahalSqWeightExponent = mahalSqWeightExponent & ! LCOV_EXCL_LINE
                                        , tightness = tightness & ! LCOV_EXCL_LINE
                                        , inclusionFraction = inclusionFraction & ! LCOV_EXCL_LINE
                                        , parVolNormed = CKParVolNormed(ic) & ! LCOV_EXCL_LINE
                                        , logNormFac = logNormFac & ! LCOV_EXCL_LINE
                                        , Point = Point(1:nd,IpStart(ic):IpEnd(ic)) & ! LCOV_EXCL_LINE
                                        , neopt = CKneopt(ic) & ! LCOV_EXCL_LINE
                                        , ClusterCenter = ClusterCenter(1:nd,icStart:icEnd) & ! LCOV_EXCL_LINE
                                        , ClusterSize = ClusterSize(icStart:icEnd) & ! LCOV_EXCL_LINE
                                        , ClusterInvCovMat = ClusterInvCovMat(1:nd,1:nd,icStart:icEnd) & ! LCOV_EXCL_LINE
                                        , ClusterChoLowCovUpp = ClusterChoLowCovUpp(1:nd,1:nd,icStart:icEnd) & ! LCOV_EXCL_LINE
                                        , ClusterChoDia = ClusterChoDia(1:nd,icStart:icEnd) & ! LCOV_EXCL_LINE
                                        , ClusterLogVolume = ClusterLogVolume(icStart:icEnd) & ! LCOV_EXCL_LINE
                                        , ClusterMembership = ClusterMembership(IpStart(ic):IpEnd(ic)) & ! LCOV_EXCL_LINE
                                        , PointIndex = PointIndex(IpStart(ic):IpEnd(ic)) & ! LCOV_EXCL_LINE
                                        , convergenceFailureCount = convergenceFailureCount & ! LCOV_EXCL_LINE
                                        )
                    if (ic==2) ClusterMembership(IpStart(ic):IpEnd(ic)) = ClusterMembership(IpStart(ic):IpEnd(ic)) + icStart - 1
                else
                    CKneopt(ic)=1
                end if

                if ( CKneopt(ic)==1 .or. CKnemax(ic)==1 ) then ! There is only one cluster here
                    CKneopt(ic) = 1
                    ClusterCenter(1:nd,icStart) = CKCenter(1:nd,ic)
                    ClusterSize(icStart) = CKSize(ic)
                    ClusterInvCovMat(1:nd,1:nd,icStart) = CKBEInvCovMat(1:nd,1:nd,ic) / ScaleFacSq(ic)
                    ClusterChoLowCovUpp(1:nd,1:nd,icStart) = CKBEChoLowCovUpp(1:nd,1:nd,ic)
                    ClusterChoDia(1:nd,icStart) = CKBEChoDia(1:nd,ic)
                    do j = 1,nd-1 ! Now rescale Cholesky factor and its diagonals, such that they describe the bounding ellipsoid of the cluster
                        ClusterChoDia(j,icStart) = ClusterChoDia(j,icStart) * scaleFactor(ic)
                        do i = j+1,nd
                            ClusterChoLowCovUpp(i,j,icStart) = ClusterChoLowCovUpp(i,j,icStart) * scaleFactor(ic)
                        end do
                    end do
                    ClusterChoDia(nd,icStart) = ClusterChoDia(nd,icStart) * scaleFactor(ic)
                end if

            end do loopSearchForChildren

            neopt = sum(CKneopt)

        else furtherClusteringCheck  ! one cluster is better

            neopt = 1
            ClusterMembership = 1
            ClusterSize(1) = np

        end if furtherClusteringCheck

        !write(*,*) 'nemax: ', nemax
        !write(*,*) 'neopt: ', neopt
        !write(*,*) 'ClusterSize: ', ClusterSize(1:neopt)
        !write(*,*) 'ClusterCenter: ', ClusterCenter(:,1:neopt)
        !write(*,*) 'Membership: ', ClusterMembership(:)
        !read(*,*)

    end subroutine runMinVolKernel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine writeMinVolPartitionOOP(self, fileUnit, nd, np, Point)
        use Constants_mod, only: IK, RK
        implicit none
        class(MinVolPartition_type) , intent(in)    :: self
        integer(IK)                 , intent(in)    :: fileUnit, nd, np
        real(RK)                    , intent(in)    :: Point(nd,np)
        call writeMinVolPartition   ( fileUnit = fileUnit & ! LCOV_EXCL_LINE
                                    , nd = nd & ! LCOV_EXCL_LINE
                                    , np = np & ! LCOV_EXCL_LINE
                                    , nc = self%neopt & ! LCOV_EXCL_LINE
                                    , ClusterSize = self%Size & ! LCOV_EXCL_LINE
                                    , ClusterCenter = self%Center & ! LCOV_EXCL_LINE
                                    , ClusterLogVolume = self%LogVolNormed & ! LCOV_EXCL_LINE
                                    , ClusterMembership = self%Membership & ! LCOV_EXCL_LINE
                                    , ClusterChoLowCovUpp = self%ChoLowCovUpp & ! LCOV_EXCL_LINE
                                    , ClusterChoDia = self%ChoDia & ! LCOV_EXCL_LINE
                                    , Point = Point & ! LCOV_EXCL_LINE
                                    )
    end subroutine writeMinVolPartitionOOP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine writeMinVolPartition ( fileUnit & ! LCOV_EXCL_LINE
                                    , nd,np,nc & ! LCOV_EXCL_LINE
                                    , ClusterSize & ! LCOV_EXCL_LINE
                                    , ClusterCenter & ! LCOV_EXCL_LINE
                                    , ClusterLogVolume & ! LCOV_EXCL_LINE
                                    , ClusterMembership & ! LCOV_EXCL_LINE
                                    , ClusterChoLowCovUpp & ! LCOV_EXCL_LINE
                                    , ClusterChoDia & ! LCOV_EXCL_LINE
                                    , Point & ! LCOV_EXCL_LINE
                                    )
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK)     , intent(in)    :: fileUnit
        integer(IK)     , intent(in)    :: nd, np, nc
        integer(IK)     , intent(in)    :: ClusterSize(nc)
        integer(IK)     , intent(in)    :: ClusterMembership(np)
        real(RK)        , intent(in)    :: ClusterLogVolume(nc)
        real(RK)        , intent(in)    :: ClusterCenter(nd,nc)
        real(RK)        , intent(in)    :: ClusterChoLowCovUpp(nd,nd,nc)
        real(RK)        , intent(in)    :: ClusterChoDia(nd,nc)
        real(RK)        , intent(in)    :: Point(nd,np)
        character(*)    , parameter     :: fileFormat = "(*(g0.8,:,','))"
        integer(IK)                     :: i, j, ic

        write(fileUnit,"(A)") "nd, np, nc"
        write(fileUnit,fileFormat) nd, np, nc

        write(fileUnit,"(A)") "Size"
        write(fileUnit,fileFormat) ClusterSize

        write(fileUnit,"(A)") "Center"
        write(fileUnit,fileFormat) ClusterCenter

        write(fileUnit,"(A)") "LogVolume"
        write(fileUnit,fileFormat) ClusterLogVolume

        write(fileUnit,"(A)") "CholeskyLower"
        write(fileUnit,fileFormat) ((ClusterChoDia(j,ic), (ClusterChoLowCovUpp(i,j,ic), i=j+1,nd), j=1,nd), ic=1,nc)

!write(*,fileFormat)
!write(*,fileFormat)
!write(*,fileFormat) ((ClusterChoDia(j,ic), (ClusterChoLowCovUpp(i,j,ic), i=j+1,nd), j=1,nd), ic=1,nc)
!write(*,fileFormat) ((ClusterChoDia(j,ic), j=1,nd), ic=1,nc)
!write(*,fileFormat) (((ClusterChoLowCovUpp(i,j,ic), i=j+1,nd), j=1,nd), ic=1,nc)
!write(*,fileFormat)
!write(*,fileFormat)
!write(*,fileFormat) (((ClusterChoLowCovUpp(i,j,ic), i=1,j), j=1,nd), ic=1,nc)
block
real(RK), allocatable :: cholow(:,:), covmat(:,:)
allocate(cholow(nd,nd), source = 0._RK)
do j = 1, nd
    cholow(j,j) = ClusterChoDia(j,1)
    do i = j+1,nd
        cholow(i,j) = ClusterChoLowCovUpp(i,j,1)
    end do
end do
covmat = matmul(cholow,transpose(cholow))
write(*,fileFormat) ((covmat(i,j), i=1,j), j=1,nd)
end block

        write(fileUnit,"(A)") "Point"
        write(fileUnit,fileFormat) Point

        write(fileUnit,"(A)") "Membership"
        write(fileUnit,fileFormat) ClusterMembership

    end subroutine writeMinVolPartition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  end module MinVolPartition_mod