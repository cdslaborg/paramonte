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
!> This file implements the bodies of the different Partitioning modules that are based on different optimization criteria.
!> \author
!> Amir Shahmoradi

#define OPTIMIZATION_ENABLED

    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type

    implicit none

#if defined MAXDEN
    character(*), parameter :: MODULE_NAME = "@PartitionMaxDen_mod"
#elif defined MINVOL
    character(*), parameter :: MODULE_NAME = "@PartitionMinVol_mod"
#else
    character(*), parameter :: MODULE_NAME = "@Partition_mod"
#endif

    !> The `Partition_type` class. Partitions an input array `Point(nd,np)` with `nd` attributes and `np` observations (points).
    type :: Partition_type
        integer(IK)                 :: nd,np                            !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        integer(IK)                 :: nemax                            !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        integer(IK)                 :: neopt                            !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        integer(IK)                 :: numRecursiveCall                 !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        integer(IK)                 :: convergenceFailureCount          !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        integer(IK)                 :: partitionMaxAllowedFailure       !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        integer(IK)                 :: partitionMaxAllowedRecursion     !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        integer(IK)                 :: partitionMaxAllowedKmeansFailure !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        logical                     :: partitionStabilizationRequested  !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        logical                     :: partitionOptimizationRequested   !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
#if !defined MAXDEN && !defined MINVOL
        real(RK)                    :: mahalSqWeightExponent            !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
#endif
        real(RK)                    :: inclusionFraction                !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        real(RK)                    :: logTightness                     !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        real(RK)                    :: parLogVol                        !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        integer(IK) , allocatable   :: Membership(:)                    !< An array of size `(np)` representing the bounding-ellipsoid membership IDs of the corresponding data points.
        integer(IK) , allocatable   :: PointIndex(:)                    !< An array of size `(np)` representing the indices of the original data points before partitioning.
        real(RK)    , allocatable   :: LogVolNormed(:)                  !< An array of size `(nemax)` representing the log-volumes of the corresponding bounding ellipsoids.
        real(RK)    , allocatable   :: ChoLowCovUpp(:,:,:)              !< An array of size `(nd,nd,nemax)` representing the Cholesky lower triangle, diagonal, and covariance matrices of the bounding ellipsoids.
        real(RK)    , allocatable   :: InvCovMat(:,:,:)                 !< An array of size `(nd,nd,nemax)` representing the full symmetric inverse covariance matrices of the bounding ellipsoids.
        real(RK)    , allocatable   :: ChoDia(:,:)                      !< An array of size `(nd,nemax)` representing the Cholesky lower triangle, diagonal, and covariance matrices of the bounding ellipsoids.
        real(RK)    , allocatable   :: Center(:,:)                      !< An array of size `(nd,nemax)` representing the centers of the bounding ellipsoids.
        integer(IK) , allocatable   :: Size(:)                          !< An array of size `(nemax)` representing the sizes of the corresponding bounding ellipsoids.
        type(Err_type)              :: Err                              !< An object of class [Err_type](@ref err_mod::err_type) containing information about error occurrence.
    contains
        procedure, pass :: write => writePartitionOOP
    end type Partition_type

    interface Partition_type
        module procedure :: constructPartition
    end interface Partition_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructPartition ( Point, nd, np, nemax & ! LCOV_EXCL_LINE
                                , partitionMaxAllowedFailure & ! LCOV_EXCL_LINE
                                , partitionMaxAllowedRecursion & ! LCOV_EXCL_LINE
                                , partitionMaxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                , partitionStabilizationRequested & ! LCOV_EXCL_LINE
                                , partitionOptimizationRequested & ! LCOV_EXCL_LINE
#if !defined MAXDEN && !defined MINVOL
                                , mahalSqWeightExponent & ! LCOV_EXCL_LINE
#endif
                                , inclusionFraction & ! LCOV_EXCL_LINE
                                , logTightness & ! LCOV_EXCL_LINE
                                , trimEnabled & ! LCOV_EXCL_LINE
                                , parLogVol & ! LCOV_EXCL_LINE
                                ) result(self)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructPartition
#endif
        use Constants_mod, only: IK, RK !, NEGINF_RK

        implicit none

        integer(IK) , intent(in)                :: nd,np
        integer(IK) , intent(in), optional      :: nemax
        real(RK)    , intent(inout)             :: Point(nd,np) ! do not move this up. This must appear after nd, np declaration.
        integer(IK) , intent(in), optional      :: partitionMaxAllowedFailure
        integer(IK) , intent(in), optional      :: partitionMaxAllowedRecursion
        integer(IK) , intent(in), optional      :: partitionMaxAllowedKmeansFailure
        logical     , intent(in), optional      :: partitionStabilizationRequested
        logical     , intent(in), optional      :: partitionOptimizationRequested
#if !defined MAXDEN && !defined MINVOL
        real(RK)    , intent(in), optional      :: mahalSqWeightExponent
#endif
        real(RK)    , intent(in), optional      :: inclusionFraction
        real(RK)    , intent(in), optional      :: logTightness
        logical     , intent(in), optional      :: trimEnabled
        real(RK)    , intent(inout), optional   :: parLogVol
        type(Partition_type)                    :: self

        logical                                 :: trimEnabledDefault

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

        self%partitionMaxAllowedFailure = 10_IK; if (present(partitionMaxAllowedFailure)) self%partitionMaxAllowedFailure = partitionMaxAllowedFailure
        self%partitionMaxAllowedRecursion = 50_IK; if (present(partitionMaxAllowedRecursion)) self%partitionMaxAllowedRecursion = partitionMaxAllowedRecursion
        self%partitionMaxAllowedKmeansFailure = 10_IK; if (present(partitionMaxAllowedKmeansFailure)) self%partitionMaxAllowedKmeansFailure = partitionMaxAllowedKmeansFailure
        self%partitionStabilizationRequested = .true.; if (present(partitionStabilizationRequested)) self%partitionStabilizationRequested = partitionStabilizationRequested
        self%partitionOptimizationRequested = .true.; if (present(partitionOptimizationRequested)) self%partitionOptimizationRequested = partitionOptimizationRequested
#if !defined MAXDEN && !defined MINVOL
        self%mahalSqWeightExponent = 0._RK; if (present(mahalSqWeightExponent)) self%mahalSqWeightExponent = mahalSqWeightExponent
#endif
        self%inclusionFraction = 0._RK; if (present(inclusionFraction)) self%inclusionFraction = inclusionFraction
        self%logTightness = log(1.3_RK); if (present(logTightness)) self%logTightness = logTightness
        self%parLogVol = -huge(1._RK); if (present(parLogVol)) self%parLogVol = parLogVol

        allocate( self%Size         (self%nemax) & ! LCOV_EXCL_LINE
                , self%Center       (nd,self%nemax) & ! LCOV_EXCL_LINE
                , self%ChoDia       (nd,self%nemax) & ! LCOV_EXCL_LINE
                , self%InvCovMat    (nd,nd,self%nemax) & ! LCOV_EXCL_LINE
                , self%ChoLowCovUpp (nd,nd,self%nemax) & ! LCOV_EXCL_LINE
                , self%LogVolNormed (self%nemax) & ! LCOV_EXCL_LINE
                , self%Membership   (np) & ! LCOV_EXCL_LINE
                , self%PointIndex   (np) & ! LCOV_EXCL_LINE
                )

        call runPartition   ( Point = Point & ! LCOV_EXCL_LINE
                            , nd = nd & ! LCOV_EXCL_LINE
                            , np = np & ! LCOV_EXCL_LINE
                            , nemax = self%nemax & ! LCOV_EXCL_LINE
                            , partitionMaxAllowedFailure = self%partitionMaxAllowedFailure & ! LCOV_EXCL_LINE
                            , partitionMaxAllowedRecursion = self%partitionMaxAllowedRecursion & ! LCOV_EXCL_LINE
                            , partitionMaxAllowedKmeansFailure = self%partitionMaxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                            , partitionStabilizationRequested = self%partitionStabilizationRequested & ! LCOV_EXCL_LINE
                            , partitionOptimizationRequested = self%partitionOptimizationRequested & ! LCOV_EXCL_LINE
#if !defined MAXDEN && !defined MINVOL
                            , mahalSqWeightExponent = self%mahalSqWeightExponent & ! LCOV_EXCL_LINE
#endif
                            , inclusionFraction = self%inclusionFraction & ! LCOV_EXCL_LINE
                            , logTightness = self%logTightness & ! LCOV_EXCL_LINE
                            , neopt = self%neopt & ! LCOV_EXCL_LINE
                            , PointIndex = self%PointIndex & ! LCOV_EXCL_LINE
                            , PartitionSize = self%Size & ! LCOV_EXCL_LINE
                            , PartitionCenter = self%Center & ! LCOV_EXCL_LINE
                            , PartitionChoDia = self%ChoDia & ! LCOV_EXCL_LINE
                            , PartitionLogVol = self%LogVolNormed & ! LCOV_EXCL_LINE
                            , PartitionInvCovMat = self%InvCovMat & ! LCOV_EXCL_LINE
                            , PartitionMembership = self%Membership & ! LCOV_EXCL_LINE
                            , PartitionChoLowCovUpp = self%ChoLowCovUpp & ! LCOV_EXCL_LINE
                            , convergenceFailureCount = self%convergenceFailureCount & ! LCOV_EXCL_LINE
                            , numRecursiveCall = self%numRecursiveCall & ! LCOV_EXCL_LINE
                            , Err = self%Err & ! LCOV_EXCL_LINE
                            , parLogVol = parLogVol & ! LCOV_EXCL_LINE
                            )

        if (trimEnabledDefault) then
            self%Size           = self%Size         (1:self%neopt)
            self%Center         = self%Center       (1:nd,1:self%neopt)
            self%ChoDia         = self%ChoDia       (1:nd,1:self%neopt)
            self%InvCovMat      = self%InvCovMat    (1:nd,1:nd,1:self%neopt)
            self%ChoLowCovUpp   = self%ChoLowCovUpp (1:nd,1:nd,1:self%neopt)
            self%LogVolNormed   = self%LogVolNormed (1:self%neopt)
        end if

    end function constructPartition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine runPartition ( Point & ! LCOV_EXCL_LINE
                            , nd,np & ! LCOV_EXCL_LINE
                            , nemax & ! LCOV_EXCL_LINE
                            , partitionMaxAllowedRecursion & ! LCOV_EXCL_LINE
                            , partitionMaxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                            , partitionMaxAllowedFailure & ! LCOV_EXCL_LINE
                            , partitionStabilizationRequested & ! LCOV_EXCL_LINE
                            , partitionOptimizationRequested & ! LCOV_EXCL_LINE
#if !defined MAXDEN && !defined MINVOL
                            , mahalSqWeightExponent & ! LCOV_EXCL_LINE
#endif
                            , inclusionFraction & ! LCOV_EXCL_LINE
                            , logTightness & ! LCOV_EXCL_LINE
                            ! output
                            , neopt & ! LCOV_EXCL_LINE
                            , PointIndex & ! LCOV_EXCL_LINE
                            , PartitionSize & ! LCOV_EXCL_LINE
                            , PartitionCenter & ! LCOV_EXCL_LINE
                            , PartitionChoDia & ! LCOV_EXCL_LINE
                            , PartitionLogVol & ! LCOV_EXCL_LINE
                            , PartitionInvCovMat & ! LCOV_EXCL_LINE
                            , PartitionMembership & ! LCOV_EXCL_LINE
                            , PartitionChoLowCovUpp & ! LCOV_EXCL_LINE
                            , convergenceFailureCount & ! LCOV_EXCL_LINE
                            , numRecursiveCall & ! LCOV_EXCL_LINE
                            , Err & ! LCOV_EXCL_LINE
                            , parLogVol & ! LCOV_EXCL_LINE
                            )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: runPartition
#endif
        use Constants_mod, only: IK, RK, NEGINF_RK
        use Matrix_mod, only: getInvMatFromCholFac
        use Matrix_mod, only: getCholeskyFactor
        use String_mod, only: num2str
        use Err_mod, only: Err_type

        implicit none

        integer(IK)     , intent(in)                :: nd,np                                ! the maximum allowed number of ellipsoidal clusters, number of dimensions, number of points
        integer(IK)     , intent(in)                :: nemax                                ! the maximum allowed number of ellipsoidal clusters, number of dimensions, number of points
        real(RK)        , intent(inout)             :: Point(nd,np)                         ! the input data (points) to be clustered. on output, it is ordered. This must appear after nd, np declaration.
        integer(IK)     , intent(in)                :: partitionMaxAllowedFailure
        integer(IK)     , intent(in)                :: partitionMaxAllowedRecursion
        integer(IK)     , intent(in)                :: partitionMaxAllowedKmeansFailure
        logical         , intent(in)                :: partitionStabilizationRequested
        logical         , intent(in)                :: partitionOptimizationRequested
#if !defined MAXDEN && !defined MINVOL
        real(RK)        , intent(in)                :: mahalSqWeightExponent
#endif
        real(RK)        , intent(in)                :: inclusionFraction                    ! the fraction of non-cluster-member points that ARE inside the cluster, to be considered in the estimation of the true volume of the cluster.
        real(RK)        , intent(in)                :: logTightness                         ! the factor deremining how large the bounding ellipsoid can be, compared to the true volume of the ellipsoid.
        integer(IK)     , intent(out)               :: neopt                                ! this variable is saved throughout recursive iterations of subroutine
        integer(IK)     , intent(out)               :: PointIndex(np)                       ! Vector containing the scrambled indices of Point elements. It represent the mapping from the input Point elements to output scrambled Point vector.
        integer(IK)     , intent(out)               :: PartitionSize(nemax)                 ! cluster sizes (nemax)
        real(RK)        , intent(out)               :: PartitionLogVol(nemax)               ! CBE: Cluster Bounding Ellipsoid - On output: Each element of the array contains the volume of the corresponding Cluster Bounding Ellipsoid.
        integer(IK)     , intent(out)               :: PartitionMembership(np)              ! the PartitionMembership id of each point (representing the cluster number to which the point belongs)
        real(RK)        , intent(out)               :: PartitionCenter(nd,nemax)            ! Centers of the clusters (nd,nemax)
        real(RK)        , intent(out)               :: PartitionChoDia(nd,nemax)            ! CBE: Cluster Bounding Ellipsoid - On output: The diagonal elements of the Cholesky factor of the covariance matrix of the input points, on output, contains the cholesky factor of all found clusters.
        real(RK)        , intent(out)               :: PartitionInvCovMat(nd,nd,nemax)      ! CBE: Cluster Bounding Ellipsoid - On output: Full symmetric Inverse Covariance Matrices of the Bounding Ellipsoids of all clusters.
        real(RK)        , intent(out)               :: PartitionChoLowCovUpp(nd,nd,nemax)   ! CBE: Cluster Bounding Ellipsoid - On output: The lower triangle contains the Cholesky factor of the covariance matrix of the input points, the upper triangle maybe the CovMat itself (not needed). On output contains in the lower triangle, the Cholesky factor of all clusters found.
        integer(IK)     , intent(out)               :: convergenceFailureCount              ! The number of times the algorithm failed to converge.
        integer(IK)     , intent(out)               :: numRecursiveCall                     ! The number of recursive calls before convergence occurs.
        type(Err_type)  , intent(out)               :: Err                                  ! Error object containing error information.
        real(RK)        , intent(inout), optional   :: parLogVol                            ! input log-estimate of the total volume of the points. Input zero or negative for the value to be ignored.

        logical                                     :: parLogVolIsPresent
        logical                                     :: boundedRegionIsTooLarge
        real(RK)                                    :: scaleFactorSqInverse
        real(RK)                                    :: NormedPoint(nd,np)
        real(RK)                                    :: scaleFactorSq
        real(RK)                                    :: scaleFactor
        real(RK)                                    :: mahalSq
        integer                                     :: i, j, ip

        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME//"@runPartition()"

        parLogVolIsPresent = present(parLogVol)

        Err%occurred = .false.

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Initialize the first cluster's center
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        PartitionCenter(1:nd,1) = 0._RK
        do ip = 1, np
          PointIndex(ip) = ip
          PartitionCenter(1:nd,1) = PartitionCenter(1:nd,1) + Point(1:nd,ip)
        end do
        PartitionCenter(1:nd,1) = PartitionCenter(1:nd,1) / real(np,kind=RK)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the covariance matrix and the Cholesky factorization of the first cluster.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Compute the covariance matrix upper.

        PartitionChoLowCovUpp(1:nd,1:nd,1) = 0._RK
        do ip = 1, np
            do j = 1, nd
                NormedPoint(j,ip) = Point(j,ip) - PartitionCenter(j,1)
                PartitionChoLowCovUpp(1:j,j,1) = PartitionChoLowCovUpp(1:j,j,1) + NormedPoint(1:j,ip) * NormedPoint(j,ip)
            end do
        end do
        PartitionChoLowCovUpp(1:nd,1:nd,1) = PartitionChoLowCovUpp(1:nd,1:nd,1) / real(np-1,kind=RK)

        ! Compute the Cholesky factor lower.

        call getCholeskyFactor(nd, PartitionChoLowCovUpp(1:nd,1:nd,1), PartitionChoDia(1:nd,1))
        if (PartitionChoDia(1,1)<0._RK) then
            Err%msg = PROCEDURE_NAME//": Singular covariance matrix detected for nemax, nd, np: "//num2str(nemax)//", "//num2str(nd)//", "//num2str(np)
            write(*,*) Err%msg
            error stop
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the inverse covariance matrix.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        PartitionInvCovMat(1:nd,1:nd,1) = getInvMatFromCholFac(nd, PartitionChoLowCovUpp(1:nd,1:nd,1), PartitionChoDia(1:nd,1))

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the Mahal distances squared and the scale factor which makes the covariance matrix a bounding ellipsoid.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        scaleFactorSq = NEGINF_RK
        do ip = 1, np
            mahalSq = dot_product( NormedPoint(1:nd,ip) , matmul(PartitionInvCovMat(1:nd,1:nd,1), NormedPoint(1:nd,ip)) )
            if (scaleFactorSq < mahalSq) scaleFactorSq = mahalSq
        end do
        scaleFactor = sqrt(scaleFactorSq)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the largest possible LogVolume of the first cluster based on the parent cluster and the input estimated volume.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Compute the volume of the bounded cluster.

        PartitionLogVol(1) = sum( log( PartitionChoDia(1:nd,1) ) ) + nd * log(scaleFactor)

        ! Rescale the scaleFactor, if needed, to enlarge the bounded cluster in the future.
        ! This must be done here, although it will be used at the very end.
        ! Otherwise, PartitionLogVol(1) will be potentially overwritten.

        if (parLogVolIsPresent) then
            boundedRegionIsTooLarge = PartitionLogVol(1) > parLogVol
            if (boundedRegionIsTooLarge) then
                boundedRegionIsTooLarge = PartitionLogVol(1) > logTightness + parLogVol
            else
                scaleFactor = scaleFactor * exp( (parLogVol - PartitionLogVol(1)) / nd )
                scaleFactorSq = scaleFactor**2
                PartitionLogVol(1) = parLogVol
            end if
        else
            boundedRegionIsTooLarge = .true.
        end if

        if (nemax > 1 .and. boundedRegionIsTooLarge) then
            numRecursiveCall = 0_IK
            convergenceFailureCount = 0_IK
            call runPartitionKernel ( Point = Point & ! LCOV_EXCL_LINE
                                    , nd = nd & ! LCOV_EXCL_LINE
                                    , np = np & ! LCOV_EXCL_LINE
                                    , nemax = nemax & ! LCOV_EXCL_LINE
                                    , partitionMaxAllowedFailure = partitionMaxAllowedFailure & ! LCOV_EXCL_LINE
                                    , partitionMaxAllowedRecursion = partitionMaxAllowedRecursion & ! LCOV_EXCL_LINE
                                    , partitionMaxAllowedKmeansFailure = partitionMaxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                    , partitionStabilizationRequested = partitionStabilizationRequested & ! LCOV_EXCL_LINE
                                    , partitionOptimizationRequested = partitionOptimizationRequested & ! LCOV_EXCL_LINE
#if !defined MAXDEN && !defined MINVOL
                                    , mahalSqWeightExponent = mahalSqWeightExponent & ! LCOV_EXCL_LINE
#endif
                                    , inclusionFraction = inclusionFraction & ! LCOV_EXCL_LINE
                                    , logTightness = logTightness & ! LCOV_EXCL_LINE
                                    , neopt = neopt & ! LCOV_EXCL_LINE
                                    , PointIndex = PointIndex & ! LCOV_EXCL_LINE
                                    , PartitionSize = PartitionSize & ! LCOV_EXCL_LINE
                                    , PartitionLogVol = PartitionLogVol & ! LCOV_EXCL_LINE
                                    , PartitionCenter = PartitionCenter & ! LCOV_EXCL_LINE
                                    , PartitionChoDia = PartitionChoDia & ! LCOV_EXCL_LINE
                                    , PartitionInvCovMat = PartitionInvCovMat & ! LCOV_EXCL_LINE
                                    , PartitionMembership = PartitionMembership & ! LCOV_EXCL_LINE
                                    , PartitionChoLowCovUpp = PartitionChoLowCovUpp & ! LCOV_EXCL_LINE
                                    , convergenceFailureCount = convergenceFailureCount & ! LCOV_EXCL_LINE
                                    , numRecursiveCall = numRecursiveCall & ! LCOV_EXCL_LINE
                                    , parLogVol = parLogVol & ! LCOV_EXCL_LINE
                                    )
        else ! There is only one cluster
            neopt = 1_IK
            PartitionSize(1) = np ! @todo : xxx this is redundant and should be removed in the future.
            PartitionMembership(1:np) = 1_IK ! @todo : xxx this is redundant and should be removed in the future.
        end if

        if (neopt==1_IK) then

            ! Now rescale Cholesky factor and its diagonals, such that they describe the bounding ellipsoid of the cluster.

            scaleFactorSqInverse = 1._RK / scaleFactorSq
            do j = 1, nd
                PartitionChoDia(j,1) = PartitionChoDia(j,1) * scaleFactor
                !PartitionChoLowCovUpp(1:j,j,1) = PartitionChoLowCovUpp(1:j,j,1) * scaleFactorSq
                PartitionInvCovMat(1:nd,j,1) = PartitionInvCovMat(1:nd,j,1) * scaleFactorSqInverse
                do i = j+1, nd
                    PartitionChoLowCovUpp(i,j,1) = PartitionChoLowCovUpp(i,j,1) * scaleFactor
                end do
            end do

        end if

    end subroutine runPartition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This code assumes that nemax > 1 on input.
    ! important input requirements:
    recursive subroutine runPartitionKernel ( Point & ! LCOV_EXCL_LINE
                                            , nd & ! LCOV_EXCL_LINE
                                            , np & ! LCOV_EXCL_LINE
                                            , nemax & ! LCOV_EXCL_LINE
                                            , partitionMaxAllowedFailure & ! LCOV_EXCL_LINE
                                            , partitionMaxAllowedRecursion & ! LCOV_EXCL_LINE
                                            , partitionMaxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                            , partitionStabilizationRequested & ! LCOV_EXCL_LINE
                                            , partitionOptimizationRequested & ! LCOV_EXCL_LINE
#if !defined MAXDEN && !defined MINVOL
                                            , mahalSqWeightExponent & ! LCOV_EXCL_LINE
#endif
                                            , inclusionFraction & ! LCOV_EXCL_LINE
                                            , logTightness & ! LCOV_EXCL_LINE
                                            , neopt & ! LCOV_EXCL_LINE
                                            , PointIndex & ! LCOV_EXCL_LINE
                                            , PartitionSize & ! LCOV_EXCL_LINE
                                            , PartitionCenter & ! LCOV_EXCL_LINE
                                            , PartitionChoDia & ! LCOV_EXCL_LINE
                                            , PartitionLogVol & ! LCOV_EXCL_LINE
                                            , PartitionInvCovMat & ! LCOV_EXCL_LINE
                                            , PartitionMembership & ! LCOV_EXCL_LINE
                                            , PartitionChoLowCovUpp & ! LCOV_EXCL_LINE
                                            , convergenceFailureCount & ! LCOV_EXCL_LINE
                                            , numRecursiveCall & ! LCOV_EXCL_LINE
                                            , parLogVol & ! LCOV_EXCL_LINE
                                            )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: runPartitionKernel
#endif
        use Constants_mod, only: IK, RK, SPR, HUGE_RK
        use Matrix_mod, only: getInvMatFromCholFac
        use Matrix_mod, only: getCholeskyFactor
        use Kmeans_mod, only: runKmeans
        use Math_mod, only: getLogSumExp
        use Err_mod, only: Err_type

        implicit none

        integer(IK) , intent(in)                :: nemax                            ! maximum number of clusters determined from the nd and np
        integer(IK) , intent(in)                :: nd,np                            ! number of dimensions, number of points
        integer(IK) , intent(in)                :: partitionMaxAllowedFailure
        integer(IK) , intent(in)                :: partitionMaxAllowedRecursion
        integer(IK) , intent(in)                :: partitionMaxAllowedKmeansFailure
        logical     , intent(in)                :: partitionStabilizationRequested
        logical     , intent(in)                :: partitionOptimizationRequested
#if !defined MAXDEN && !defined MINVOL
        real(RK)    , intent(in)                :: mahalSqWeightExponent
#endif
        real(RK)    , intent(in)                :: inclusionFraction                    ! the fraction of non-cluster-member points that ARE inside the cluster, to be considered in the estimation of the true volume of the cluster.
        real(RK)    , intent(in)                :: logTightness                         ! the factor deremining how large the bounding ellipsoid can be, compared to the true volume of the ellipsoid.
        integer(IK) , intent(out)               :: neopt                                ! The predicted optimal number of clusters identified in the input data.
        integer(IK) , intent(out)               :: PartitionSize(nemax)                 ! cluster sizes (nemax)
        integer(IK) , intent(inout)             :: convergenceFailureCount              ! returns the number of times the algorithm failed to converge.
        integer(IK) , intent(inout)             :: PartitionMembership(np)              ! the PartitionMembership id of each point (representing the cluster number to which the point belongs)
        integer(IK) , intent(inout)             :: PointIndex(np)                       ! Vector containing the scrambled indices of Point elements. It represent the mapping from the input Point elements to output scrambled Point vector.
        real(RK)    , intent(inout)             :: Point(nd,np)                         ! the input data (points) to be clustered.
        real(RK)    , intent(inout)             :: PartitionCenter(nd,nemax)            ! Centers of the clusters (nd,nemax)
        real(RK)    , intent(inout)             :: PartitionChoLowCovUpp(nd,nd,nemax)   ! CBE: Cluster Bounding Ellipsoid - On input: The lower triangle contains the Cholesky factor of the covariance matrix of the input points, the upper triangle maybe the CovMat itself (not needed). On output contains in the lower triangle, the Cholesky factor of all clusters found.
        real(RK)    , intent(inout)             :: PartitionInvCovMat(nd,nd,nemax)      ! CBE: Cluster Bounding Ellipsoid - On input: The lower triangle contains the Cholesky factor of the covariance matrix of the input points, the upper triangle maybe the CovMat itself (not needed). On output contains in the lower triangle, the Cholesky factor of all clusters found.
        real(RK)    , intent(inout)             :: PartitionChoDia(nd,nemax)            ! CBE: Cluster Bounding Ellipsoid - On input: The diagonal elements of the Cholesky factor of the covariance matrix of the input points, on output, contains the cholesky factor of all found clusters.
        real(RK)    , intent(inout)             :: PartitionLogVol(nemax)               ! CBE: Cluster Bounding Ellipsoid - On input: The volume of the Bounding Ellipsoid of the parent cluster. On output: The Bounding volume of all child clusters.
        integer(IK) , intent(inout)             :: numRecursiveCall                     ! The number of recursive calls before convergence occurs.
        real(RK)    , intent(inout), optional   :: parLogVol                            ! parent's sqrt(det(invCovMat)), where invCovMat is the inverse covariance matrix of the bounding ellipsoid of the input points.

        real(RK)        :: KmeansChoLowCovUpp(nd,nd,2)  , KmeansChoLowCovUppOld(nd,nd,2)    ! Covariance Matrices of the bounding ellipsoid of the Kmeans' 2 clusters, reported in the upper triangle. The lower triangle contains the Cholesky factors
        real(RK)        :: KmeansInvCovMat(nd,nd,2)     , KmeansInvCovMatOld(nd,nd,2)       ! Inverse of the Covariance Matrices of the bounding ellipsoid of the Kmeans' 2 clusters.
        real(RK)        :: KmeansScaleFactorSq(2)       , KmeansScaleFactorSqOld(2)         ! The factor by which determinant must be multiplied in order for the ellipsoid to become a bounding ellipsoid.
        real(RK)        :: KmeansScaleFactor(2)         , KmeansScaleFactorOld(2)           ! The factor by which determinant must be multiplied in order for the ellipsoid to become a bounding ellipsoid.
        real(RK)        :: KmeansMahalSq(np,2)          , KmeansMahalSqOld(np,2)            ! Vector of Maximum Mahalanobis distances for each cluster
        real(RK)        :: KmeansCenter(nd,2)           , KmeansCenterOld(nd,2)             ! Cluster Centers for nc=2 Kmeans algorithm
        real(RK)        :: KmeansChoDia(nd,2)           , KmeansChoDiaOld(nd,2)             ! determinants of the inverse covariance matrices of the bounding ellipsoids of the Kmeans' 2 clusters.
        real(RK)        :: KmeansLogVol(2)              , KmeansLogVolOld(2)                ! New Cluster Bounding Ellisoid Volumes
        real(RK)        :: KmeansPoint(nd,np)           ! Cluster Points in 2-means partition result
        real(RK)        :: NormedPoint(nd,np,2)         ! Vector of Maximum Mahalanobis distances for each cluster
        real(RK)        :: KmeansNormedPoint(nd,np,2)   ! Vector of Maximum Mahalanobis distances for each cluster
        real(RK)        :: KmeansLogVolEstimate(2)      ! New Cluster Volume estimates
        real(RK)        :: KmeansMinDistanceSq(np)      ! Vector of Maximum Mahalanobis distances for each cluster
        integer(IK)     :: KmeansNumPointInside(2)      ! cluster sizes for nc=2 Kmeans algorithm
        integer(IK)     :: KmeansMemberCounter(2)       ! cluster member counter dummy variable.
        integer(IK)     :: KmeansPointIndex(np)         ! Vector containing the scrambled indices of Point elements. It represent the mapping from the input Point elements to output scrambled Point vector.
        integer(IK)     :: KmeansNemax(2)               ! maximum allowed number of clusters in each of the two K-means clusters.
        integer(IK)     :: KmeansNeopt(2)               ! optimal number of clusters in each of the child clusters
        integer(IK)     :: KmeansSize(2)                ! cluster sizes for nc=2 Kmeans algorithm
        integer(IK)     :: zcsIteration                 ! Counts the number of times zero sized clusters are found, zcs stands for zero cluster size
        integer(IK)     :: IpStart(2)
        integer(IK)     :: IpEnd(2)
#if !defined MINVOL
        real(RK)        :: MahalSqWeight(2)             ! MahalSqWeight definition of Lu et al 2007
#endif

        integer(IK)     :: i, j, ip, ic, icOther, iteration, ndEffectivePlusOne
        integer(IK)     :: icStart,icEnd,recursionCounter,kmeansFailureCount,partitionFailureCount
        integer(IK)     :: minClusterSize
        logical         :: reclusteringNeeded
        logical         :: parLogVolIsPresent
        logical         :: boundedRegionIsTooLarge
        logical         :: isZeroInclusionFraction
        real(RK)        :: scaleFactorSqInverse
        real(RK)        :: potential                    ! The factor by which determinant must be multiplied in order for the ellipsoid to become a bounding ellipsoid.
        type(Err_type)  :: Err

        numRecursiveCall = numRecursiveCall + 1
        minClusterSize = nd + 1_IK
        ndEffectivePlusOne = 1_IK + np / (nemax + 1_IK)
        isZeroInclusionFraction = abs(inclusionFraction) < 1.e-6_RK

        parLogVolIsPresent = present(parLogVol)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Perform the initial Kmeans clustering.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (kmeansClusteringFailed()) return

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the scaled Mahalanobis distances and the bounded regions properties.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call computeBoundedResgion()

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Find two approximately minimum-volume sub-clusters and if requested, store a copy of the results for recovery.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined OPTIMIZATION_ENABLED

        partitionFailureCount = 0
        recursionCounter = 0

        loopRecursiveRunPartitionKernel: do

            ! If the maximum allowed number of recursion has reached, quit

            recursionCounter = recursionCounter + 1_IK
            if(recursionCounter > partitionMaxAllowedRecursion) then
                write(*,*) "recursionCounter > partitionMaxAllowedRecursion: ", recursionCounter
                neopt = 1_IK
                PartitionSize(1) = np
                PartitionMembership = 1_IK
                convergenceFailureCount = convergenceFailureCount + 1_IK
                return
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Ensure the stability of the volume shrinkage, if requested, by caching the previous estimates.
            ! This stabilization is mathematically sound and safe, as no new alien point is being added to the bounded clusters.
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            blockStabilityCheck: if (partitionOptimizationRequested .and. partitionStabilizationRequested) then

                if ( KmeansLogVol(ic) > KmeansLogVolOld(ic) ) then ! This must never occur on the first try

                    KmeansChoLowCovUpp(1:nd,1:nd,ic) = KmeansChoLowCovUppOld(1:nd,1:nd,ic)
                    KmeansInvCovMat(1:nd,1:nd,ic) = KmeansInvCovMatOld(1:nd,1:nd,ic)
                    KmeansScaleFactorSq(ic) = KmeansScaleFactorSqOld(ic)
                    KmeansScaleFactor(ic) = KmeansScaleFactorOld(ic)
                    KmeansCenter(1:nd,ic) = KmeansCenterOld(1:nd,ic)
                    KmeansChoDia(1:nd,ic) = KmeansChoDiaOld(1:nd,ic)
                    KmeansLogVol(ic) = KmeansLogVolOld(ic)

                    ! MahalSq must be handled separately, since the order of the points has changed.

                    KmeansMemberCounter(1) = 0
                    KmeansMemberCounter(2) = KmeansSize(1)
                    do ip = 1, np
                        KmeansMemberCounter(PartitionMembership(ip)) = KmeansMemberCounter(PartitionMembership(ip)) + 1_IK
                        KmeansMahalSq(KmeansMemberCounter(PartitionMembership(ip)),ic) = KmeansMahalSqOld(ip,ic)
                    end do
                    KmeansMahalSqOld(1:np,ic) = KmeansMahalSq(1:np,ic)

                else  ! no stability correction is required. Everything just fine...

                    KmeansChoLowCovUppOld(1:nd,1:nd,ic) = KmeansChoLowCovUpp(1:nd,1:nd,ic)
                    KmeansInvCovMatOld(1:nd,1:nd,ic)    = KmeansInvCovMat(1:nd,1:nd,ic)
                    KmeansScaleFactorSqOld(ic)          = KmeansScaleFactorSq(ic)
                    KmeansMahalSqOld(1:np,ic)           = KmeansMahalSq(1:np,ic)
                    KmeansScaleFactorOld(ic)            = KmeansScaleFactor(ic)
                    KmeansCenterOld(1:nd,ic)            = KmeansCenter(1:nd,ic)
                    KmeansChoDiaOld(1:nd,ic)            = KmeansChoDia(1:nd,ic)
                    KmeansLogVolOld(ic)                 = KmeansLogVol(ic)

                end if

            end if blockStabilityCheck

            !if (KmeansLogVol(ic)<KmeansLogVolEstimate(ic)) then  ! enlarge the bounding ellipsoid
            !  !write(*,*) 'enlargement happening...'
            !  KmeansScaleFactorSq(ic) = KmeansScaleFactorSq(ic) * (KmeansLogVolEstimate(ic)/KmeansLogVol(ic))**(2._RK/nd)
            !  KmeansScaleFactor(ic) = sqrt(KmeansScaleFactorSq(ic))
            !  KmeansLogVol(ic) = KmeansLogVolEstimate(ic)
            !end if

#if defined DEBUG_ENABLED && false
            debugVolNorm: block
                real(RK) :: dummy
                dummy = log( real(KmeansSize(ic),kind=RK) / real(np,kind=RK) ) + parLogVol
                if ( KmeansLogVol(ic) < dummy ) then  ! enlarge the bounding ellipsoid
                    write(*,*) 'enlargement happening...'
                    write(*,*) 'enlargement happening...'
                    write(*,*) 'enlargement happening...'
                    write(*,*) 'enlargement happening...'
                    KmeansScaleFactorSq(ic) = KmeansScaleFactorSq(ic) * (dummy/KmeansLogVol(ic))**(2._RK/nd)
                    KmeansScaleFactor(ic) = sqrt(KmeansScaleFactorSq(ic))
                    KmeansLogVol(ic) = dummy
                end if
            end block debugVolNorm
#endif
            !end if blockStabilityCheck

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Reassign points to clusters if needed, based on the partitioning criterion: minimum-volume, maximum-density, ...
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            blockPartitionOptimization: if (partitionOptimizationRequested .and. recursionCounter < 2_IK) then

#if !defined MAXDEN && !defined MINVOL
                MahalSqWeight = ( KmeansLogVol / KmeansNumPointInside )**mahalSqWeightExponent
#elif defined MAXDEN
                MahalSqWeight = ( KmeansLogVol / KmeansNumPointInside )**2
#endif
                KmeansCenter(1:nd,1) = KmeansCenter(1:nd,1) * KmeansSize(1)  ! Now this is sum instead of mean
                KmeansCenter(1:nd,2) = KmeansCenter(1:nd,2) * KmeansSize(2)  ! Now this is sum instead of mean

                reclusteringNeeded = .false.

                ic = 1
                icOther = 2
                do ip = IpStart(ic), IpEnd(ic)
                !ip = IpStart(ic)
                !do
                    PartitionMembership(ip) = ic    ! This has not been assigned correctly in the past, here is the first attempt
#if defined MINVOL
                    if ( KmeansMahalSq(ip,ic) >= KmeansMahalSq(ip,icOther) ) then
#else
                    if ( MahalSqWeight(ic)*KmeansMahalSq(ip,ic) >= MahalSqWeight(icOther)*KmeansMahalSq(ip,icOther) ) then
#endif
                        PartitionMembership(ip) = icOther
                        KmeansSize(ic) = KmeansSize(ic) - 1
                        KmeansSize(icOther) = KmeansSize(icOther) + 1
                        KmeansCenter(1:nd,ic) = KmeansCenter(1:nd,ic) - Point(1:nd,ip)
                        KmeansCenter(1:nd,icOther) = KmeansCenter(1:nd,icOther) + Point(1:nd,ip)
                        !PointDum(1:nd) = Point(1:nd,ip)
                        !Point(1:nd,ip) = Point(1:nd,IpEnd(ic))
                        !Point(1:nd,IpEnd(ic)) = PointDum(1:nd)
                        !IpEnd(ic) = IpEnd(ic) - 1
                        reclusteringNeeded = .true.
                    end if
                    !if (ip == IpEnd(ic)) exit
                    !ip = ip + 1_IK
                end do

                ic = 2
                icOther = 1
                do ip = IpStart(ic), IpEnd(ic)
                    PartitionMembership(ip) = ic    ! This has not been assigned correctly in the past, here is the first attempt
#if defined MINVOL
                    if ( KmeansMahalSq(ip,ic) > KmeansMahalSq(ip,icOther) ) then
#else
                    if ( MahalSqWeight(ic)*KmeansMahalSq(ip,ic) > MahalSqWeight(icOther)*KmeansMahalSq(ip,icOther) ) then
#endif
                        PartitionMembership(ip) = icOther
                        KmeansSize(ic) = KmeansSize(ic) - 1
                        KmeansSize(icOther) = KmeansSize(icOther) + 1
                        KmeansCenter(1:nd,ic) = KmeansCenter(1:nd,ic) - Point(1:nd,ip)
                        KmeansCenter(1:nd,icOther) = KmeansCenter(1:nd,icOther) + Point(1:nd,ip)
                        reclusteringNeeded = .true.

                    end if
                end do

                ! Restart the process if anything has changed

write(*,*) "reclusteringNeeded, recursionCounter", reclusteringNeeded, recursionCounter
                blockReclusteringNeeded: if (reclusteringNeeded) then ! perform reassignment

                    blockIllegalClusterSize: if ( any(KmeansSize < minClusterSize) ) then

                        partitionFailureCount = partitionFailureCount + 1
                        if (partitionFailureCount > partitionMaxAllowedFailure) then
                            neopt = 1_IK
                            PartitionSize(1) = np
                            PartitionMembership = 1_IK
                            convergenceFailureCount = convergenceFailureCount + 1_IK
                            !debugpartitionFailureCount: block
                            !  write(*,*) 'partitionFailureCount: ', partitionFailureCount
                            !end block debugpartitionFailureCount
                            return
                        end if

                        if (kmeansClusteringFailed()) return

                        cycle loopRecursiveRunPartitionKernel

                    end if blockIllegalClusterSize

                    KmeansCenter(1:nd,1) = KmeansCenter(1:nd,1) / KmeansSize(1)  ! Now it is mean instead of sum
                    KmeansCenter(1:nd,2) = KmeansCenter(1:nd,2) / KmeansSize(2)  ! Now it is mean instead of sum

                    ! Reorder Point based on the identified clusters.

                    KmeansMemberCounter(1) = 0_IK
                    KmeansMemberCounter(2) = KmeansSize(1)
                    do ip = 1, np
                        KmeansMemberCounter(PartitionMembership(ip)) = KmeansMemberCounter(PartitionMembership(ip)) + 1_IK
                        KmeansPointIndex(KmeansMemberCounter(PartitionMembership(ip))) = PointIndex(ip)
                        KmeansPoint(1:nd,KmeansMemberCounter(PartitionMembership(ip))) = Point(1:nd,ip)
                    end do
                    Point = KmeansPoint             ! Point is now ordered
                    PointIndex = KmeansPointIndex   ! PointIndex is now ordered

                    do ic = 1, 2
                        !do ip = IpStart(ic), IpEnd(ic)
                        !    NormedPoint(1:nd,ip,ic) = Point(1:nd,ip) - KmeansCenter(1:nd,ic)
                        do ip = 1, np
                            NormedPoint(1:nd,ip,ic) = Point(1:nd,ip) - KmeansCenter(1:nd,PartitionMembership(ip))
                        end do
                    end do

                    cycle loopRecursiveRunPartitionKernel

                end if blockReclusteringNeeded

                KmeansCenter(1:nd,1) = KmeansCenter(1:nd,1) / KmeansSize(1)  ! Now it is mean instead of sum
                KmeansCenter(1:nd,2) = KmeansCenter(1:nd,2) / KmeansSize(2)  ! Now it is mean instead of sum

            end if blockPartitionOptimization

            exit loopRecursiveRunPartitionKernel ! mission accomplished, current volumes minimized.

        end do loopRecursiveRunPartitionKernel

#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Rescale the volumes, if needed and for as long as needed, based on the user-input volume estimate
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (parLogVolIsPresent) then
            do ic = 1, 2
#if defined DEBUG_ENABLED
                if (KmeansNumPointInside(ic) > np) then
                    write(*,*) "FATAL: Internal error occurred: KmeansNumPointInside(ic) > np", KmeansNumPointInside(ic), np
                    error stop
                end if
#endif
                KmeansLogVolEstimate(ic) = parLogVol + log( KmeansNumPointInside(ic) / real(np, kind = SPR) )
                if (KmeansLogVol(ic) < KmeansLogVolEstimate(ic)) then
                    KmeansScaleFactor(ic) = KmeansScaleFactor(ic) * exp( (KmeansLogVolEstimate(ic) - KmeansLogVol(ic)) / nd )
                    KmeansScaleFactorSq(ic) = KmeansScaleFactor(ic)**2
                    KmeansLogVol(ic) = KmeansLogVolEstimate(ic)
                end if
            end do
            boundedRegionIsTooLarge = PartitionLogVol(1) > logTightness + parLogVol
        else
            boundedRegionIsTooLarge = .false.
        end if

        ! @todo: There is room for improvement here. When inclusionFraction /= 0, KmeansLogVolEstimate could be recursively checked for potential rescaling.

        !blockScaleFactorRescale: do
        !    do ic = 1, 2
        !        KmeansLogVolEstimate(ic) = parLogVol + real( log( real(KmeansNumPointInside(ic), kind = SPR) / real(np, kind = SPR) ) , kind = RK )
        !    end do
        !    KmeansNumPointInside(ic) = KmeansSize(ic) + nint( inclusionFraction * count(KmeansMahalSq(IpStart(icOther):IpEnd(icOther),ic)<KmeansScaleFactorSq(ic)), IK )
        !end blockScaleFactorRescale

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Check if further clustering is warranted.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !write(*,*) logTightness, parLogVol, PartitionLogVol(1), getLogSumExp(2, KmeansLogVol)
        furtherClusteringCheck: if  (logTightness + getLogSumExp(2, KmeansLogVol) < PartitionLogVol(1) .or. boundedRegionIsTooLarge) then

            ! At least two clusters is better than one. now search for more child clusters.

            KmeansNemax(1) = max( KmeansSize(1) / ndEffectivePlusOne , 1_IK )
            KmeansNemax(2) = max( KmeansSize(2) / ndEffectivePlusOne , 1_IK )

            IpStart(1) = 1                  ; IpEnd(1) = KmeansSize(1)
            IpStart(2) = KmeansSize(1) + 1  ; IpEnd(2) = np

            loopSearchForChildren: do ic = 1, 2   ! Search for more grand-child clusters

                if (ic==1_IK) then
                    icStart = 1_IK
                    icEnd = KmeansNemax(1)
                else
                    icStart = KmeansNeopt(1) + 1_IK
                    icEnd = KmeansNeopt(1) + KmeansNemax(2)
                end if

                PartitionLogVol(icStart) = KmeansLogVol(ic)
                if (KmeansNemax(ic) > 1_IK) then ! .and. KmeansLogVol(ic)>1.1_RK*KmeansLogVolEstimate(ic)) then

                    if (parLogVolIsPresent) parLogVol = KmeansLogVolEstimate(ic)

                    call runPartitionKernel ( nemax = KmeansNemax(ic) & ! LCOV_EXCL_LINE
                                            , nd = nd & ! LCOV_EXCL_LINE
                                            , np = KmeansSize(ic) & ! LCOV_EXCL_LINE
                                            , partitionMaxAllowedFailure = partitionMaxAllowedFailure & ! LCOV_EXCL_LINE
                                            , partitionMaxAllowedRecursion = partitionMaxAllowedRecursion & ! LCOV_EXCL_LINE
                                            , partitionMaxAllowedKmeansFailure = partitionMaxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                            , partitionStabilizationRequested = partitionStabilizationRequested & ! LCOV_EXCL_LINE
                                            , partitionOptimizationRequested = partitionOptimizationRequested & ! LCOV_EXCL_LINE
#if !defined MAXDEN && !defined MINVOL
                                            , mahalSqWeightExponent = mahalSqWeightExponent & ! LCOV_EXCL_LINE
#endif
                                            , inclusionFraction = inclusionFraction & ! LCOV_EXCL_LINE
                                            , logTightness = logTightness & ! LCOV_EXCL_LINE
                                            , parLogVol = parLogVol & ! LCOV_EXCL_LINE
                                            , Point = Point(1:nd,IpStart(ic):IpEnd(ic)) & ! LCOV_EXCL_LINE
                                            , neopt = KmeansNeopt(ic) & ! LCOV_EXCL_LINE
                                            , PointIndex = PointIndex(IpStart(ic):IpEnd(ic)) & ! LCOV_EXCL_LINE
                                            , PartitionSize = PartitionSize(icStart:icEnd) & ! LCOV_EXCL_LINE
                                            , PartitionLogVol = PartitionLogVol(icStart:icEnd) & ! LCOV_EXCL_LINE
                                            , PartitionCenter = PartitionCenter(1:nd,icStart:icEnd) & ! LCOV_EXCL_LINE
                                            , PartitionChoDia = PartitionChoDia(1:nd,icStart:icEnd) & ! LCOV_EXCL_LINE
                                            , PartitionInvCovMat = PartitionInvCovMat(1:nd,1:nd,icStart:icEnd) & ! LCOV_EXCL_LINE
                                            , PartitionMembership = PartitionMembership(IpStart(ic):IpEnd(ic)) & ! LCOV_EXCL_LINE
                                            , PartitionChoLowCovUpp = PartitionChoLowCovUpp(1:nd,1:nd,icStart:icEnd) & ! LCOV_EXCL_LINE
                                            , convergenceFailureCount = convergenceFailureCount & ! LCOV_EXCL_LINE
                                            , numRecursiveCall = numRecursiveCall & ! LCOV_EXCL_LINE
                                            )
                    PartitionMembership(IpStart(ic):IpEnd(ic)) = PartitionMembership(IpStart(ic):IpEnd(ic)) + icStart - 1_IK

                else

                    KmeansNeopt(ic) = 1_IK
                    PartitionMembership(IpStart(ic):IpEnd(ic)) = icStart

                end if

                if (KmeansNeopt(ic) == 1_IK) then ! .or. KmeansNemax(ic) == 1_IK ) then ! There is only one cluster here

                    KmeansNeopt(ic) = 1_IK
                    PartitionSize(icStart) = KmeansSize(ic)
                    PartitionChoLowCovUpp(1:nd,1:nd,icStart) = KmeansChoLowCovUpp(1:nd,1:nd,ic)

                    ! Rescale the Cholesky factor and its diagonals such that they describe the bounded shape of the cluster

                    scaleFactorSqInverse = 1._RK / KmeansScaleFactorSq(ic)
                    do j = 1, nd
                        PartitionCenter(j,icStart) = KmeansCenter(j,ic)
                        PartitionChoDia(j,icStart) = KmeansChoDia(j,ic) * KmeansScaleFactor(ic)
                        PartitionInvCovMat(1:nd,j,icStart) = KmeansInvCovMat(1:nd,j,ic) * scaleFactorSqInverse
                        do i = j + 1, nd
                            PartitionChoLowCovUpp(i,j,icStart) = PartitionChoLowCovUpp(i,j,icStart) * KmeansScaleFactor(ic)
                        end do
                    end do

#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
                    block
                        integer :: ip
                        logical :: isInside
                        real(RK) :: mahalSq
                        real(RK), allocatable :: NormedPoint(:)
                        do ip = IpStart(ic), IpEnd(ic)
                            NormedPoint = Point(:,ip) - PartitionCenter(:,icStart)
                            mahalSq = dot_product(NormedPoint,matmul(PartitionInvCovMat(:,:,icStart),NormedPoint))
                            isInside = mahalSq - 1._RK <= 1.e-6_RK
                            if (.not. isInside) then
                                write(*,"(*(g0.15,:,' '))") new_line("a"), "FATAL - POINT NOT INSIDE!, MAHAL = ", sqrt(mahalSq), new_line("a")
                                write(*,"(60(g0,:,','))") "numRecursiveCall", numRecursiveCall
                                write(*,"(60(g0,:,','))") "KmeansNemax(icStart), icStart", KmeansNemax(icStart), icStart
                                write(*,"(60(g0,:,','))") "ip, IpStart, IpEnd, size(PartitionMembership)", ip, IpStart(icStart), IpEnd(icStart), size(PartitionMembership(IpStart(icStart):IpEnd(icStart)))
                                write(*,"(60(g0,:,','))") "PartitionMembership", PartitionMembership(IpStart(icStart):IpEnd(icStart))
                                write(*,"(60(g0,:,','))") "PartitionChoLowCovUpp", PartitionChoLowCovUpp(:,:,icStart)
                                write(*,"(60(g0,:,','))") "PartitionInvCovMat", PartitionInvCovMat(:,:,icStart)
                                write(*,"(60(g0,:,','))") "PartitionChoDia", PartitionChoDia(:,icStart)
                                write(*,"(60(g0,:,','))") "NormedPoint", NormedPoint
                                error stop
                            end if
                        end do
                    end block
#endif

                end if

            end do loopSearchForChildren

            neopt = sum(KmeansNeopt)

        else furtherClusteringCheck  ! one cluster is better

            neopt = 1_IK
            PartitionSize(1) = np
            PartitionMembership = 1_IK
            return

        end if furtherClusteringCheck

#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
        block
            integer :: ic, ip
            logical :: isInside
            real(RK) :: mahalSq
            real(RK), allocatable :: NormedPoint(:)
            do ip = 1, np
                ic = PartitionMembership(ip)
                NormedPoint = Point(:,ip) - PartitionCenter(:,ic)
                mahalSq = dot_product(NormedPoint,matmul(PartitionInvCovMat(:,:,ic),NormedPoint))
                isInside = mahalSq - 1._RK <= 1.e-6_RK
                if (.not. isInside) then
                    write(*,"(*(g0.15,:,' '))") new_line("a"), "FATAL - POINT NOT INSIDE!, MAHALSQ = ", mahalSq, new_line("a")
                    write(*,"(60(g0,:,','))") "numRecursiveCall, ic, size(PartitionMembership)", numRecursiveCall, ic, size(PartitionMembership)
                    write(*,"(60(g0,:,','))") "PartitionMembership", PartitionMembership
                    write(*,"(60(g0,:,','))") "PartitionChoLowCovUpp", PartitionChoLowCovUpp(:,:,:) ! ic)
                    write(*,"(60(g0,:,','))") "PartitionInvCovMat", PartitionInvCovMat(:,:,ic)
                    write(*,"(60(g0,:,','))") "NormedPoint", NormedPoint
                    error stop
                end if
            end do
            !write(*,*) 'nemax: ', nemax
            !write(*,*) 'neopt: ', neopt
            !write(*,*) 'Membership: ', PartitionMembership(:)
            !write(*,*) 'PartitionSize: ', PartitionSize(1:neopt)
            !write(*,*) 'PartitionCenter: ', PartitionCenter(:,1:neopt)
            !read(*,*)
        end block
#endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    contains

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !> This is an internal function. 
        !> Compute the Kmeans clustering of the points, rearrange the points. 
        !> Return `failed = .false.` if the procedure fails to converge.
        !> The following quantities are computed by this procedure.
        !> KmeansSize
        !> PartitionMembership
        !> potential
        !> zcsIteration
        !> NormedPoint
        function kmeansClusteringFailed() result(failed)

            implicit none
            logical :: failed

            kmeansFailureCount = 0_IK
            loopRunKmeans: do
                call runKmeans  ( nd = nd & ! LCOV_EXCL_LINE
                                , np = np & ! LCOV_EXCL_LINE
                                , nc = 2_IK & ! LCOV_EXCL_LINE
                                , Point = Point & ! LCOV_EXCL_LINE
                                , Size = KmeansSize & ! LCOV_EXCL_LINE
                                , Center = KmeansCenter & ! LCOV_EXCL_LINE
                                , NormedPoint = NormedPoint & ! LCOV_EXCL_LINE
                                , Membership = PartitionMembership & ! LCOV_EXCL_LINE
                                , MinDistanceSq = KmeansMinDistanceSq & ! LCOV_EXCL_LINE
                                , potential = potential & ! LCOV_EXCL_LINE
                                , nzsci = zcsIteration & ! LCOV_EXCL_LINE
                                , niter = iteration & ! LCOV_EXCL_LINE
                                , Err = Err & ! LCOV_EXCL_LINE
                                !, niterMax & ! LCOV_EXCL_LINE
                                !, nzsciMax & ! LCOV_EXCL_LINE
                                !, relTol & ! LCOV_EXCL_LINE
                                )
                failed = Err%occurred .or. any(KmeansSize < minClusterSize)
                if (failed) then
                    kmeansFailureCount = kmeansFailureCount + 1
                    failed = kmeansFailureCount > partitionMaxAllowedKmeansFailure
                    if (failed) then
                        neopt = 1_IK
                        PartitionSize(1) = np
                        PartitionMembership = 1_IK
                        convergenceFailureCount = convergenceFailureCount + 1_IK
                        !debugKmeansFailureCounter: block
                        !  write(*,*) 'kmeansFailureCount: ', kmeansFailureCount
                        !end block debugKmeansFailureCounter
                        return
                    end if
                    cycle loopRunKmeans
                end if
                exit loopRunKmeans
            end do loopRunKmeans

            ! Reorder Point based on the identified clusters.

            KmeansMemberCounter(1) = 0_IK
            KmeansMemberCounter(2) = KmeansSize(1)
            do ip = 1, np
                KmeansMemberCounter(PartitionMembership(ip)) = KmeansMemberCounter(PartitionMembership(ip)) + 1_IK
                KmeansNormedPoint(1:nd,KmeansMemberCounter(PartitionMembership(ip)),1:2) = NormedPoint(1:nd,ip,1:2)
                KmeansPointIndex(KmeansMemberCounter(PartitionMembership(ip))) = PointIndex(ip)
                KmeansPoint(1:nd,KmeansMemberCounter(PartitionMembership(ip))) = Point(1:nd,ip)
            end do
            Point = KmeansPoint             ! Point is now ordered
            PointIndex = KmeansPointIndex   ! PointIndex is now ordered
            NormedPoint = KmeansNormedPoint ! NormedPoint is now ordered
            KmeansLogVolOld = HUGE_RK       ! This is important for the first iteration of loopRecursiveRunPartitionKernel.

        end function kmeansClusteringFailed

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !> This is an internal procedure. Compute the following properties implicitly and return nothing.
        !> IpEnd
        !> IpStart
        !> KmeansChoDia
        !> KmeansCenter
        !> KmeansLogVol
        !> KmeansMahalSq
        !> KmeansInvCovMat
        !> KmeansScaleFactor
        !> KmeansChoLowCovUpp
        !> KmeansScaleFactorSq
        !> KmeansNumPointInside
        subroutine computeBoundedResgion()

            implicit none

            IpStart(1) = 1_IK                   ; IpEnd(1) = KmeansSize(1)
            IpStart(2) = KmeansSize(1) + 1_IK   ; IpEnd(2) = np

            loopComputeScaleFactor: do ic = 1, 2

                icOther = 3 - ic

                ! Compute the upper covariance matrix of the cluster covariance matrices.

                KmeansChoLowCovUpp(1:nd,1:nd,ic) = 0._RK
                do ip = IpStart(ic), IpEnd(ic)
                    do j = 1, nd
                        KmeansChoLowCovUpp(1:j,j,ic) = KmeansChoLowCovUpp(1:j,j,ic) + NormedPoint(1:j,ip,ic) * NormedPoint(j,ip,ic)
                    end do
                end do
                KmeansChoLowCovUpp(1:nd,1:nd,ic) = KmeansChoLowCovUpp(1:nd,1:nd,ic) / real(KmeansSize(ic)-1, kind = RK)

                ! Compute the Cholesky Factor of the cluster covariance matrices.

                call getCholeskyFactor(nd, KmeansChoLowCovUpp(1:nd,1:nd,ic), KmeansChoDia(1:nd,ic))
                if (KmeansChoDia(1,ic)<0._RK) then
                    ! LCOV_EXCL_START
                    write(*,*) "    FATAL: runPartitionKernel() @ getSamCholFac() @ getCholeskyFactor() failed."
                    write(*,*) "           Cluster id, dimension, size:"
                    write(*,*) "           ", ic, nd, KmeansSize(ic)
                    error stop
                    ! LCOV_EXCL_STOP
                end if

#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
                block
                    use Statistics_mod, only: getSamCholFac
                    real(RK), parameter :: TOLERANCE = 1.e-2_RK
                    real(RK) :: ChoLowDiff(nd,nd), ChoDiaDiff(nd)
                    real(RK) :: ChoLow(nd,nd), ChoDia(nd)
                    call getSamCholFac  ( nd = nd & ! LCOV_EXCL_LINE
                                        , np = KmeansSize(ic) & ! LCOV_EXCL_LINE
                                        , Mean = KmeansCenter(1:nd,ic) & ! LCOV_EXCL_LINE
                                        , Point = Point(1:nd,IpStart(ic):IpEnd(ic)) & ! LCOV_EXCL_LINE
                                        , CholeskyLower = ChoLow & ! LCOV_EXCL_LINE
                                        , CholeskyDiago = ChoDia & ! LCOV_EXCL_LINE
                                        )
                    ChoLowDiff = 2 * abs(KmeansChoLowCovUpp(1:nd,1:nd,ic)-ChoLow) / abs(KmeansChoLowCovUpp(1:nd,1:nd,ic)+ChoLow)
                    if ( any(ChoLowDiff > TOLERANCE) ) then
                    ! LCOV_EXCL_START
                        write(*,*)
                        write(*,*) "KmeansChoLowCovUpp(1:nd,1:nd,ic)-ChoLow) > TOLERANCE"
                        write(*,*) "TOLERANCE"
                        write(*,*) TOLERANCE
                        write(*,*) "ChoLowDiff"
                        write(*,*) ChoLowDiff
                        write(*,*) "ChoLow"
                        write(*,*) ChoLow
                        write(*,*) "KmeansChoLowCovUpp(1:nd,1:nd,ic)"
                        write(*,*) KmeansChoLowCovUpp(1:nd,1:nd,ic)
                        write(*,*)
                        !error stop
                    end if
                    ! LCOV_EXCL_STOP
                    ChoDiaDiff = 2 * abs(KmeansChoDia(1:nd,ic)-ChoDia) / abs(KmeansChoDia(1:nd,ic)+ChoDia)
                    if ( any(ChoDiaDiff > TOLERANCE) ) then
                    ! LCOV_EXCL_START
                        write(*,*)
                        write(*,*) "KmeansChoDia(1:nd,ic)-ChoDia) > TOLERANCE"
                        write(*,*) "TOLERANCE"
                        write(*,*) TOLERANCE
                        write(*,*) "ChoDiaDiff"
                        write(*,*) ChoDiaDiff
                        write(*,*) "ChoDia"
                        write(*,*) ChoDia
                        write(*,*) "KmeansChoDia(1:nd,1:nd,ic)"
                        write(*,*) KmeansChoDia(1:nd,ic)
                        write(*,*)
                        !error stop
                    end if
                    ! LCOV_EXCL_STOP
                end block
#endif

                ! Compute the inverse of the cluster covariance matrices.

                KmeansInvCovMat(1:nd,1:nd,ic) = getInvMatFromCholFac( nd, KmeansChoLowCovUpp(1:nd,1:nd,ic), KmeansChoDia(1:nd,ic) )

                ! Compute the MahalSq of as many points as needed.

                do concurrent(ip = 1:np)
                    KmeansMahalSq(ip,ic) = dot_product( KmeansNormedPoint(1:nd,ip,ic) , matmul(KmeansInvCovMat(1:nd,1:nd,ic), KmeansNormedPoint(1:nd,ip,ic)) )
                end do

#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
                do ip = 1, np
                    KmeansMahalSq(ip,ic) = dot_product( KmeansNormedPoint(1:nd,ip,ic) , matmul(KmeansInvCovMat(1:nd,1:nd,ic), KmeansNormedPoint(1:nd,ip,ic)) )
                    if (KmeansMahalSq(ip,ic)<0._RK) then
                    ! LCOV_EXCL_START
                        KmeansMahalSq(1,ic) = -1._RK
                        write(*,*) "KmeansMahalSq(ip)<0._RK", ip, KmeansMahalSq(ip,ic)
                        error stop
                    end if
                    ! LCOV_EXCL_STOP
                    block
                        use Statistics_mod, only: getMahalSq
                        real(RK), parameter :: TOLERANCE = 1.e-2_RK
                        real(RK) :: dum, diff
                        dum = getMahalSq( nd, KmeansCenter(1:nd,ic), KmeansInvCovMat(1:nd,1:nd,ic), Point(1:nd,ip) )
                        diff = 2 * abs(dum-KmeansMahalSq(ip,ic)) / abs(dum+KmeansMahalSq(ip,ic))
                        if (diff > TOLERANCE) then
                            ! LCOV_EXCL_START
                            write(*,*)
                            write(*,*) "dum/=KmeansMahalSq(ip,ic)", ip, KmeansMahalSq(ip,ic), dum
                            write(*,*) "TOLERANCE"
                            write(*,*) TOLERANCE
                            write(*,*) "diff"
                            write(*,*) diff
                            write(*,*)
                            !error stop
                            ! LCOV_EXCL_STOP
                        end if
                    end block
                end do
#endif
                ! Compute the scaleFcator of the bounding region.

                KmeansScaleFactorSq(ic) = maxval(KmeansMahalSq(IpStart(ic):IpEnd(ic),ic)); KmeansScaleFactor(ic) = sqrt(KmeansScaleFactorSq(ic))
                KmeansMahalSq(1:np,ic) = KmeansMahalSq(1:np,ic) / KmeansScaleFactorSq(ic)
                KmeansLogVol(ic) = nd * log(KmeansScaleFactor(ic)) + sum( log(KmeansChoDia(1:nd,ic)) )

            end do loopComputeScaleFactor

            if (isZeroInclusionFraction) then
                KmeansNumPointInside(ic) = KmeansSize(ic)
            else
                KmeansNumPointInside(ic) = KmeansSize(ic) + nint( inclusionFraction * count(KmeansMahalSq(IpStart(icOther):IpEnd(icOther),ic)<1._RK), IK )
            end if

        end subroutine computeBoundedResgion

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end subroutine runPartitionKernel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine writePartitionOOP(self, fileUnit, nd, np, Point)
        use Constants_mod, only: IK, RK
        implicit none
        class(Partition_type)   , intent(in)    :: self
        integer(IK)             , intent(in)    :: fileUnit, nd, np
        real(RK)                , intent(in)    :: Point(nd,np)
        call writePartition ( fileUnit = fileUnit & ! LCOV_EXCL_LINE
                            , nd = nd & ! LCOV_EXCL_LINE
                            , np = np & ! LCOV_EXCL_LINE
                            , nc = self%neopt & ! LCOV_EXCL_LINE
                            , PartitionSize = self%Size & ! LCOV_EXCL_LINE
                            , PartitionCenter = self%Center & ! LCOV_EXCL_LINE
                            , PartitionLogVol = self%LogVolNormed & ! LCOV_EXCL_LINE
                            , PartitionMembership = self%Membership & ! LCOV_EXCL_LINE
                            , PartitionChoLowCovUpp = self%ChoLowCovUpp & ! LCOV_EXCL_LINE
                            , PartitionChoDia = self%ChoDia & ! LCOV_EXCL_LINE
                            , Point = Point & ! LCOV_EXCL_LINE
                            )
    end subroutine writePartitionOOP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine writePartition   ( fileUnit & ! LCOV_EXCL_LINE
                                , nd,np,nc & ! LCOV_EXCL_LINE
                                , PartitionSize & ! LCOV_EXCL_LINE
                                , PartitionCenter & ! LCOV_EXCL_LINE
                                , PartitionLogVol & ! LCOV_EXCL_LINE
                                , PartitionMembership & ! LCOV_EXCL_LINE
                                , PartitionChoLowCovUpp & ! LCOV_EXCL_LINE
                                , PartitionChoDia & ! LCOV_EXCL_LINE
                                , Point & ! LCOV_EXCL_LINE
                                )
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK)     , intent(in)    :: fileUnit
        integer(IK)     , intent(in)    :: nd, np, nc
        integer(IK)     , intent(in)    :: PartitionSize(nc)
        integer(IK)     , intent(in)    :: PartitionMembership(np)
        real(RK)        , intent(in)    :: PartitionLogVol(nc)
        real(RK)        , intent(in)    :: PartitionCenter(nd,nc)
        real(RK)        , intent(in)    :: PartitionChoLowCovUpp(nd,nd,nc)
        real(RK)        , intent(in)    :: PartitionChoDia(nd,nc)
        real(RK)        , intent(in)    :: Point(nd,np)
        character(*)    , parameter     :: fileFormat = "(*(g0.8,:,','))"
        integer(IK)                     :: i, j, ic

        write(fileUnit,"(A)") "nd, np, nc"
        write(fileUnit,fileFormat) nd, np, nc

        write(fileUnit,"(A)") "Size"
        write(fileUnit,fileFormat) PartitionSize

        write(fileUnit,"(A)") "Center"
        write(fileUnit,fileFormat) PartitionCenter

        write(fileUnit,"(A)") "LogVolume"
        write(fileUnit,fileFormat) PartitionLogVol

        write(fileUnit,"(A)") "CholeskyLower"
        write(fileUnit,fileFormat) ((PartitionChoDia(j,ic), (PartitionChoLowCovUpp(i,j,ic), i=j+1,nd), j=1,nd), ic=1,nc)

!write(*,fileFormat)
!write(*,fileFormat)
!write(*,fileFormat) ((PartitionChoDia(j,ic), (PartitionChoLowCovUpp(i,j,ic), i=j+1,nd), j=1,nd), ic=1,nc)
!write(*,fileFormat) ((PartitionChoDia(j,ic), j=1,nd), ic=1,nc)
!write(*,fileFormat) (((PartitionChoLowCovUpp(i,j,ic), i=j+1,nd), j=1,nd), ic=1,nc)
!write(*,fileFormat)
!write(*,fileFormat)
!write(*,fileFormat) (((PartitionChoLowCovUpp(i,j,ic), i=1,j), j=1,nd), ic=1,nc)
!block
!real(RK), allocatable :: cholow(:,:), covmat(:,:)
!allocate(cholow(nd,nd), source = 0._RK)
!do j = 1, nd
!    cholow(j,j) = PartitionChoDia(j,1)
!    do i = j+1,nd
!        cholow(i,j) = PartitionChoLowCovUpp(i,j,1)
!    end do
!end do
!covmat = matmul(cholow,transpose(cholow))
!write(*,fileFormat) ((covmat(i,j), i=1,j), j=1,nd)
!end block

        write(fileUnit,"(A)") "Point"
        write(fileUnit,fileFormat) Point

        write(fileUnit,"(A)") "Membership"
        write(fileUnit,fileFormat) PartitionMembership

    end subroutine writePartition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#undef OPTIMIZATION_ENABLED