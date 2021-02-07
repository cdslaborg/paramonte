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

    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type

    implicit none

#if defined KMEANS
    character(*), parameter :: MODULE_NAME = "@PartitionKmeans_mod"
#elif defined MAXDEN
    character(*), parameter :: MODULE_NAME = "@PartitionMaxDen_mod"
#elif defined MINVOL
    character(*), parameter :: MODULE_NAME = "@PartitionMinVol_mod"
#else
    character(*), parameter :: MODULE_NAME = "@Partition_mod"
#endif

    !> The `Partition_type` class. Partitions an input array `Point(nd,np)` with `nd` attributes and `np` observations (points).
    type :: Partition_type
        integer(IK)                 :: nd,np                                !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        integer(IK)                 :: nc,nt                                !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        integer(IK)                 :: nemax                                !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        integer(IK)                 :: neopt                                !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        integer(IK)                 :: numRecursiveCall                     !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        integer(IK)                 :: convergenceFailureCount              !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        integer(IK)                 :: partitionMaxAllowedKmeansFailure     !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        integer(IK)                 :: partitionMaxAllowedKmeansRecursion   !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
#if !defined KMEANS
        integer(IK)                 :: partitionMaxAllowedRecursion         !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
#endif
#if !(defined KMEANS || defined MAXDEN || defined MINVOL)
        real(RK)                    :: mahalSqWeightExponent                !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
#endif
        real(RK)                    :: partitionKmeansRelTol                !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        real(RK)                    :: inclusionFraction                    !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        real(RK)                    :: parLogVolNormed                      !< The logarithm of the volume of the input set of points normalized to the volume of the unit nd-ball.
        real(RK)                    :: logTightness                         !< See the interface of [runPartitionKernel()](@ref runPartitionkernel).
        integer(IK) , allocatable   :: Membership(:)                        !< An array of size `(np)` representing the bounding-ellipsoid membership IDs of the corresponding data points.
        integer(IK) , allocatable   :: PointIndex(:)                        !< An array of size `(np)` representing the indices of the original data points before partitioning.
        real(RK)    , allocatable   :: LogVolNormed(:)                      !< An array of size `(nemax)` representing the log-volumes of the corresponding bounding ellipsoids.
        real(RK)    , allocatable   :: ChoLowCovUpp(:,:,:)                  !< An array of size `(nd,nd,nemax)` representing the Cholesky lower triangle, diagonal, and covariance matrices of the bounding ellipsoids.
        real(RK)    , allocatable   :: InvCovMat(:,:,:)                     !< An array of size `(nd,nd,nemax)` representing the full symmetric inverse covariance matrices of the bounding ellipsoids.
        real(RK)    , allocatable   :: ChoDia(:,:)                          !< An array of size `(nd,nemax)` representing the Cholesky lower triangle, diagonal, and covariance matrices of the bounding ellipsoids.
        real(RK)    , allocatable   :: Center(:,:)                          !< An array of size `(nd,nemax)` representing the centers of the bounding ellipsoids.
        integer(IK) , allocatable   :: Size(:)                              !< An array of size `(nemax)` representing the sizes of the corresponding bounding ellipsoids.
        type(Err_type)              :: Err                                  !< An object of class [Err_type](@ref err_mod::err_type) containing information about error occurrence.
    contains
        procedure, pass :: write => writePartitionOOP
    end type Partition_type

    !interface Partition_type ! stupid gfortran complains about this interface
    !    module procedure :: constructPartition
    !end interface Partition_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getPartition   ( nd, np, Point & ! LCOV_EXCL_LINE
                            , nc, nt, nemax & ! LCOV_EXCL_LINE
#if !defined KMEANS
                            , partitionMaxAllowedRecursion & ! LCOV_EXCL_LINE
#endif
#if !(defined KMEANS || defined MAXDEN || defined MINVOL)
                            , mahalSqWeightExponent & ! LCOV_EXCL_LINE
#endif
                            , inclusionFraction & ! LCOV_EXCL_LINE
                            , logTightness & ! LCOV_EXCL_LINE
                            , trimEnabled & ! LCOV_EXCL_LINE
                            , parLogVolNormed & ! LCOV_EXCL_LINE
                            , partitionKmeansRelTol & ! LCOV_EXCL_LINE
                            , partitionMaxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                            , partitionMaxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                            ) result(Partition)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getPartition
#endif
        use Constants_mod, only: IK, RK !, NEGINF_RK

        implicit none

        integer(IK) , intent(in)                :: nd,np
        real(RK)    , intent(inout)             :: Point(nd,np) ! do not move this up. This must appear after nd, np declaration.
        integer(IK) , intent(in), optional      :: nc, nt, nemax
#if !defined KMEANS
        integer(IK) , intent(in), optional      :: partitionMaxAllowedRecursion
#endif
#if !(defined KMEANS || defined MAXDEN || defined MINVOL)
        real(RK)    , intent(in), optional      :: mahalSqWeightExponent
#endif
        real(RK)    , intent(in), optional      :: inclusionFraction
        real(RK)    , intent(in), optional      :: logTightness
        logical     , intent(in), optional      :: trimEnabled
        real(RK)    , intent(inout), optional   :: parLogVolNormed
        real(RK)    , intent(in), optional      :: partitionKmeansRelTol
        integer(IK) , intent(in), optional      :: partitionMaxAllowedKmeansFailure
        integer(IK) , intent(in), optional      :: partitionMaxAllowedKmeansRecursion
        type(Partition_type)                    :: Partition

        logical                                 :: trimEnabledDefault

        Partition%nd = nd
        Partition%np = np
        Partition%nc = 2_IK; if (present(nc)) Partition%nc = nc
        Partition%nt = 1_IK; if (present(nt)) Partition%nt = nt
        Partition%Err%occurred = .false.
        trimEnabledDefault = .false.; if (present(trimEnabled)) trimEnabledDefault = trimEnabled

        if (present(nemax)) then
            Partition%nemax = nemax
        else
            Partition%nemax = np / (nd + 1)
        end if
        Partition%nemax = max(Partition%nc, Partition%nemax)

        !if (present(nemax)) then
            !if (nemax > Partition%nemax) then
            !    Partition%Err%occurred = .true.
            !    Partition%Err%msg = "The input value for `nemax` cannot be larger than `np / (nd + 1)`."
            !    return
            !else
            !    Partition%nemax = nemax
            !    !if (Partition%nemax<nemax) Partition%nemax = nemax
            !end if
        !end if

        Partition%partitionMaxAllowedKmeansRecursion = 300_IK; if (present(partitionMaxAllowedKmeansRecursion)) Partition%partitionMaxAllowedKmeansRecursion = partitionMaxAllowedKmeansRecursion
        Partition%partitionMaxAllowedKmeansFailure = 10_IK; if (present(partitionMaxAllowedKmeansFailure)) Partition%partitionMaxAllowedKmeansFailure = partitionMaxAllowedKmeansFailure
#if !defined KMEANS
        Partition%partitionMaxAllowedRecursion = 50_IK; if (present(partitionMaxAllowedRecursion)) Partition%partitionMaxAllowedRecursion = partitionMaxAllowedRecursion
#endif
#if !(defined KMEANS || defined MAXDEN || defined MINVOL)
        Partition%mahalSqWeightExponent = 0._RK; if (present(mahalSqWeightExponent)) Partition%mahalSqWeightExponent = mahalSqWeightExponent
#endif
        Partition%partitionKmeansRelTol = 1.e-4_RK; if (present(partitionKmeansRelTol)) Partition%partitionKmeansRelTol = partitionKmeansRelTol
        Partition%inclusionFraction = 0._RK; if (present(inclusionFraction)) Partition%inclusionFraction = inclusionFraction
        Partition%logTightness = log(1.3_RK); if (present(logTightness)) Partition%logTightness = logTightness
        Partition%parLogVolNormed = -huge(1._RK); if (present(parLogVolNormed)) Partition%parLogVolNormed = parLogVolNormed

        allocate( Partition%Size         (Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%Center       (nd,Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%ChoDia       (nd,Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%InvCovMat    (nd,nd,Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%ChoLowCovUpp (nd,nd,Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%LogVolNormed (Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%Membership   (np) & ! LCOV_EXCL_LINE
                , Partition%PointIndex   (np) & ! LCOV_EXCL_LINE
                )

        call runPartition   ( Point = Point & ! LCOV_EXCL_LINE
                            , nd = nd & ! LCOV_EXCL_LINE
                            , np = np & ! LCOV_EXCL_LINE
                            , nc = Partition%nc & ! LCOV_EXCL_LINE
                            , nt = Partition%nt & ! LCOV_EXCL_LINE
                            , nemax = Partition%nemax & ! LCOV_EXCL_LINE
#if !defined KMEANS
                            , partitionMaxAllowedRecursion = Partition%partitionMaxAllowedRecursion & ! LCOV_EXCL_LINE
#endif
#if !(defined KMEANS || defined MAXDEN || defined MINVOL)
                            , mahalSqWeightExponent = Partition%mahalSqWeightExponent & ! LCOV_EXCL_LINE
#endif
                            , inclusionFraction = Partition%inclusionFraction & ! LCOV_EXCL_LINE
                            , logTightness = Partition%logTightness & ! LCOV_EXCL_LINE
                            , neopt = Partition%neopt & ! LCOV_EXCL_LINE
                            , PointIndex = Partition%PointIndex & ! LCOV_EXCL_LINE
                            , PartitionSize = Partition%Size & ! LCOV_EXCL_LINE
                            , PartitionCenter = Partition%Center & ! LCOV_EXCL_LINE
                            , PartitionChoDia = Partition%ChoDia & ! LCOV_EXCL_LINE
                            , PartitionInvCovMat = Partition%InvCovMat & ! LCOV_EXCL_LINE
                            , PartitionMembership = Partition%Membership & ! LCOV_EXCL_LINE
                            , PartitionChoLowCovUpp = Partition%ChoLowCovUpp & ! LCOV_EXCL_LINE
                            , PartitionLogVolNormed = Partition%LogVolNormed & ! LCOV_EXCL_LINE
                            , convergenceFailureCount = Partition%convergenceFailureCount & ! LCOV_EXCL_LINE
                            , numRecursiveCall = Partition%numRecursiveCall & ! LCOV_EXCL_LINE
                            , Err = Partition%Err & ! LCOV_EXCL_LINE
                            , parLogVolNormed = parLogVolNormed & ! LCOV_EXCL_LINE
                            , partitionKmeansRelTol = Partition%partitionKmeansRelTol & ! LCOV_EXCL_LINE
                            , partitionMaxAllowedKmeansFailure = Partition%partitionMaxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                            , partitionMaxAllowedKmeansRecursion = Partition%partitionMaxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                            )

        if (trimEnabledDefault) then
            Partition%Size           = Partition%Size         (1:Partition%neopt)
            Partition%Center         = Partition%Center       (1:nd,1:Partition%neopt)
            Partition%ChoDia         = Partition%ChoDia       (1:nd,1:Partition%neopt)
            Partition%InvCovMat      = Partition%InvCovMat    (1:nd,1:nd,1:Partition%neopt)
            Partition%ChoLowCovUpp   = Partition%ChoLowCovUpp (1:nd,1:nd,1:Partition%neopt)
            Partition%LogVolNormed   = Partition%LogVolNormed (1:Partition%neopt)
        end if

    end function getPartition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \warning
    !> On input, `nemax` must be always larger than or equal to `nc`.
    subroutine runPartition ( Point & ! LCOV_EXCL_LINE
                            , nd,np & ! LCOV_EXCL_LINE
                            , nc,nt & ! LCOV_EXCL_LINE
                            , nemax & ! LCOV_EXCL_LINE
#if !defined KMEANS
                            , partitionMaxAllowedRecursion & ! LCOV_EXCL_LINE
#endif
#if !(defined KMEANS || defined MAXDEN || defined MINVOL)
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
                            , PartitionInvCovMat & ! LCOV_EXCL_LINE
                            , PartitionMembership & ! LCOV_EXCL_LINE
                            , PartitionChoLowCovUpp & ! LCOV_EXCL_LINE
                            , PartitionLogVolNormed & ! LCOV_EXCL_LINE
                            , convergenceFailureCount & ! LCOV_EXCL_LINE
                            , numRecursiveCall & ! LCOV_EXCL_LINE
                            , Err & ! LCOV_EXCL_LINE
                            , parLogVolNormed & ! LCOV_EXCL_LINE
                            , partitionKmeansRelTol & ! LCOV_EXCL_LINE
                            , partitionMaxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                            , partitionMaxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
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
        integer(IK)     , intent(in)                :: nc,nt                                ! the number of subclusters in each partitioning step. The number of Kmeans clustering tries.
        integer(IK)     , intent(in)                :: nemax                                ! the maximum allowed number of ellipsoidal clusters, number of dimensions, number of points
        real(RK)        , intent(inout)             :: Point(nd,np)                         ! the input data (points) to be clustered. on output, it is ordered. This must appear after nd, np declaration.
#if !defined KMEANS
        integer(IK)     , intent(in)                :: partitionMaxAllowedRecursion
#endif
#if !(defined KMEANS || defined MAXDEN || defined MINVOL)
        real(RK)        , intent(in)                :: mahalSqWeightExponent
#endif
        real(RK)        , intent(in)                :: inclusionFraction                    ! the fraction of non-cluster-member points that ARE inside the cluster, to be considered in the estimation of the true volume of the cluster.
        real(RK)        , intent(in)                :: logTightness                         ! the factor deremining how large the bounding ellipsoid can be, compared to the true volume of the ellipsoid.
        integer(IK)     , intent(out)               :: neopt                                ! this variable is saved throughout recursive iterations of subroutine
        integer(IK)     , intent(out)               :: PointIndex(np)                       ! Vector containing the scrambled indices of Point elements. It represent the mapping from the input Point elements to output scrambled Point vector.
        integer(IK)     , intent(out)               :: PartitionSize(nemax)                 ! cluster sizes (nemax)
        integer(IK)     , intent(out)               :: PartitionMembership(np)              ! the PartitionMembership id of each point (representing the cluster number to which the point belongs)
        real(RK)        , intent(out)               :: PartitionCenter(nd,nemax)            ! Centers of the clusters (nd,nemax)
        real(RK)        , intent(out)               :: PartitionChoDia(nd,nemax)            ! CBE: Cluster Bounding Ellipsoid - On output: The diagonal elements of the Cholesky factor of the covariance matrix of the input points, on output, contains the cholesky factor of all found clusters.
        real(RK)        , intent(out)               :: PartitionLogVolNormed(nemax)         ! CBE: Cluster Bounding Ellipsoid - On output: Each element of the array contains the volume of the corresponding Cluster Bounding Ellipsoid.
        real(RK)        , intent(out)               :: PartitionInvCovMat(nd,nd,nemax)      ! CBE: Cluster Bounding Ellipsoid - On output: Full symmetric Inverse Covariance Matrices of the Bounding Ellipsoids of all clusters.
        real(RK)        , intent(out)               :: PartitionChoLowCovUpp(nd,nd,nemax)   ! CBE: Cluster Bounding Ellipsoid - On output: The lower triangle contains the Cholesky factor of the covariance matrix of the input points, the upper triangle maybe the CovMat itself (not needed). On output contains in the lower triangle, the Cholesky factor of all clusters found.
        integer(IK)     , intent(out)               :: convergenceFailureCount              ! The number of times the algorithm failed to converge.
        integer(IK)     , intent(out)               :: numRecursiveCall                     ! The number of recursive calls before convergence occurs.
        type(Err_type)  , intent(out)               :: Err                                  ! Error object containing error information.
        real(RK)        , intent(inout), optional   :: parLogVolNormed                            ! input log-estimate of the total volume of the points. Input zero or negative for the value to be ignored.
        real(RK)        , intent(in)   , optional   :: partitionKmeansRelTol                ! The `relTol` input argument to the `Kmeans_type()` constructor.
        integer(IK)     , intent(in)   , optional   :: partitionMaxAllowedKmeansFailure     ! The `nfailMax` input argument to the `Kmeans_type()` constructor.
        integer(IK)     , intent(in)   , optional   :: partitionMaxAllowedKmeansRecursion   ! The `niterMax` input argument to the `Kmeans_type()` constructor.

        logical                                     :: parLogVolNormedIsPresent
        logical                                     :: boundedRegionIsTooLarge
        real(RK)                                    :: scaleFactorSqInverse
        real(RK)                                    :: NormedPoint(nd,np)
        real(RK)                                    :: scaleFactorSq
        real(RK)                                    :: scaleFactor
        real(RK)                                    :: mahalSq
        integer(IK)                                 :: i, j, ip

        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME//"@runPartition()"

        parLogVolNormedIsPresent = present(parLogVolNormed)

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

        PartitionLogVolNormed(1) = sum( log( PartitionChoDia(1:nd,1) ) ) + nd * log(scaleFactor)

        ! Rescale the scaleFactor, if needed, to enlarge the bounded cluster in the future.
        ! This must be done here, although it will be used at the very end.
        ! Otherwise, PartitionLogVolNormed(1) will be potentially overwritten.

        if (parLogVolNormedIsPresent) then
            boundedRegionIsTooLarge = PartitionLogVolNormed(1) > parLogVolNormed
            if (boundedRegionIsTooLarge) then
                boundedRegionIsTooLarge = PartitionLogVolNormed(1) > logTightness + parLogVolNormed
            else
                scaleFactor = scaleFactor * exp( (parLogVolNormed - PartitionLogVolNormed(1)) / nd )
                scaleFactorSq = scaleFactor**2
                PartitionLogVolNormed(1) = parLogVolNormed
            end if
        else
            boundedRegionIsTooLarge = .true.
        end if

        if (nemax > 1_IK .and. boundedRegionIsTooLarge) then
            numRecursiveCall = 0_IK
            convergenceFailureCount = 0_IK
            call runPartitionKernel ( Point = Point & ! LCOV_EXCL_LINE
                                    , nd = nd & ! LCOV_EXCL_LINE
                                    , np = np & ! LCOV_EXCL_LINE
                                    , nc = min(nc,nemax) & ! LCOV_EXCL_LINE
                                    , nt = nt & ! LCOV_EXCL_LINE
                                    , nemax = nemax & ! LCOV_EXCL_LINE
#if !defined KMEANS
                                    , partitionMaxAllowedRecursion = partitionMaxAllowedRecursion & ! LCOV_EXCL_LINE
#endif
#if !(defined KMEANS || defined MAXDEN || defined MINVOL)
                                    , mahalSqWeightExponent = mahalSqWeightExponent & ! LCOV_EXCL_LINE
#endif
                                    , inclusionFraction = inclusionFraction & ! LCOV_EXCL_LINE
                                    , logTightness = logTightness & ! LCOV_EXCL_LINE
                                    , neopt = neopt & ! LCOV_EXCL_LINE
                                    , PointIndex = PointIndex & ! LCOV_EXCL_LINE
                                    , PartitionSize = PartitionSize & ! LCOV_EXCL_LINE
                                    , PartitionCenter = PartitionCenter & ! LCOV_EXCL_LINE
                                    , PartitionChoDia = PartitionChoDia & ! LCOV_EXCL_LINE
                                    , PartitionInvCovMat = PartitionInvCovMat & ! LCOV_EXCL_LINE
                                    , PartitionMembership = PartitionMembership & ! LCOV_EXCL_LINE
                                    , PartitionChoLowCovUpp = PartitionChoLowCovUpp & ! LCOV_EXCL_LINE
                                    , PartitionLogVolNormed = PartitionLogVolNormed & ! LCOV_EXCL_LINE
                                    , convergenceFailureCount = convergenceFailureCount & ! LCOV_EXCL_LINE
                                    , numRecursiveCall = numRecursiveCall & ! LCOV_EXCL_LINE
                                    , parLogVolNormed = parLogVolNormed & ! LCOV_EXCL_LINE
                                    , partitionKmeansRelTol = partitionKmeansRelTol & ! LCOV_EXCL_LINE
                                    , partitionMaxAllowedKmeansFailure = partitionMaxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                    , partitionMaxAllowedKmeansRecursion = partitionMaxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
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
                                            , nc & ! LCOV_EXCL_LINE
                                            , nt & ! LCOV_EXCL_LINE
                                            , nemax & ! LCOV_EXCL_LINE
#if !defined KMEANS
                                            , partitionMaxAllowedRecursion & ! LCOV_EXCL_LINE
#endif
#if !(defined KMEANS || defined MAXDEN || defined MINVOL)
                                            , mahalSqWeightExponent & ! LCOV_EXCL_LINE
#endif
                                            , inclusionFraction & ! LCOV_EXCL_LINE
                                            , logTightness & ! LCOV_EXCL_LINE
                                            , neopt & ! LCOV_EXCL_LINE
                                            , PointIndex & ! LCOV_EXCL_LINE
                                            , PartitionSize & ! LCOV_EXCL_LINE
                                            , PartitionCenter & ! LCOV_EXCL_LINE
                                            , PartitionChoDia & ! LCOV_EXCL_LINE
                                            , PartitionInvCovMat & ! LCOV_EXCL_LINE
                                            , PartitionMembership & ! LCOV_EXCL_LINE
                                            , PartitionChoLowCovUpp & ! LCOV_EXCL_LINE
                                            , PartitionLogVolNormed & ! LCOV_EXCL_LINE
                                            , convergenceFailureCount & ! LCOV_EXCL_LINE
                                            , numRecursiveCall & ! LCOV_EXCL_LINE
                                            , parLogVolNormed & ! LCOV_EXCL_LINE
                                            , partitionKmeansRelTol & ! LCOV_EXCL_LINE
                                            , partitionMaxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                            , partitionMaxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                                            )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: runPartitionKernel
#endif
        use Constants_mod, only: IK, RK, SPR, HUGE_IK, HUGE_RK
        use Kmeans_mod, only: Kmeans_type !, getKmeans, runKmeans
        use Math_mod, only: getLogSumExp

        implicit none

        integer(IK) , intent(in)                :: nemax                            ! maximum number of clusters determined from the nd and np
        integer(IK) , intent(in)                :: nd, np, nc, nt                   ! number of dimensions, number of points, number of clusters, number of tries
#if !defined KMEANS
        integer(IK) , intent(in)                :: partitionMaxAllowedRecursion
#endif
#if !(defined KMEANS || defined MAXDEN || defined MINVOL)
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
        real(RK)    , intent(inout)             :: PartitionLogVolNormed(nemax)         ! CBE: Cluster Bounding Ellipsoid - On input: The volume of the Bounding Ellipsoid of the parent cluster. On output: The Bounding volume of all child clusters.
        real(RK)    , intent(inout)             :: PartitionChoDia(nd,nemax)            ! CBE: Cluster Bounding Ellipsoid - On input: The diagonal elements of the Cholesky factor of the covariance matrix of the input points, on output, contains the cholesky factor of all found clusters.
        integer(IK) , intent(inout)             :: numRecursiveCall                     ! The number of recursive calls before convergence occurs.
        real(RK)    , intent(inout), optional   :: parLogVolNormed                      ! parent's sqrt(det(invCovMat)), where invCovMat is the inverse covariance matrix of the bounding ellipsoid of the input points.
        real(RK)    , intent(in)   , optional   :: partitionKmeansRelTol
        integer(IK) , intent(in)   , optional   :: partitionMaxAllowedKmeansFailure
        integer(IK) , intent(in)   , optional   :: partitionMaxAllowedKmeansRecursion

        type(Kmeans_type)       :: Kmeans
        real(RK)                :: NormedPoint(nd)                      ! Vector of Maximum Mahalanobis distances for each cluster.
        real(RK)                :: KmeansLogVolEstimate(nc)             ! New Cluster Volume estimates.
        real(RK)                :: kmeansScaleFactor                    ! New Cluster Bounding Ellipsoid Volumes
        integer(IK)             :: KmeansNemax(nc)                      ! maximum allowed number of clusters in each of the two K-means clusters.
        integer(IK)             :: KmeansNeopt(0:nc)                    ! optimal number of clusters in each of the child clusters.
#if !defined KMEANS
        type(Kmeans_type)       :: KmeansBackup
#endif
#if !(defined KMEANS || defined MINVOL)
        real(RK)                :: MahalSqWeight(nc)                    ! MahalSqWeight definition of Lu et al 2007.
#endif
        integer(IK)             :: i, j, ip, ic, ipstart, ipend !, ndEffectivePlusOne
        integer(IK)             :: icstart,icend,recursionCounter
        integer(IK)             :: minClusterSize
        integer(IK)             :: nemaxRemained
        integer(IK)             :: icmin
        real(RK)                :: scaleFactorSqInverse
        logical                 :: boundedRegionIsTooLarge
        logical                 :: reclusteringNeeded
        logical                 :: parLogVolNormedIsPresent

        minClusterSize = nd + 1_IK
        numRecursiveCall = numRecursiveCall + 1

        parLogVolNormedIsPresent = present(parLogVolNormed)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Perform the initial Kmeans clustering.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Kmeans = Kmeans_type( nd = nd & ! LCOV_EXCL_LINE
                            , np = np & ! LCOV_EXCL_LINE
                            , nc = nc & ! LCOV_EXCL_LINE
                            , nt = nt & ! LCOV_EXCL_LINE
                            , Point = Point & ! LCOV_EXCL_LINE
                            , relTol = partitionKmeansRelTol & ! LCOV_EXCL_LINE
                            , niterMax = partitionMaxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                            , nfailMax = partitionMaxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                            , inclusionFraction = inclusionFraction & ! LCOV_EXCL_LINE
                            , minSize = minClusterSize & ! LCOV_EXCL_LINE
                            , propEnabled = .true. & ! LCOV_EXCL_LINE
                            , Index = PointIndex & ! LCOV_EXCL_LINE
                            )
        if (Kmeans%Err%occurred) then
            neopt = 1_IK
            PartitionSize(1) = np
            PartitionMembership = 1_IK
            return
        end if

#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
        loopComputeScaleFactor: do ic = 1, nc
                block
                    use Statistics_mod, only: getSamCholFac, getMahalSq
                    real(RK), parameter :: TOLERANCE = 1.e-2_RK
                    real(RK) :: ChoLowDiff(nd,nd), ChoDiaDiff(nd)
                    real(RK) :: ChoLow(nd,nd), ChoDia(nd)
                    real(RK) :: dum, diff
                    call getSamCholFac  ( nd = nd & ! LCOV_EXCL_LINE
                                        , np = Kmeans%Size(ic) & ! LCOV_EXCL_LINE
                                        , Mean = Kmeans%Center(1:nd,ic) & ! LCOV_EXCL_LINE
                                        , Point = Point(1:nd,Kmeans%Prop%CumSumSize(ic-1)+1:Kmeans%Prop%CumSumSize(ic)) & ! LCOV_EXCL_LINE
                                        , CholeskyLower = ChoLow & ! LCOV_EXCL_LINE
                                        , CholeskyDiago = ChoDia & ! LCOV_EXCL_LINE
                                        )
                    ChoLowDiff = 2 * abs(Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic)-ChoLow) / abs(Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic)+ChoLow)
                    if ( any(ChoLowDiff > TOLERANCE) ) then
                        ! LCOV_EXCL_START
                        write(*,*)
                        write(*,*) "Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic)-ChoLow) > TOLERANCE"
                        write(*,*) "TOLERANCE"
                        write(*,*) TOLERANCE
                        write(*,*) "ChoLowDiff"
                        write(*,*) ChoLowDiff
                        write(*,*) "ChoLow"
                        write(*,*) ChoLow
                        write(*,*) "Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic)"
                        write(*,*) Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic)
                        write(*,*)
                        !error stop
                        ! LCOV_EXCL_STOP
                    end if
                    ChoDiaDiff = 2 * abs(Kmeans%Prop%ChoDia(1:nd,ic)-ChoDia) / abs(Kmeans%Prop%ChoDia(1:nd,ic)+ChoDia)
                    if ( any(ChoDiaDiff > TOLERANCE) ) then
                        ! LCOV_EXCL_START
                        write(*,*)
                        write(*,*) "Kmeans%Prop%ChoDia(1:nd,ic)-ChoDia) > TOLERANCE"
                        write(*,*) "TOLERANCE"
                        write(*,*) TOLERANCE
                        write(*,*) "ChoDiaDiff"
                        write(*,*) ChoDiaDiff
                        write(*,*) "ChoDia"
                        write(*,*) ChoDia
                        write(*,*) "Kmeans%Prop%ChoDia(1:nd,1:nd,ic)"
                        write(*,*) Kmeans%Prop%ChoDia(1:nd,ic)
                        write(*,*)
                        !error stop
                        ! LCOV_EXCL_STOP
                    end if

                    do ip = 1, np
                        NormedPoint(1:nd) = Point(1:nd,ip) - Kmeans%Center(1:nd,ic)
                        Kmeans%Prop%MahalSq(ip,ic) = dot_product( NormedPoint , matmul(Kmeans%Prop%InvCovMat(1:nd,1:nd,ic), NormedPoint) )
                        if (Kmeans%Prop%MahalSq(ip,ic)<0._RK) then
                        ! LCOV_EXCL_START
                            Kmeans%Prop%MahalSq(1,ic) = -1._RK
                            write(*,*) "Kmeans%Prop%MahalSq(ip)<0._RK", ip, Kmeans%Prop%MahalSq(ip,ic)
                            error stop
                        end if
                        ! LCOV_EXCL_STOP
                        dum = getMahalSq( nd, Kmeans%Center(1:nd,ic), Kmeans%Prop%InvCovMat(1:nd,1:nd,ic), Point(1:nd,ip) )
                        diff = 2 * abs(dum-Kmeans%Prop%MahalSq(ip,ic)) / abs(dum+Kmeans%Prop%MahalSq(ip,ic))
                        if (diff > TOLERANCE) then
                            ! LCOV_EXCL_START
                            write(*,*)
                            write(*,*) "dum/=Kmeans%Prop%MahalSq(ip,ic)", ip, Kmeans%Prop%MahalSq(ip,ic), dum
                            write(*,*) "TOLERANCE"
                            write(*,*) TOLERANCE
                            write(*,*) "diff"
                            write(*,*) diff
                            write(*,*)
                            !error stop
                            ! LCOV_EXCL_STOP
                        end if
                    end do
                end block
                !Kmeans%Prop%MahalSq(1:np,ic) = Kmeans%Prop%MahalSq(1:np,ic) / Kmeans%Prop%ScaleFactorSq(ic)
            end do loopComputeScaleFactor
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! If optimization is enabled, reassign points to clusters based on minimum-volume, maximum-density, ...
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if !defined KMEANS

        recursionCounter = 0_IK
        !call KmeansBackup%Prop%allocate(nd, np, nc)
        !KmeansBackup%Prop%LogVolNormed = 2 * Kmeans%Prop%LogVolNormed

        loopRecursivePartitionOptimization: do

            ! If the maximum allowed number of recursion has reached, quit.

            recursionCounter = recursionCounter + 1_IK
            if(recursionCounter > partitionMaxAllowedRecursion) then
                write(*,*) "recursionCounter > partitionMaxAllowedRecursion: ", recursionCounter
                neopt = 1_IK
                PartitionSize(1) = np
                PartitionMembership = 1_IK
                convergenceFailureCount = convergenceFailureCount + 1_IK
                !error stop
                return
            end if

            ! Ensure the stability of the minimum bounding regions.

            !if ( sum(Kmeans%Prop%LogVolNormed) < sum(KmeansBackup%Prop%LogVolNormed) ) then ! This must never occur on the first try
            !    KmeansBackup = Kmeans ! Create a backup of the current clustering properties.
            !else
            !    Kmeans = KmeansBackup ! Restore the last clustering properties.
            !    exit loopRecursivePartitionOptimization
            !end if

#if !(defined KMEANS || defined MAXDEN || defined MINVOL)
            MahalSqWeight = ( Kmeans%Prop%LogVolNormed / Kmeans%Prop%EffectiveSize )**mahalSqWeightExponent
#elif defined MAXDEN
            MahalSqWeight = ( Kmeans%Prop%LogVolNormed / Kmeans%Prop%EffectiveSize )**2
#endif

            ! Convert Center to sum(Point).

            do concurrent(ic = 1:nc)
                Kmeans%Center(1:nd,ic) = Kmeans%Center(1:nd,ic) * Kmeans%Size(ic)
            end do

            ! Optimize the Kmeans clusters.

            reclusteringNeeded = .false.
            do ip = 1, np

                ! Determine the cluster to which the point is the closest.

                icmin = Kmeans%Membership(ip)
                do ic = 1, nc
                    if (ic/=icmin) then
#if defined MINVOL
                        if ( Kmeans%Prop%MahalSq(ip,ic) < Kmeans%Prop%MahalSq(ip,icmin) ) icmin = ic
#else
                        if ( Kmeans%Prop%MahalSq(ip,ic) * MahalSqWeight(ic) < Kmeans%Prop%MahalSq(ip,icmin) * MahalSqWeight(icmin) ) icmin = ic
#endif
                    end if
                end do

                ! Switch the membership of point if needed.

                if (icmin /= Kmeans%Membership(ip) .and. Kmeans%Size(Kmeans%Membership(ip)) > minClusterSize) then
                    Kmeans%Size(icmin) = Kmeans%Size(icmin) + 1
                    Kmeans%Size(Kmeans%Membership(ip)) = Kmeans%Size(Kmeans%Membership(ip)) - 1
                    Kmeans%Center(1:nd,icmin) = Kmeans%Center(1:nd,icmin) + Point(1:nd,ip)
                    Kmeans%Center(1:nd,Kmeans%Membership(ip)) = Kmeans%Center(1:nd,Kmeans%Membership(ip)) - Point(1:nd,ip)
                    Kmeans%Membership(ip) = icmin
                    reclusteringNeeded = .true.
                end if

            end do

            ! Convert sum(Point) back to center of clusters.

            do concurrent(ic = 1:nc)
                Kmeans%Center(1:nd,ic) = Kmeans%Center(1:nd,ic) / Kmeans%Size(ic)
            end do

            ! Restart the process if anything has changed.

            blockReclusteringNeeded: if (reclusteringNeeded) then ! perform reassignment

                ! Reorder Point based on the identified clusters and recompute the cluster properties.

                call Kmeans%getProp(nd, np, Point, PointIndex)
                if (Kmeans%Err%occurred) then ! This must revert back to the original Kmeans cluster and continue instead of return
                    neopt = 1_IK
                    PartitionSize(1) = np
                    PartitionMembership = 1_IK
                    return
                end if

                cycle loopRecursivePartitionOptimization

            end if blockReclusteringNeeded

            exit loopRecursivePartitionOptimization ! mission accomplished, current volumes minimized.

        end do loopRecursivePartitionOptimization

#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Rescale the volumes, if needed and for as long as needed, based on the user-input volume estimate
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (parLogVolNormedIsPresent) then
            do ic = 1, nc
#if defined DEBUG_ENABLED
                if (Kmeans%Prop%EffectiveSize(ic) > np) then
                    write(*,*) "FATAL: Internal error occurred: Kmeans%EffectiveSize(ic) > np", Kmeans%Prop%EffectiveSize(ic), np
                    error stop
                end if
#endif
                KmeansLogVolEstimate(ic) = parLogVolNormed + log( Kmeans%Prop%EffectiveSize(ic) / real(np, kind = SPR) )
                if (Kmeans%Prop%LogVolNormed(ic) < KmeansLogVolEstimate(ic)) then
                    !kmeansScaleFactor(ic) = kmeansScaleFactor(ic) * exp( (KmeansLogVolEstimate(ic) - Kmeans%Prop%LogVolNormed(ic)) / nd )
                    !Kmeans%ScaleFactorSq(ic) = kmeansScaleFactor(ic)**2
                    Kmeans%Prop%ScaleFactorSq(ic) = Kmeans%Prop%ScaleFactorSq(ic) * exp( 2 * (KmeansLogVolEstimate(ic) - Kmeans%Prop%LogVolNormed(ic)) / nd )
                    Kmeans%Prop%LogVolNormed(ic) = KmeansLogVolEstimate(ic)
                end if
            end do
            boundedRegionIsTooLarge = PartitionLogVolNormed(1) > logTightness + parLogVolNormed
        else
            boundedRegionIsTooLarge = .false.
        end if

        ! @todo: There is room for improvement here. When inclusionFraction /= 0, KmeansLogVolEstimate could be recursively checked for potential rescaling.

        !blockScaleFactorRescale: do
        !    do ic = 1, nc
        !        KmeansLogVolEstimate(ic) = parLogVolNormed + real( log( real(Kmeans%Prop%EffectiveSize(ic), kind = SPR) / real(np, kind = SPR) ) , kind = RK )
        !    end do
        !    Kmeans%Prop%EffectiveSize(ic) = Kmeans%Size(ic) + nint( inclusionFraction * count(Kmeans%Prop%MahalSq(IpStart(icOther):IpEnd(icOther),ic)<Kmeans%Prop%ScaleFactorSq(ic)), IK )
        !end blockScaleFactorRescale

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Check if further clustering is warranted.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!write(*,*) logTightness, parLogVolNormed, PartitionLogVolNormed(1), getLogSumExp(2, Kmeans%Prop%LogVolNormed)
!write(*,"(*(g0,:,','))")
!write(*,"(*(g0,:,','))")
!write(*,"(*(g0,:,','))") "numRecursiveCall", numRecursiveCall

        furtherClusteringCheck: if  (logTightness + getLogSumExp(nc, Kmeans%Prop%LogVolNormed) < PartitionLogVolNormed(1) .or. boundedRegionIsTooLarge) then

            ! At least nc sub-clusters is better than one cluster. Now search for more sub-sub-clusters.

            !ndEffectivePlusOne = 1_IK + np / (nemax + 1_IK)
            !KmeansNemax(1:nc) = max( Kmeans%Size(1:nc) / ndEffectivePlusOne , 1_IK )
            !KmeansNemax(1:nc) = min( nc , max( Kmeans%Size(1:nc) / ndEffectivePlusOne , 1_IK ) )

            !KmeansNemax(1) = min( nemax , max( Kmeans%Size(1:nc) / minClusterSize , 1_IK )
            !nemaxRemained = nemax - nc + 1_IK
            nemaxRemained = nemax
            KmeansNeopt(0) = 0_IK
            icstart = 1_IK

            ! Search for more subclusters.

            loopSearchForChildren: do ic = 1, nc

                nemaxRemained = nemaxRemained - KmeansNeopt(ic-1)
                !ndEffectivePlusOne = 1_IK + np / (nemaxRemained + 1_IK)
                KmeansNemax(ic) = min(nemaxRemained, Kmeans%Size(ic) / minClusterSize)
                !KmeansNemax(ic) = max( Kmeans%Size(ic) / ndEffectivePlusOne , 1_IK )
                !KmeansNemax(ic) = min(nemaxRemained, Kmeans%Size(ic) / minClusterSize)
                !KmeansNemax(ic) = Kmeans%Size(ic) / minClusterSize

                ipstart = Kmeans%Prop%CumSumSize(ic-1) + 1
                ipend = Kmeans%Prop%CumSumSize(ic)

                icstart = icstart + KmeansNeopt(ic-1)
                icend = icstart + KmeansNemax(ic) - 1

!write(*,"(*(g0,:,','))") "icstart, icend", icstart, icend
!write(*,"(*(g0,:,','))") "KmeansNeopt(ic-1)", KmeansNeopt(ic-1)
!write(*,"(*(g0,:,','))") "nc, nemax, KmeansNemax(ic)", nc, nemax, KmeansNemax(ic)
!write(*,"(*(g0,:,','))") "Kmeans%Size(ic)", Kmeans%Size(ic)
!write(*,"(*(g0,:,','))") "nemaxRemained", nemaxRemained

                if (KmeansNemax(ic) > 1_IK) then ! .and. Kmeans%Prop%LogVolNormed(ic)>1.1_RK*KmeansLogVolEstimate(ic)) then

                    PartitionLogVolNormed(icstart) = Kmeans%Prop%LogVolNormed(ic)

                    if (parLogVolNormedIsPresent) parLogVolNormed = KmeansLogVolEstimate(ic)

                    call runPartitionKernel ( nemax = KmeansNemax(ic) & ! LCOV_EXCL_LINE
                                            , nd = nd & ! LCOV_EXCL_LINE
                                            , np = Kmeans%Size(ic) & ! LCOV_EXCL_LINE
                                            , nc = min(nc, KmeansNemax(ic)) & ! LCOV_EXCL_LINE
                                            , nt = nt & ! LCOV_EXCL_LINE
#if !defined KMEANS
                                            , partitionMaxAllowedRecursion = partitionMaxAllowedRecursion & ! LCOV_EXCL_LINE
#endif
#if !(defined KMEANS || defined MAXDEN || defined MINVOL)
                                            , mahalSqWeightExponent = mahalSqWeightExponent & ! LCOV_EXCL_LINE
#endif
                                            , inclusionFraction = inclusionFraction & ! LCOV_EXCL_LINE
                                            , logTightness = logTightness & ! LCOV_EXCL_LINE
                                            , Point = Point(1:nd,ipstart:ipend) & ! LCOV_EXCL_LINE
                                            , neopt = KmeansNeopt(ic) & ! LCOV_EXCL_LINE
                                            , PointIndex = PointIndex(ipstart:ipend) & ! LCOV_EXCL_LINE
                                            , PartitionSize = PartitionSize(icstart:icend) & ! LCOV_EXCL_LINE
                                            , PartitionCenter = PartitionCenter(1:nd,icstart:icend) & ! LCOV_EXCL_LINE
                                            , PartitionChoDia = PartitionChoDia(1:nd,icstart:icend) & ! LCOV_EXCL_LINE
                                            , PartitionInvCovMat = PartitionInvCovMat(1:nd,1:nd,icstart:icend) & ! LCOV_EXCL_LINE
                                            , PartitionMembership = PartitionMembership(ipstart:ipend) & ! LCOV_EXCL_LINE
                                            , PartitionLogVolNormed = PartitionLogVolNormed(icstart:icend) & ! LCOV_EXCL_LINE
                                            , PartitionChoLowCovUpp = PartitionChoLowCovUpp(1:nd,1:nd,icstart:icend) & ! LCOV_EXCL_LINE
                                            , convergenceFailureCount = convergenceFailureCount & ! LCOV_EXCL_LINE
                                            , numRecursiveCall = numRecursiveCall & ! LCOV_EXCL_LINE
                                            , parLogVolNormed = parLogVolNormed & ! LCOV_EXCL_LINE
                                            , partitionKmeansRelTol = partitionKmeansRelTol & ! LCOV_EXCL_LINE
                                            , partitionMaxAllowedKmeansFailure = partitionMaxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                            , partitionMaxAllowedKmeansRecursion = partitionMaxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                                            )
                    PartitionMembership(ipstart:ipend) = PartitionMembership(ipstart:ipend) + icstart - 1_IK

                elseif (KmeansNemax(ic) == 1_IK) then

                    KmeansNeopt(ic) = 1_IK
                    PartitionMembership(ipstart:ipend) = icstart

                else

                    KmeansNeopt(ic:nc) = 0_IK
                    exit loopSearchForChildren

                end if

                if (KmeansNeopt(ic) == 1_IK) then ! .or. KmeansNemax(ic) == 1_IK ) then ! There is only one cluster here

                    KmeansNeopt(ic) = 1_IK
                    PartitionSize(icstart) = Kmeans%Size(ic)
                    PartitionChoLowCovUpp(1:nd,1:nd,icstart) = Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic)

                    ! Rescale the Cholesky factor and its diagonals such that they describe the bounded shape of the cluster

                    kmeansScaleFactor = sqrt(Kmeans%Prop%ScaleFactorSq(ic))
                    scaleFactorSqInverse = 1._RK / Kmeans%Prop%ScaleFactorSq(ic)
                    do j = 1, nd
                        PartitionCenter(j,icstart) = Kmeans%Center(j,ic)
                        PartitionChoDia(j,icstart) = Kmeans%Prop%ChoDia(j,ic) * kmeansScaleFactor
                        PartitionInvCovMat(1:nd,j,icstart) = Kmeans%Prop%InvCovMat(1:nd,j,ic) * scaleFactorSqInverse
                        do i = j + 1, nd
                            PartitionChoLowCovUpp(i,j,icstart) = PartitionChoLowCovUpp(i,j,icstart) * kmeansScaleFactor
                        end do
                    end do

#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
                    block
                        integer :: ip
                        logical :: isInside
                        real(RK) :: mahalSq
                        real(RK), allocatable :: NormedPoint(:)
                        do ip = ipstart, ipend
                            NormedPoint = Point(:,ip) - PartitionCenter(:,icstart)
                            mahalSq = dot_product(NormedPoint,matmul(PartitionInvCovMat(:,:,icstart),NormedPoint))
                            isInside = mahalSq - 1._RK <= 1.e-6_RK
                            if (.not. isInside) then
                                write(*,"(*(g0.15,:,' '))") new_line("a"), "FATAL - POINT NOT INSIDE!, MAHAL = ", sqrt(mahalSq), new_line("a")
                                write(*,"(60(g0,:,','))") "numRecursiveCall", numRecursiveCall
                                write(*,"(60(g0,:,','))") "KmeansNemax(icstart), icstart", KmeansNemax(icstart), icstart
                                write(*,"(60(g0,:,','))") "ip, IpStart, IpEnd, size(PartitionMembership)", ip, ipstart, ipend, size(PartitionMembership(ipstart:ipend))
                                write(*,"(60(g0,:,','))") "PartitionMembership", PartitionMembership(ipstart:ipend)
                                write(*,"(60(g0,:,','))") "PartitionChoLowCovUpp", PartitionChoLowCovUpp(:,:,icstart)
                                write(*,"(60(g0,:,','))") "PartitionInvCovMat", PartitionInvCovMat(:,:,icstart)
                                write(*,"(60(g0,:,','))") "PartitionChoDia", PartitionChoDia(:,icstart)
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
                            , PartitionMembership = self%Membership & ! LCOV_EXCL_LINE
                            , PartitionChoLowCovUpp = self%ChoLowCovUpp & ! LCOV_EXCL_LINE
                            , PartitionLogVolNormed = self%LogVolNormed & ! LCOV_EXCL_LINE
                            , PartitionChoDia = self%ChoDia & ! LCOV_EXCL_LINE
                            , Point = Point & ! LCOV_EXCL_LINE
                            )
    end subroutine writePartitionOOP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine writePartition   ( fileUnit & ! LCOV_EXCL_LINE
                                , nd,np,nc & ! LCOV_EXCL_LINE
                                , PartitionSize & ! LCOV_EXCL_LINE
                                , PartitionCenter & ! LCOV_EXCL_LINE
                                , PartitionMembership & ! LCOV_EXCL_LINE
                                , PartitionChoLowCovUpp & ! LCOV_EXCL_LINE
                                , PartitionLogVolNormed & ! LCOV_EXCL_LINE
                                , PartitionChoDia & ! LCOV_EXCL_LINE
                                , Point & ! LCOV_EXCL_LINE
                                )
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK)     , intent(in)    :: fileUnit
        integer(IK)     , intent(in)    :: nd, np, nc
        integer(IK)     , intent(in)    :: PartitionSize(nc)
        real(RK)        , intent(in)    :: PartitionCenter(nd,nc)
        integer(IK)     , intent(in)    :: PartitionMembership(np)
        real(RK)        , intent(in)    :: PartitionLogVolNormed(nc)
        real(RK)        , intent(in)    :: PartitionChoLowCovUpp(nd,nd,nc)
        real(RK)        , intent(in)    :: PartitionChoDia(nd,nc)
        real(RK)        , intent(in)    :: Point(nd,np)
        character(*)    , parameter     :: fileFormat = "(*(g0.8,:,','))"
        integer(IK)                     :: i, j, ic
        
        write(fileUnit,"(*(g0,:,'"//new_line("a")//"'))") "nd", nd, "np", np, "nc", nc

        write(fileUnit,"(A)") "Size"
        write(fileUnit,fileFormat) PartitionSize

        write(fileUnit,"(A)") "Center"
        write(fileUnit,fileFormat) PartitionCenter

        write(fileUnit,"(A)") "LogVolume"
        write(fileUnit,fileFormat) PartitionLogVolNormed

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
