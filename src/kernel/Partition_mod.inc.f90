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
!> This module contains procedures and routines for finding the minimum-volume bounding ellipsoid partitioning
!> of a given set of points by minimizing the volumes of the ellipsoids and comparing them against their parents.
!> \author
!> Amir Shahmoradi

    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type

    implicit none

#if defined MINVOL
    character(*), parameter :: MODULE_NAME = "@PartitionMinVol_mod"
#elif defined MAXDEN
    character(*), parameter :: MODULE_NAME = "@PartitionMaxDen_mod"
#endif

    !> \brief
    !> The `Partition_type` class. Partitions an input array `Point(nd,np)` with `nd` attributes and `np` observations (points).
    type :: Partition_type
        integer(IK)                 :: nsim                         !< See the interface of [constructPartition()](@ref constructpartition).
        integer(IK)                 :: nd,np                        !< See the interface of [constructPartition()](@ref constructpartition).
        integer(IK)                 :: nc,nt                        !< See the interface of [constructPartition()](@ref constructpartition).
        integer(IK)                 :: nemax                        !< See the interface of [constructPartition()](@ref constructpartition).
        integer(IK)                 :: minSize                      !< See the interface of [constructPartition()](@ref constructpartition).
        integer(IK)                 :: maxAllowedKmeansFailure      !< See the interface of [constructPartition()](@ref constructpartition).
        integer(IK)                 :: maxAllowedKmeansRecursion    !< See the interface of [constructPartition()](@ref constructpartition).
        integer(IK)                 :: maxAllowedKvolumeRecursion   !< See the interface of [constructPartition()](@ref constructpartition).
        real(RK)                    :: parentLogVolNormed           !< See the interface of [constructPartition()](@ref constructpartition).
        real(RK)                    :: pointLogVolNormed            !< The estimated volume of a single point based on the input value of `parentLogVolNormed`.
        real(RK)                    :: inclusionFraction            !< See the interface of [constructPartition()](@ref constructpartition).
!#if defined MAXDEN
!        real(RK)                    :: poisModeLogPMF               !< The log-probability mass function of the Poisson distribution at mode given a log-density of `poisLogDen`.
!        real(RK)                    :: poisLogDen                   !< The inferred log(density) from `parentLogVolNormed` and `np`.
!        real(RK)                    :: poisMode                     !< The (lower) mode of the Poisson distribution with a log-mean of `poisLogDen`.
!#endif
        real(RK)                    :: logExpansion                 !< See the interface of [constructPartition()](@ref constructpartition).
        real(RK)                    :: kmeansRelTol                 !< See the interface of [constructPartition()](@ref constructpartition).
        logical                     :: stanEnabled                  !< See the interface of [constructPartition()](@ref constructpartition).
        integer(IK)                 :: neopt                        !< The predicted optimal number of clusters identified in the input data.
        integer(IK)                 :: numRecursiveCall             !< The total number of recursive calls to the partitioning algorithm.
        integer(IK)                 :: convergenceFailureCount      !< The number of times the partitioning algorithm has failed to converge.
        integer(IK) , allocatable   :: Membership(:)                !< An array of size `(np)` representing the bounding-ellipsoid membership IDs of the corresponding data points.
        integer(IK) , allocatable   :: PointIndex(:)                !< An array of size `(np)` of indices such that the input `Point(1:nd,Index(ip))` is saved at the output `Point(1:nd,ip)`.
        real(RK)    , allocatable   :: LogVolNormed(:)              !< An array of size `(nemax)` representing the log-volumes of the corresponding bounding ellipsoids.
        real(RK)    , allocatable   :: ChoLowCovUpp(:,:,:)          !< An array of size `(nd,nd,nemax)` representing the Cholesky lower triangle, diagonal, and covariance matrices of the bounding ellipsoids.
        real(RK)    , allocatable   :: InvCovMat(:,:,:)             !< An array of size `(nd,nd,nemax)` representing the full symmetric inverse covariance matrices of the bounding ellipsoids.
        real(RK)    , allocatable   :: ChoDia(:,:)                  !< An array of size `(nd,nemax)` representing the Cholesky lower triangle, diagonal, and covariance matrices of the bounding ellipsoids.
        real(RK)    , allocatable   :: Center(:,:)                  !< An array of size `(nd,nemax)` representing the centers of the bounding ellipsoids.
        integer(IK) , allocatable   :: Size(:)                      !< An array of size `(nemax)` representing the sizes of the corresponding bounding ellipsoids.
        type(Err_type)              :: Err                          !< An object of class [Err_type](@ref err_mod::err_type) containing information about error occurrence.
    contains
        procedure, pass             :: write => writePartitionMethod
        procedure, pass             :: run => runPartition
    end type Partition_type

    interface Partition_type
        module procedure :: constructPartition
    end interface Partition_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This function can serve as the constructor for `Partition_type` while also performing the partitioning for the first time.
    !> Construct a hierarchical partitioning of the input `Point` using a combination of Kmeans, Kvolume, and density maximization algorithms.
    !>
    !> \param[in]   Point                       :   An array of size `(nd = #dimensions, np = #observations)` containing the input set of points to partition.
    !> \param[in]   nc                          :   The number of clusters at each level of partitioning (**optional**, default = `2`).
    !> \param[in]   nt                          :   The number of Kmeans clustering tries before choosing the optimal clustering among all (**optional**, default = `1`).
    !> \param[in]   nsim                        :   The number of Monte Carlo simulation tries to estimate the sum of volumes of partitions excluding overlaps (**optional**, default = `0`).
    !> \param[in]   nemax                       :   The maximum number of partitions to be expected or identified in the input `Point` (**optional**, default = `max(nc, 1 + np / (nd + 1))`).
    !> \param[in]   minSize                     :   The minimum allowed size of a partition, that is, the minimum allowed number of points in a partition (**optional**, default = `0`).
    !> \param[in]   stanEnabled                 :   A logical value, indicating weather `Point` should be standardized prior to each Kmeans clustering (**optional**, default = `.true.`).
    !> \param[in]   trimEnabled                 :   A logical value, indicating weather all allocatable components have been trimmed to the minimum size (**optional**, default = `.false.`).
    !> \param[in]   kmeansRelTol                :   The relative tolerance below which the Kmeans clustering is assumed to have converged (**optional**, default = `1.e-4`).
    !> \param[in]   logExpansion                :   The logarithm of the factor by which the volumes of the children bounding ellipsoids are enlarged before comparing their sums to their parents (**optional**, default = `0.`).
    !> \param[in]   inclusionFraction           :   The fraction of non-member points that **are** inside the partition, to be used in estimating of the true volumes of the partitions (**optional**, default = `1.`).
    !> \param[in]   parentLogVolNormed          :   The logarithm of the volume of the input set of points normalized to the volume of the unit `nd`-ball (**optional**, default = `-infinity`).
    !> \param[in]   maxAllowedKmeansFailure     :   The maximum number of times Kmeans clustering is allowed to fail (**optional**, default = `10`).
    !> \param[in]   maxAllowedKmeansRecursion   :   The maximum number of times Kmeans clustering is allowed to refine the cluster centers via the algorithm of Lloyd (**optional**, default = `300`).
    !> \param[in]   maxAllowedKvolumeRecursion  :   The maximum number of times Kvolume algorithm is run to minimize the partition volumes at a given partitioning level (**optional**, default = `3`).
    !>
    !> \warning
    !> On output, `Point` will be reordered.
    !>
    !> \warning
    !> When `trimEnabled = .false.`, only the first `1:neopt` elements of all allocatable arrays contain meaningful information.
    !> As such, the rest of the trailing elements must be properly ignored when using the allocatable Partition components.
    !> For example, upon exit, only `Partition%LogVolNormed(1:Partition%neopt)` contain meaningful information.
    !>
    !> \warning
    !> When `trimEnabled = .true.`, all allocatable arrays will be trimmed from size `nemax` to the size `neopt`.
    !> In such a case, calling the method `Partition%tun()` after constructing the `Partition` object for the first time
    !> will likely lead to segmentation fault error.
    !>
    !> \remark
    !> When `nsim = 0`, the effective sum of volumes will be computed and used assuming that no partitions overlap.
    !> This is essentially equivalent to summing over all partition volumes.
    !>
    !> \remark
    !> If no input value if provided for `parentLogVolNormed`, both `parentLogVolNormed` and `pointLogVolNormed` will (must) be set to `NEGINF_RK`.
    !>
    !> \remark
    !> If `parentLogVolNormed` is provided as input, the algorithm will attempt to keep the volumes of all partitions above
    !> the limit of `logExpansion + pointLogVolNormed + log(np)` where `np` is the number of points in the given partition.
    !>
    !> \author
    !> Amir Shahmoradi, April 03, 2017, 2:16 AM, ICES, University of Texas at Austin
    function constructPartition ( Point & ! LCOV_EXCL_LINE
                                , nc,nt & ! LCOV_EXCL_LINE
                                , nsim & ! LCOV_EXCL_LINE
                                , nemax & ! LCOV_EXCL_LINE
                                , minSize & ! LCOV_EXCL_LINE
                                , stanEnabled & ! LCOV_EXCL_LINE
                                , trimEnabled & ! LCOV_EXCL_LINE
                                , kmeansRelTol & ! LCOV_EXCL_LINE
                                , logExpansion & ! LCOV_EXCL_LINE
                                , inclusionFraction & ! LCOV_EXCL_LINE
                                , parentLogVolNormed & ! LCOV_EXCL_LINE
                                , maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                , maxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                                , maxAllowedKvolumeRecursion & ! LCOV_EXCL_LINE
                                ) result(Partition)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructPartition
#endif
        use Constants_mod, only: IK, RK, NEGINF_RK
        !use Math_mod, only: getLogVolUnitBall

        implicit none

        real(RK)    , intent(inout)             :: Point(:,:)
        integer(IK) , intent(in)    , optional  :: nc,nt
        integer(IK) , intent(in)    , optional  :: nsim
        integer(IK) , intent(in)    , optional  :: nemax
        integer(IK) , intent(in)    , optional  :: minSize
        logical     , intent(in)    , optional  :: stanEnabled
        logical     , intent(in)    , optional  :: trimEnabled
        real(RK)    , intent(in)    , optional  :: kmeansRelTol
        real(RK)    , intent(in)    , optional  :: logExpansion
        real(RK)    , intent(in)    , optional  :: inclusionFraction
        real(RK)    , intent(in)    , optional  :: parentLogVolNormed
        integer(IK) , intent(in)    , optional  :: maxAllowedKmeansFailure
        integer(IK) , intent(in)    , optional  :: maxAllowedKmeansRecursion
        integer(IK) , intent(in)    , optional  :: maxAllowedKvolumeRecursion
        type(Partition_type)                    :: Partition

#if defined MAXDEN
        real(RK), parameter :: LOG_EXPANSION = 0._RK
#elif defined MINVOL
        real(RK), parameter :: LOG_EXPANSION = log(1.3_RK)
#endif

        Partition%nd = size(Point(:,1))
        Partition%np = size(Point(1,:))

        !Partition%negLogVolUnitBall = -getLogVolUnitBall(Partition%nd)

        Partition%nc = 2_IK; if (present(nc)) Partition%nc = max(Partition%nc, nc)
        Partition%nt = 1_IK; if (present(nt)) Partition%nt = max(Partition%nt, nt)
        Partition%nsim = 0_IK; if (present(nsim)) Partition%nsim = max(Partition%nsim, nsim)
        Partition%minSize = Partition%nd + 1_IK; if (present(minSize)) Partition%minSize = max(1_IK, minSize)
        Partition%nemax = Partition%np / Partition%minSize; if (present(nemax)) Partition%nemax = max(1_IK, nemax)
!write(*,*) "minSize, Partition%np, nemax, Partition%minSize, Partition%nemax", minSize, Partition%np, nemax, Partition%minSize, Partition%nemax
       !Partition%nemax = 1_IK + Partition%np / Partition%nc; if (present(nemax)) Partition%nemax = nemax; Partition%nemax = max(Partition%nc, Partition%nemax)
        Partition%maxAllowedKvolumeRecursion = 3_IK; if (present(maxAllowedKvolumeRecursion)) Partition%maxAllowedKvolumeRecursion = maxAllowedKvolumeRecursion
        Partition%maxAllowedKmeansRecursion = 300_IK; if (present(maxAllowedKmeansRecursion)) Partition%maxAllowedKmeansRecursion = maxAllowedKmeansRecursion
        Partition%maxAllowedKmeansFailure = 10_IK; if (present(maxAllowedKmeansFailure)) Partition%maxAllowedKmeansFailure = maxAllowedKmeansFailure
        Partition%inclusionFraction = 0._RK; if (present(inclusionFraction)) Partition%inclusionFraction = inclusionFraction
        Partition%logExpansion = LOG_EXPANSION; if (present(logExpansion)) Partition%logExpansion = logExpansion
        Partition%kmeansRelTol = 1.e-4_RK; if (present(kmeansRelTol)) Partition%kmeansRelTol = kmeansRelTol
        Partition%stanEnabled = .true.; if (present(stanEnabled)) Partition%stanEnabled = stanEnabled

        allocate( Partition%Size         (Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%Center       (Partition%nd,Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%ChoDia       (Partition%nd,Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%InvCovMat    (Partition%nd,Partition%nd,Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%ChoLowCovUpp (Partition%nd,Partition%nd,Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%LogVolNormed (Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%Membership   (Partition%np) & ! LCOV_EXCL_LINE
                , Partition%PointIndex   (Partition%np) & ! LCOV_EXCL_LINE
                )

        call Partition%run(Point, parentLogVolNormed)

        if (present(trimEnabled)) then
            if (trimEnabled) then
                Partition%Size           = Partition%Size         (1:Partition%neopt)
                Partition%Center         = Partition%Center       (1:Partition%nd,1:Partition%neopt)
                Partition%ChoDia         = Partition%ChoDia       (1:Partition%nd,1:Partition%neopt)
                Partition%InvCovMat      = Partition%InvCovMat    (1:Partition%nd,1:Partition%nd,1:Partition%neopt)
                Partition%ChoLowCovUpp   = Partition%ChoLowCovUpp (1:Partition%nd,1:Partition%nd,1:Partition%neopt)
                Partition%LogVolNormed   = Partition%LogVolNormed (1:Partition%neopt)
            end if
        end if

    end function constructPartition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a method of the class [Partition_type](@ref partition_type).
    !> Perform recursive clustering of the input `Point(nd,np)`.
    !>
    !> \param[inout]    Partition   :   An object of class [Partition_type](@ref partition_type) that contains the partitioning configuration on input and the partition on output.
    !> \param[inout]    Point       :   An array of size `(nd,np)` containing the input set of points to partition.
    !>
    !> \warning
    !> On input, `nemax` must be always larger than or equal to `nc`.
    !>
    !> \warning
    !> On output, `Point` will be reordered.
    subroutine runPartition(Partition, Point, parentLogVolNormed)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: runPartition
#endif
        use Constants_mod, only: IK, RK, NEGINF_RK
        use Matrix_mod, only: getInvMatFromCholFac
        use Matrix_mod, only: getCholeskyFactor
        use String_mod, only: num2str
        use Sort_mod, only: indexArray

        implicit none

        class(Partition_type)   , intent(inout) :: Partition
        real(RK)                , intent(inout) :: Point(Partition%nd,Partition%np) ! This must appear after Partition declaration.
#if defined MAXDEN
        real(RK)                , intent(in)    :: parentLogVolNormed
#elif defined MINVOL
        real(RK)    , optional  , intent(in)    :: parentLogVolNormed
#endif

        integer(IK)                             :: i, j, ip, ic
        integer(IK) , allocatable               :: PointIndexSorted(:)
        real(RK)                                :: NormedPoint(Partition%nd,Partition%np)
        logical                                 :: boundedRegionIsTooLarge
        real(RK)                                :: PointMean(Partition%nd)
        real(RK)                                :: PointStd(Partition%nd)
        real(RK)                                :: scaleFactorSqInverse
        real(RK)                                :: scaleFactorSq
        real(RK)                                :: scaleFactor
        real(RK)                                :: mahalSq
#if defined MAXDEN
        real(RK)                                :: unifrnd
        real(RK)                                :: logVolRatio
        real(RK)                                :: probPoisDiff
#endif
        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME//"@runPartition()"

        Partition%Err%occurred = .false.

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Initialize the center of the first cluster
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Partition%Center(1:Partition%nd,1) = 0._RK
        do ip = 1, Partition%np
          Partition%PointIndex(ip) = ip
          Partition%Center(1:Partition%nd,1) = Partition%Center(1:Partition%nd,1) + Point(1:Partition%nd,ip)
        end do
        Partition%Center(1:Partition%nd,1) = Partition%Center(1:Partition%nd,1) / Partition%np

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the covariance matrix and the Cholesky factorization of the first cluster.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Compute the covariance matrix upper.

        Partition%ChoLowCovUpp(1:Partition%nd,1:Partition%nd,1) = 0._RK
        do ip = 1, Partition%np
            do j = 1, Partition%nd
                NormedPoint(j,ip) = Point(j,ip) - Partition%Center(j,1)
                Partition%ChoLowCovUpp(1:j,j,1) = Partition%ChoLowCovUpp(1:j,j,1) + NormedPoint(1:j,ip) * NormedPoint(j,ip)
            end do
        end do
        Partition%ChoLowCovUpp(1:Partition%nd,1:Partition%nd,1) = Partition%ChoLowCovUpp(1:Partition%nd,1:Partition%nd,1) / (Partition%np - 1)

        ! Compute Std if requested.

        if (Partition%stanEnabled) then
            do concurrent(i = 1:Partition%nd)
                PointStd(i) = sqrt(Partition%ChoLowCovUpp(i,i,1))
            end do
        end if

        ! Compute the Cholesky factor lower.

        call getCholeskyFactor(Partition%nd, Partition%ChoLowCovUpp(1:Partition%nd,1:Partition%nd,1), Partition%ChoDia(1:Partition%nd,1))
        if (Partition%ChoDia(1,1)<0._RK) then
            ! LCOV_EXCL_START
            Partition%Err%msg = PROCEDURE_NAME//": Singular covariance matrix detected for nemax, nd, np: "//num2str(Partition%nemax)//", "//num2str(Partition%nd)//", "//num2str(Partition%np)
            write(*,*) Partition%Err%msg ! xxx This needs to be fixed.
            error stop
            ! LCOV_EXCL_STOP
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the inverse covariance matrix.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Partition%InvCovMat(1:Partition%nd,1:Partition%nd,1) = getInvMatFromCholFac ( Partition%nd & ! LCOV_EXCL_LINE
                                                                                    , Partition%ChoLowCovUpp(1:Partition%nd,1:Partition%nd,1) & ! LCOV_EXCL_LINE
                                                                                    , Partition%ChoDia(1:Partition%nd,1) & ! LCOV_EXCL_LINE
                                                                                    )

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the Mahal distances squared and the scale factor which makes the covariance matrix a bounding ellipsoid.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        scaleFactorSq = NEGINF_RK
        do ip = 1, Partition%np
            mahalSq = dot_product( NormedPoint(1:Partition%nd,ip) , matmul(Partition%InvCovMat(1:Partition%nd,1:Partition%nd,1), NormedPoint(1:Partition%nd,ip)) )
            if (scaleFactorSq < mahalSq) scaleFactorSq = mahalSq
        end do
        scaleFactor = sqrt(scaleFactorSq)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the largest possible LogVolume of the first cluster based on the parent cluster and the input estimated volume.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ! Compute the volume of the bounded cluster.

        Partition%LogVolNormed(1) = sum( log( Partition%ChoDia(1:Partition%nd,1) ) ) + Partition%nd * log(scaleFactor)

        ! Rescale the scaleFactor, if needed, to enlarge the bounded cluster in the future.
        ! This must be done here, although it will be used at the very end.
        ! Otherwise, Partition%LogVolNormed(1) will be potentially overwritten.

#if defined MAXDEN
        Partition%parentLogVolNormed = Partition%logExpansion + parentLogVolNormed
        logVolRatio = Partition%LogVolNormed(1) - Partition%parentLogVolNormed
        Partition%pointLogVolNormed = Partition%parentLogVolNormed - log(real(Partition%np,RK))
        probPoisDiff = getProbPoisDiff(Partition%np, logVolRatio)
        call random_number(unifrnd)
        if (probPoisDiff < unifrnd) then
            boundedRegionIsTooLarge = logVolRatio > 0._RK
            if (.not. boundedRegionIsTooLarge) then
                scaleFactor = scaleFactor * exp( (Partition%parentLogVolNormed - Partition%LogVolNormed(1)) / Partition%nd )
                Partition%LogVolNormed(1) = Partition%parentLogVolNormed
                scaleFactorSq = scaleFactor**2
            end if
        else
            boundedRegionIsTooLarge = .false.
        end if
#elif defined MINVOL
        if (present(parentLogVolNormed)) then
            Partition%parentLogVolNormed = parentLogVolNormed
            Partition%pointLogVolNormed = Partition%parentLogVolNormed - log(real(Partition%np,RK))
            if (Partition%LogVolNormed(1) > Partition%parentLogVolNormed) then
                boundedRegionIsTooLarge = Partition%LogVolNormed(1) > Partition%logExpansion + Partition%parentLogVolNormed
            else
                boundedRegionIsTooLarge = .false.
                scaleFactor = scaleFactor * exp( (Partition%parentLogVolNormed - Partition%LogVolNormed(1)) / Partition%nd )
                Partition%LogVolNormed(1) = Partition%parentLogVolNormed
                scaleFactorSq = scaleFactor**2
            end if
        else
            Partition%parentLogVolNormed = NEGINF_RK
            Partition%pointLogVolNormed = NEGINF_RK
            boundedRegionIsTooLarge = .true.
        end if
#endif

        blockSubclusterSearch: if (Partition%nemax > 1_IK .and. boundedRegionIsTooLarge) then

            Partition%numRecursiveCall = 0_IK
            Partition%convergenceFailureCount = 0_IK

            blockStandardization: if (Partition%stanEnabled) then

                PointMean(1:Partition%nd) = Partition%Center(1:Partition%nd,1)
                call runRecursivePartition  ( nd    = Partition%nd & ! LCOV_EXCL_LINE
                                            , np    = Partition%np & ! LCOV_EXCL_LINE
                                            , nc    = min(Partition%nc, Partition%nemax) & ! LCOV_EXCL_LINE
                                            , nt    = Partition%nt & ! LCOV_EXCL_LINE
                                            , nsim  = Partition%nsim & ! LCOV_EXCL_LINE
                                            , nemax = Partition%nemax & ! LCOV_EXCL_LINE
                                            , minSize = Partition%minSize & ! LCOV_EXCL_LINE
                                            , kmeansRelTol = Partition%kmeansRelTol & ! LCOV_EXCL_LINE
                                            , logExpansion = Partition%logExpansion & ! LCOV_EXCL_LINE
                                            , inclusionFraction = Partition%inclusionFraction & ! LCOV_EXCL_LINE
                                            , pointLogVolNormed = Partition%pointLogVolNormed & ! LCOV_EXCL_LINE
                                            , parentLogVolNormed = Partition%parentLogVolNormed & ! LCOV_EXCL_LINE
                                            , maxAllowedKmeansFailure = Partition%maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                            , maxAllowedKmeansRecursion = Partition%maxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                                            , maxAllowedKvolumeRecursion = Partition%maxAllowedKvolumeRecursion & ! LCOV_EXCL_LINE
                                            , Point = NormedPoint & ! LCOV_EXCL_LINE
                                            , Size = Partition%Size & ! LCOV_EXCL_LINE
                                            , neopt = Partition%neopt & ! LCOV_EXCL_LINE
                                            , Center = Partition%Center & ! LCOV_EXCL_LINE
                                            , ChoDia = Partition%ChoDia & ! LCOV_EXCL_LINE
                                            , InvCovMat = Partition%InvCovMat & ! LCOV_EXCL_LINE
                                            !, Membership = NormedPointMembership & ! LCOV_EXCL_LINE
                                            , Membership = Partition%Membership & ! LCOV_EXCL_LINE
                                            , PointIndex = Partition%PointIndex & ! LCOV_EXCL_LINE
                                            , ChoLowCovUpp = Partition%ChoLowCovUpp & ! LCOV_EXCL_LINE
                                            , LogVolNormed = Partition%LogVolNormed & ! LCOV_EXCL_LINE
                                            , numRecursiveCall = Partition%numRecursiveCall & ! LCOV_EXCL_LINE
                                            , convergenceFailureCount = Partition%convergenceFailureCount & ! LCOV_EXCL_LINE
                                            , PointStd = PointStd & ! LCOV_EXCL_LINE
                                            )
                if (Partition%neopt > 1_IK) then
                    do concurrent(ic = 1:Partition%neopt)
                        Partition%Center(1:Partition%nd,ic) = Partition%Center(1:Partition%nd,ic) + PointMean(1:Partition%nd)
                    end do
                end if

                !block
                !use Timer_mod, only: Timer_type
                !type(Timer_type) :: Timer
                !call Timer%tic()
                ! The following step takes about 1 microseconds (in low dims) to 1 ms (in high dims).
                ! Sorting PointIndex at dimensions > 20 with 500-2000 tends to be superior to sorting Point.
                if (Partition%nd < 20_IK) then
                    NormedPoint(1:Partition%nd, 1:Partition%np) = Point(1:Partition%nd, Partition%PointIndex(1:Partition%np))
                    Point(1:Partition%nd, 1:Partition%np) = NormedPoint(1:Partition%nd, 1:Partition%np)
                else
                    if (.not. allocated(PointIndexSorted)) allocate(PointIndexSorted(Partition%np))
                    call indexArray(Partition%np, Partition%PointIndex, PointIndexSorted, Partition%Err)
                    if (Partition%Err%occurred) return ! LCOV_EXCL_LINE
                    Partition%Membership(1:Partition%np) = Partition%Membership(PointIndexSorted)
                end if
                !call Timer%toc()
                !write(*,*) "delta = ", Timer%Time%delta
                !end block

            else blockStandardization

                call runRecursivePartition  ( nd    = Partition%nd & ! LCOV_EXCL_LINE
                                            , np    = Partition%np & ! LCOV_EXCL_LINE
                                            , nc    = min(Partition%nc, Partition%nemax) & ! LCOV_EXCL_LINE
                                            , nt    = Partition%nt & ! LCOV_EXCL_LINE
                                            , nsim  = Partition%nsim & ! LCOV_EXCL_LINE
                                            , nemax = Partition%nemax & ! LCOV_EXCL_LINE
                                            , minSize = Partition%minSize & ! LCOV_EXCL_LINE
                                            , kmeansRelTol = Partition%kmeansRelTol & ! LCOV_EXCL_LINE
                                            , logExpansion = Partition%logExpansion & ! LCOV_EXCL_LINE
                                            , inclusionFraction = Partition%inclusionFraction & ! LCOV_EXCL_LINE
                                            , pointLogVolNormed = Partition%pointLogVolNormed & ! LCOV_EXCL_LINE
                                            , parentLogVolNormed = Partition%parentLogVolNormed & ! LCOV_EXCL_LINE
                                            , maxAllowedKmeansFailure = Partition%maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                            , maxAllowedKmeansRecursion = Partition%maxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                                            , maxAllowedKvolumeRecursion = Partition%maxAllowedKvolumeRecursion & ! LCOV_EXCL_LINE
                                            , Point = Point & ! LCOV_EXCL_LINE
                                            , Size = Partition%Size & ! LCOV_EXCL_LINE
                                            , neopt = Partition%neopt & ! LCOV_EXCL_LINE
                                            , Center = Partition%Center & ! LCOV_EXCL_LINE
                                            , ChoDia = Partition%ChoDia & ! LCOV_EXCL_LINE
                                            , InvCovMat = Partition%InvCovMat & ! LCOV_EXCL_LINE
                                            , Membership = Partition%Membership & ! LCOV_EXCL_LINE
                                            , PointIndex = Partition%PointIndex & ! LCOV_EXCL_LINE
                                            , ChoLowCovUpp = Partition%ChoLowCovUpp & ! LCOV_EXCL_LINE
                                            , LogVolNormed = Partition%LogVolNormed & ! LCOV_EXCL_LINE
                                            , numRecursiveCall = Partition%numRecursiveCall & ! LCOV_EXCL_LINE
                                            , convergenceFailureCount = Partition%convergenceFailureCount & ! LCOV_EXCL_LINE
                                            )

            end if blockStandardization

        else blockSubclusterSearch ! There is only one cluster

            Partition%neopt = 1_IK
            Partition%Size(1) = Partition%np ! @todo : xxx this is redundant and should be removed in the future.
            Partition%Membership(1:Partition%np) = 1_IK ! @todo : xxx this is redundant and should be removed in the future.

        end if blockSubclusterSearch

        ! Set the single cluster properties.

        if (Partition%neopt == 1_IK) then

            ! Now rescale Cholesky factor and its diagonals, such that they describe the bounding ellipsoid of the cluster.

            scaleFactorSqInverse = 1._RK / scaleFactorSq
            do j = 1, Partition%nd
                Partition%ChoDia(j,1) = Partition%ChoDia(j,1) * scaleFactor
                !Partition%ChoLowCovUpp(1:j,j,1) = Partition%ChoLowCovUpp(1:j,j,1) * scaleFactorSq
                Partition%InvCovMat(1:Partition%nd,j,1) = Partition%InvCovMat(1:Partition%nd,j,1) * scaleFactorSqInverse
                do i = j+1, Partition%nd
                    Partition%ChoLowCovUpp(i,j,1) = Partition%ChoLowCovUpp(i,j,1) * scaleFactor
                end do
            end do

#if defined DEBUG_ENABLED || defined TESTING_ENABLED || defined CODECOVE_ENABLED
        elseif (Partition%neopt < 1_IK) then

            write(*,*) PROCEDURE_NAME//": Internal error occurred : Partition%neopt < 1_IK : ", Partition%neopt
            error stop
#endif

        end if

    end subroutine runPartition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    ! Recursively partition the input `Point`.
    !>
    !> \warning
    !> if `parentLogVolNormed > NEGINF_RK` , then `pointLogVolNormed > NEGINF_RK` must also hold.
    recursive subroutine runRecursivePartition  ( nd & ! LCOV_EXCL_LINE
                                                , np & ! LCOV_EXCL_LINE
                                                , nc & ! LCOV_EXCL_LINE
                                                , nt & ! LCOV_EXCL_LINE
                                                , nsim & ! LCOV_EXCL_LINE
                                                , nemax & ! LCOV_EXCL_LINE
                                                , minSize & ! LCOV_EXCL_LINE
                                                , kmeansRelTol & ! LCOV_EXCL_LINE
                                                , logExpansion & ! LCOV_EXCL_LINE
                                                , inclusionFraction & ! LCOV_EXCL_LINE
                                                , pointLogVolNormed & ! LCOV_EXCL_LINE
                                                , parentLogVolNormed & ! LCOV_EXCL_LINE
                                                , maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                                , maxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                                                , maxAllowedKvolumeRecursion & ! LCOV_EXCL_LINE
                                                , Point & ! LCOV_EXCL_LINE
                                                , Size & ! LCOV_EXCL_LINE
                                                , neopt & ! LCOV_EXCL_LINE
                                                , Center & ! LCOV_EXCL_LINE
                                                , ChoDia & ! LCOV_EXCL_LINE
                                                , InvCovMat & ! LCOV_EXCL_LINE
                                                , Membership & ! LCOV_EXCL_LINE
                                                , PointIndex & ! LCOV_EXCL_LINE
                                                , ChoLowCovUpp & ! LCOV_EXCL_LINE
                                                , LogVolNormed & ! LCOV_EXCL_LINE
                                                , numRecursiveCall & ! LCOV_EXCL_LINE
                                                , convergenceFailureCount & ! LCOV_EXCL_LINE
                                                , PointStd & ! LCOV_EXCL_LINE
                                                )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: runRecursivePartition
#endif
!#if defined DEBUG_ENABLED
!        use, intrinsic :: ieee_arithmetic, only: ieee_support_underflow_control
!        use, intrinsic :: ieee_arithmetic, only: ieee_set_underflow_mode
!#endif
        use Constants_mod, only: IK, RK, NEGINF_RK !, SPR
        use Kmeans_mod, only: Kmeans_type !, getKmeans, runKmeans
        use Math_mod, only: getLogSumExp

        implicit none

        integer(IK) , intent(in)                :: nd, np, nc, nt, nsim, nemax, minSize
        real(RK)    , intent(in)                :: logExpansion
        real(RK)    , intent(in)                :: kmeansRelTol
        real(RK)    , intent(in)                :: inclusionFraction
        real(RK)    , intent(in)                :: pointLogVolNormed
        real(RK)    , intent(in)                :: parentLogVolNormed
        integer(IK) , intent(in)                :: maxAllowedKmeansFailure
        integer(IK) , intent(in)                :: maxAllowedKmeansRecursion
        integer(IK) , intent(in)                :: maxAllowedKvolumeRecursion
        integer(IK) , intent(out)               :: neopt
        integer(IK) , intent(out)               :: Size(nemax)
        integer(IK) , intent(inout)             :: numRecursiveCall
        integer(IK) , intent(inout)             :: convergenceFailureCount
        integer(IK) , intent(inout)             :: Membership(np)
        integer(IK) , intent(inout)             :: PointIndex(np)
        real(RK)    , intent(inout)             :: Point(nd,np)
        real(RK)    , intent(inout)             :: Center(nd,nemax)
        real(RK)    , intent(inout)             :: ChoLowCovUpp(nd,nd,nemax)
        real(RK)    , intent(inout)             :: InvCovMat(nd,nd,nemax)
        real(RK)    , intent(inout)             :: LogVolNormed(nemax)
        real(RK)    , intent(inout)             :: ChoDia(nd,nemax)
        real(RK)    , intent(inout), optional   :: PointStd(nd)

        type(Kmeans_type)                       :: Kmeans
        real(RK)    , allocatable               :: InvStd(:)
        real(RK)                                :: ndInverse
        real(RK)                                :: StanPoint(nd,np)
        real(RK)                                :: KmeansLogVolEstimate(nc) ! New Cluster Volume estimates.
        real(RK)                                :: scaleFactorSqInverse
        real(RK)                                :: scaleFactor
        real(RK)                                :: maxLogVolNormed
        integer(IK)                             :: KmeansNemax(nc)      ! maximum allowed number of clusters in each of the two K-means clusters.
        integer(IK)                             :: KmeansNeopt(0:nc)    ! optimal number of clusters in each of the child clusters.
        integer(IK)                             :: i, j, ip, ic, jc, ipstart, ipend
        integer(IK)                             :: icstart, icend, recursionCounter
        integer(IK)                             :: icmin, nemaxRemained
        logical                                 :: boundedRegionIsTooLarge
        logical                                 :: reclusteringNeeded
#if defined MAXDEN
        real(RK)                                :: unifrnd
        real(RK)                                :: logVolRatio
        real(RK)                                :: probPoisDiff
#endif

        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME//"@runRecursivePartition()"

        ndInverse = 1._RK / nd
        numRecursiveCall = numRecursiveCall + 1_IK

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Perform the initial Kmeans clustering.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        blockKmeansClustering: if (present(PointStd)) then

            InvStd = 1._RK / PointStd(1:nd)
            do concurrent(ip = 1:np)
                StanPoint(1:nd,ip) = Point(1:nd,ip) * InvStd(1:nd)
            end do

            Kmeans = Kmeans_type( nd = nd & ! LCOV_EXCL_LINE
                                , np = np & ! LCOV_EXCL_LINE
                                , nc = nc & ! LCOV_EXCL_LINE
                                , nt = nt & ! LCOV_EXCL_LINE
                                , Point = StanPoint & ! LCOV_EXCL_LINE
                                , relTol = kmeansRelTol & ! LCOV_EXCL_LINE
                                , minSize = minSize & ! LCOV_EXCL_LINE
                                , nfailMax = maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                , niterMax = maxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                                )
            if (Kmeans%Err%occurred) then
                Membership = 1
                Size(1) = np
                neopt = 1
                return
            end if

            ! Shift centers back.

            do concurrent(ic = 1:nc)
                Kmeans%Center(1:nd,ic) = Kmeans%Center(1:nd,ic) * PointStd(1:nd) !+ PointMean(1:nd)
            end do

        else blockKmeansClustering

            Kmeans = Kmeans_type( nd = nd & ! LCOV_EXCL_LINE
                                , np = np & ! LCOV_EXCL_LINE
                                , nc = nc & ! LCOV_EXCL_LINE
                                , nt = nt & ! LCOV_EXCL_LINE
                                , Point = Point & ! LCOV_EXCL_LINE
                                , relTol = kmeansRelTol & ! LCOV_EXCL_LINE
                                , minSize = minSize & ! LCOV_EXCL_LINE
                                , nfailMax = maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                , niterMax = maxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                                )
            if (Kmeans%Err%occurred) then
                Membership = 1
                Size(1) = np
                neopt = 1
                return
            end if

        end if blockKmeansClustering

        ! Compute cluster properties.

        call Kmeans%getProp ( nd = nd & ! LCOV_EXCL_LINE
                            , np = np & ! LCOV_EXCL_LINE
                            , Point = Point & ! LCOV_EXCL_LINE
                            , Index = PointIndex & ! LCOV_EXCL_LINE
                            , inclusionFraction = inclusionFraction & ! LCOV_EXCL_LINE
                            , pointLogVolNormed = pointLogVolNormed & ! LCOV_EXCL_LINE
                            )
        if (Kmeans%Err%occurred) then
            ! LCOV_EXCL_START
            Membership = 1
            Size(1) = np
            neopt = 1
            return
            ! LCOV_EXCL_STOP
        end if

#if false && (defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED)
        loopComputeScaleFactor: do ic = 1, nc
                block
                    use Statistics_mod, only: getSamCholFac, getMahalSq
                    real(RK), parameter :: TOLERANCE = 1.e-2_RK
                    real(RK) :: ChoLowDiff(nd,nd), ChoDiaDiff(nd)
                    real(RK) :: ChoLow(nd,nd), ChoDia(nd)
                    real(RK) :: NormedPoint(nd)
                    real(RK) :: dum, diff, scaleFac
                    call getSamCholFac  ( nd = nd & ! LCOV_EXCL_LINE
                                        , np = Kmeans%Size(ic) & ! LCOV_EXCL_LINE
                                        , Mean = Kmeans%Center(1:nd,ic) & ! LCOV_EXCL_LINE
                                        , Point = Point(1:nd,Kmeans%Prop%CumSumSize(ic-1)+1:Kmeans%Prop%CumSumSize(ic)) & ! LCOV_EXCL_LINE
                                        , CholeskyLower = ChoLow & ! LCOV_EXCL_LINE
                                        , CholeskyDiago = ChoDia & ! LCOV_EXCL_LINE
                                        )
                    scaleFac = sqrt(Kmeans%Prop%ScaleFactorSq(ic))
                    ChoDia = ChoDia * scaleFac
                    do j = 1, nd
                        do i = j+1, nd
                            ChoLow(i,j) = ChoLow(i,j) * scaleFac
                        end do
                    end do
                    ChoLowDiff = 2 * abs(Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic)-ChoLow) / abs(Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic)+ChoLow)
                    if ( any(ChoLowDiff > TOLERANCE) ) then
                        ! LCOV_EXCL_START
                        write(*,*)
                        write(*,*) PROCEDURE_NAME
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
                        error stop
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
                            write(*,*) PROCEDURE_NAME
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
        ! If optimization is enabled (i.e., not kmeans), reassign points to clusters based on minimum-volume, maximum-density, ...
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !KmeansBackup = Kmeans_type(nd, np, nc)
        !call KmeansBackup%Prop%allocate(nd, np, nc)
        !KmeansBackup%Prop%LogVolNormed = 2 * Kmeans%Prop%LogVolNormed

        loopRecursivePartitionOptimization: do recursionCounter = 1, maxAllowedKvolumeRecursion
!write(*,*) "We should not be here"; error stop ! xxx remove
            ! Ensure the stability of the minimum bounding regions by making a backup.

            !if (Kmeans%Prop%logSumVolNormed < KmeansBackup%Prop%logSumVolNormed) then ! Create a backup of the current clustering properties. This must never occur on the first try.
            !    !KmeansBackup = Kmeans
            !    BackupPoint(1:nd,1:np)                          = Point(1:nd,1:np)
            !    BackupPointIndex(1:np)                          = PointIndex(1:np)
            !    KmeansBackup%Membership(1:np)                   = Kmeans%Membership(1:np)
            !    KmeansBackup%Center(1:nd,1:nc)                  = Kmeans%Center(1:nd,1:nc)
            !    KmeansBackup%Prop%ChoDia(1:nd,1:nc)             = Kmeans%Prop%ChoDia(1:nd,1:nc)
            !    KmeansBackup%Prop%MahalSq(1:np,1:nc)            = Kmeans%Prop%MahalSq(1:np,1:nc)
            !    KmeansBackup%Prop%ScaleFactorSq(1:nc)           = Kmeans%Prop%ScaleFactorSq(1:nc)
            !    KmeansBackup%Prop%InvCovMat(1:nd,1:nd,1:nc)     = Kmeans%Prop%InvCovMat(1:nd,1:nd,1:nc)
            !    KmeansBackup%Prop%ChoLowCovUpp(1:nd,1:nd,1:nc)  = Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,1:nc)
            !    KmeansBackup%Prop%EffectiveSize(1:nc)           = Kmeans%Prop%EffectiveSize(1:nc)
            !    KmeansBackup%Prop%LogVolNormed(1:nc)            = Kmeans%Prop%LogVolNormed(1:nc)
            !    KmeansBackup%Prop%LogDenNormed(1:nc)            = Kmeans%Prop%LogDenNormed(1:nc)
            !    KmeansBackup%Prop%CumSumSize(1:nc)              = Kmeans%Prop%CumSumSize(1:nc)
            !    KmeansBackup%Prop%logSumVolNormed               = Kmeans%Prop%logSumVolNormed
            !else ! Restore the last clustering properties.
            !    !Kmeans = KmeansBackup
            !    Point(1:nd,1:np)                                = BackupPoint(1:nd,1:np)
            !    PointIndex(1:np)                                = BackupPointIndex(1:np)
            !    Kmeans%Membership(1:np)                         = KmeansBackup%Membership(1:np)
            !    Kmeans%Center(1:nd,1:nc)                        = KmeansBackup%Center(1:nd,1:nc)
            !    Kmeans%Prop%ChoDia(1:nd,1:nc)                   = KmeansBackup%Prop%ChoDia(1:nd,1:nc)
            !    Kmeans%Prop%MahalSq(1:np,1:nc)                  = KmeansBackup%Prop%MahalSq(1:np,1:nc)
            !    Kmeans%Prop%ScaleFactorSq(1:nc)                 = KmeansBackup%Prop%ScaleFactorSq(1:nc)
            !    Kmeans%Prop%InvCovMat(1:nd,1:nd,1:nc)           = KmeansBackup%Prop%InvCovMat(1:nd,1:nd,1:nc)
            !    Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,1:nc)        = KmeansBackup%Prop%ChoLowCovUpp(1:nd,1:nd,1:nc)
            !    Kmeans%Prop%EffectiveSize(1:nc)                 = KmeansBackup%Prop%EffectiveSize(1:nc)
            !    Kmeans%Prop%LogVolNormed(1:nc)                  = KmeansBackup%Prop%LogVolNormed(1:nc)
            !    Kmeans%Prop%LogDenNormed(1:nc)                  = KmeansBackup%Prop%LogDenNormed(1:nc)
            !    Kmeans%Prop%CumSumSize(1:nc)                    = KmeansBackup%Prop%CumSumSize(1:nc)
            !    Kmeans%Prop%logSumVolNormed                     = KmeansBackup%Prop%logSumVolNormed
            !    exit loopRecursivePartitionOptimization
            !end if

            ! Convert Mean(Point) to sum(Point) per cluster.

            do concurrent(ic = 1:nc)
                Kmeans%Center(1:nd,ic) = Kmeans%Center(1:nd,ic) * Kmeans%Size(ic)
            end do

            ! Optimize the Kmeans clusters.

            reclusteringNeeded = .false.
            loopReassignPoint: do ip = 1, np

                ! Determine the cluster to which the point is the closest.

                icmin = Kmeans%Membership(ip)
                do ic = 1, nc
                    if (ic/=icmin .and. Kmeans%Prop%MahalSq(ip,ic) < Kmeans%Prop%MahalSq(ip,icmin)) icmin = ic
                end do

                ! Switch the point membership if needed.

                if (icmin /= Kmeans%Membership(ip)) then ! .and. Kmeans%Size(Kmeans%Membership(ip)) > minSize) then
                !if (icmin /= Kmeans%Membership(ip) .and. Kmeans%Size(Kmeans%Membership(ip)) > minSize) then
                    Kmeans%Size(icmin) = Kmeans%Size(icmin) + 1_IK
                    Kmeans%Size(Kmeans%Membership(ip)) = Kmeans%Size(Kmeans%Membership(ip)) - 1_IK
                    Kmeans%Center(1:nd,icmin) = Kmeans%Center(1:nd,icmin) + Point(1:nd,ip)
                    Kmeans%Center(1:nd,Kmeans%Membership(ip)) = Kmeans%Center(1:nd,Kmeans%Membership(ip)) - Point(1:nd,ip)
                    Kmeans%Membership(ip) = icmin
                    reclusteringNeeded = .true.
                end if

            end do loopReassignPoint

#if defined DEBUG_ENABLED
            !if (any(Kmeans%Size<minSize)) then
            !    !write(*,*) "Kmeans%Size < minSize : minSize, Kmeans%Size = ", minSize, Kmeans%Size
            !    !error stop ! @attn: This is NOT an error anymore as Size can become less than minSize.
            !end if
#endif

            ! Convert sum(Point) back to center of clusters.

            !do concurrent(ic = 1:nc) ! xxx revive this line
            do ic = 1, nc
                !write(*,*) "Kmeans%Center(1:nd,ic), Kmeans%Size", Kmeans%Center(1:nd,ic), Kmeans%Size ! xxx
                if (Kmeans%Size(ic) > 0_IK) Kmeans%Center(1:nd,ic) = Kmeans%Center(1:nd,ic) / Kmeans%Size(ic)
            end do

            ! Restart the process if anything has changed.

            blockReclusteringNeeded: if (reclusteringNeeded) then ! perform reassignment

                ! Reorder Point based on the identified clusters and recompute the cluster properties.

                call Kmeans%getProp(nd = nd, np = np, Point = Point, Index = PointIndex, inclusionFraction = inclusionFraction, pointLogVolNormed = pointLogVolNormed)
                if (Kmeans%Err%occurred) then
                    ! LCOV_EXCL_START
                    Membership = 1
                    Size(1) = np
                    neopt = 1
                    return
                    ! LCOV_EXCL_STOP
                end if

                cycle loopRecursivePartitionOptimization

            end if blockReclusteringNeeded

            exit loopRecursivePartitionOptimization ! mission accomplished, current volumes minimized.

        end do loopRecursivePartitionOptimization

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Rescale the volumes if needed and for as long as needed, based on the user-input volume estimate
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined MINVOL
        maxLogVolNormed = NEGINF_RK ! indicator for rescaling.
        if (parentLogVolNormed > NEGINF_RK) then

            do ic = 1, nc

#if defined DEBUG_ENABLED
                if (Kmeans%Prop%EffectiveSize(ic) > np) then
                    write(*,*) PROCEDURE_NAME
                    write(*,*) "FATAL: Internal error occurred: Kmeans%EffectiveSize(ic) > np", Kmeans%Prop%EffectiveSize(ic), np
                    error stop
                end if
#endif
                !KmeansLogVolEstimate(ic) = pointLogVolNormed + log(real(Kmeans%Prop%EffectiveSize(ic),RK)) ! @attn: xxx how is Kmeans%Prop%EffectiveSize(ic) > 0 ensured?
                if (Kmeans%Size(ic) > 0_IK) then

                    KmeansLogVolEstimate(ic) = pointLogVolNormed + log(real(Kmeans%Prop%EffectiveSize(ic),RK)) ! @attn: Kmeans%Prop%EffectiveSize(ic) >= Kmeans%Size(ic) > 0

                    if (Kmeans%Prop%LogVolNormed(ic) < KmeansLogVolEstimate(ic)) then

                        Kmeans%Prop%LogVolNormed(ic) = KmeansLogVolEstimate(ic)
                        if (Kmeans%Prop%LogVolNormed(ic) > maxLogVolNormed) maxLogVolNormed = Kmeans%Prop%LogVolNormed(ic)

                        ! Rescale the Cholesky factor and its diagonals such that they describe the bounded shape of the cluster
                        ! @todo: xxx This will have to be moved to a later location for better efficiency.
                        ! pay attention to the correct value of Kmeans%Prop%ScaleFactorSq(ic).

                        scaleFactor = exp( ndInverse * (KmeansLogVolEstimate(ic) - Kmeans%Prop%LogVolNormed(ic)) )
                        Kmeans%Prop%ScaleFactorSq(ic) = scaleFactor**2
                        scaleFactorSqInverse = 1._RK / Kmeans%Prop%ScaleFactorSq(ic)
                        do j = 1, nd
                            Kmeans%Prop%ChoDia(j,ic) = Kmeans%Prop%ChoDia(j,ic) * scaleFactor
                            Kmeans%Prop%InvCovMat(1:nd,j,ic) = Kmeans%Prop%InvCovMat(1:nd,j,ic) * scaleFactorSqInverse
                            do i = j + 1, nd
                                !write(*,*) "ChoLowCovUpp(i,j,ic), scaleFactor, Kmeans%Size(ic)", ChoLowCovUpp(i,j,ic), scaleFactor, Kmeans%Size(ic) ! xxx
                                Kmeans%Prop%ChoLowCovUpp(i,j,ic) = Kmeans%Prop%ChoLowCovUpp(i,j,ic) * scaleFactor
                            end do
                        end do

                    end if

#if defined DEBUG_ENABLED
                elseif (Kmeans%Size(ic) == 0_IK .and. Kmeans%Prop%LogVolNormed(ic) /= NEGINF_RK) then
                    write(*,*) PROCEDURE_NAME//": Kmeans%Size(ic) == 0_IK .and. Kmeans%Prop%LogVolNormed(ic) /= NEGINF_RK", ic, Kmeans%Prop%LogVolNormed(ic)
                    error stop
#endif

                end if

            end do
            boundedRegionIsTooLarge = LogVolNormed(1) > logExpansion + parentLogVolNormed

        else

            KmeansLogVolEstimate(1:nc) = NEGINF_RK
            boundedRegionIsTooLarge = .false.

        end if

        if (maxLogVolNormed > NEGINF_RK) Kmeans%Prop%logSumVolNormed = getLogSumExp(nc, Kmeans%Prop%LogVolNormed, maxLogVolNormed)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Check if further clustering is warranted.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !blockSubclusterSearch: if  (logExpansion + getLogSumExp(nc, Kmeans%Prop%LogVolNormed) < LogVolNormed(1) .or. boundedRegionIsTooLarge) then
        !blockSubclusterSearch: if  (getLogSumExp(nc, Kmeans%Prop%LogVolNormed) < LogVolNormed(1) .or. boundedRegionIsTooLarge) then
        blockSubclusterSearch: if  (boundedRegionIsTooLarge .or. LogVolNormed(1) > logExpansion + Kmeans%Prop%logSumVolNormed) then
#endif

            ! At least nc sub-clusters is better than one cluster. Now search for more sub-sub-clusters.

            nemaxRemained = nemax
            KmeansNeopt(0) = 0_IK
            icstart = 1_IK

            ! Search for more subclusters.

            loopSubclusterSearch: do ic = 1, nc

                nemaxRemained = nemaxRemained - KmeansNeopt(ic-1)
                KmeansNemax(ic) = min(nemaxRemained - nc + ic, Kmeans%Size(ic) / minSize)
                boundedRegionIsTooLarge = KmeansNemax(ic) > 1_IK

#if defined MAXDEN
                if (Kmeans%size(ic) > 0_IK) then
                    KmeansLogVolEstimate(ic) = pointLogVolNormed + log(real(Kmeans%Prop%EffectiveSize(ic),RK)) ! @attn: Kmeans%Prop%EffectiveSize(ic) >= Kmeans%Size(ic) > 0
                    logVolRatio = Kmeans%Prop%LogVolNormed(ic) - KmeansLogVolEstimate(ic)
                    probPoisDiff = getProbPoisDiff(Kmeans%Prop%EffectiveSize(ic), logVolRatio)
                    call random_number(unifrnd)
                    boundedRegionIsTooLarge = boundedRegionIsTooLarge .and. probPoisDiff < unifrnd .and. logVolRatio > 0._RK
                end if
#endif

#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
                if (nemaxRemained < KmeansNemax(ic)) then
                    write(*,"(*(g0,:,', '))") PROCEDURE_NAME//": nemaxRemained < KmeansNemax(ic) : numRecursiveCall, np, ic, minSize, nemaxRemained, KmeansNemax(ic), nemax" &
                                            , numRecursiveCall, np, ic, minSize, nemaxRemained, KmeansNemax(ic), nemax
                    error stop
                end if
#endif

                ipstart = Kmeans%Prop%CumSumSize(ic-1) + 1
                ipend = Kmeans%Prop%CumSumSize(ic)

                icstart = icstart + KmeansNeopt(ic-1)
                icend = icstart + KmeansNemax(ic) - 1

                blockSubclusterCheck: if (boundedRegionIsTooLarge) then ! .and. Kmeans%Prop%LogVolNormed(ic)>1.1_RK*KmeansLogVolEstimate(ic)) then

                    LogVolNormed(icstart) = Kmeans%Prop%LogVolNormed(ic)

                    if (present(PointStd)) then

                        do concurrent(i = 1:nd)
                            PointStd(i) = sqrt(Kmeans%Prop%ChoLowCovUpp(i,i,ic))
                        end do

                        do concurrent(ip = ipstart:ipend)
                            Point(1:nd,ip) = Point(1:nd,ip) - Kmeans%Center(1:nd,ic)
                        end do

                    end if

                    call runRecursivePartition  ( nd = nd & ! LCOV_EXCL_LINE
                                                , np = Kmeans%Size(ic) & ! LCOV_EXCL_LINE
                                                , nc = min(nc, KmeansNemax(ic)) & ! LCOV_EXCL_LINE
                                                , nt = nt & ! LCOV_EXCL_LINE
                                                , nsim = nsim & ! LCOV_EXCL_LINE
                                                , nemax = KmeansNemax(ic) & ! LCOV_EXCL_LINE
                                                , minSize = minSize & ! LCOV_EXCL_LINE
                                                , PointStd = PointStd & ! LCOV_EXCL_LINE
                                                , kmeansRelTol = kmeansRelTol & ! LCOV_EXCL_LINE
                                                , logExpansion = logExpansion & ! LCOV_EXCL_LINE
                                                , inclusionFraction = inclusionFraction & ! LCOV_EXCL_LINE
                                                , pointLogVolNormed = pointLogVolNormed & ! LCOV_EXCL_LINE
                                                , parentLogVolNormed = KmeansLogVolEstimate(ic) & ! LCOV_EXCL_LINE
                                                , maxAllowedKmeansFailure = maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                                , maxAllowedKmeansRecursion = maxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                                                , maxAllowedKvolumeRecursion = maxAllowedKvolumeRecursion & ! LCOV_EXCL_LINE
                                                , Point = Point(1:nd,ipstart:ipend) & ! LCOV_EXCL_LINE
                                                , Size = Size(icstart:icend) & ! LCOV_EXCL_LINE
                                                , neopt = KmeansNeopt(ic) & ! LCOV_EXCL_LINE
                                                , Center = Center(1:nd,icstart:icend) & ! LCOV_EXCL_LINE
                                                , ChoDia = ChoDia(1:nd,icstart:icend) & ! LCOV_EXCL_LINE
                                                , InvCovMat = InvCovMat(1:nd,1:nd,icstart:icend) & ! LCOV_EXCL_LINE
                                                , Membership = Membership(ipstart:ipend) & ! LCOV_EXCL_LINE
                                                , PointIndex = PointIndex(ipstart:ipend) & ! LCOV_EXCL_LINE
                                                , ChoLowCovUpp = ChoLowCovUpp(1:nd,1:nd,icstart:icend) & ! LCOV_EXCL_LINE
                                                , LogVolNormed = LogVolNormed(icstart:icend) & ! LCOV_EXCL_LINE
                                                , numRecursiveCall = numRecursiveCall & ! LCOV_EXCL_LINE
                                                , convergenceFailureCount = convergenceFailureCount & ! LCOV_EXCL_LINE
                                                )

                    if (present(PointStd) .and. KmeansNeopt(ic) > 1_IK) then
                        do concurrent(jc = icstart : icstart + KmeansNeopt(ic) - 1)
                            Center(1:nd,jc) = Center(1:nd,jc) + Kmeans%Center(1:nd,ic)
                        end do
                    end if

#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
                    if (present(PointStd)) then
                        do concurrent(ip = ipstart:ipend)
                            Point(1:nd,ip) = Point(1:nd,ip) + Kmeans%Center(1:nd,ic)
                        end do
                    end if

                    if (any(Membership(ipstart:ipend) > KmeansNemax(ic))) then
                        write(*,*) PROCEDURE_NAME
                        write(*,*) "KmeansNemax(ic) : ", KmeansNemax(ic)
                        write(*,*) "Membership(ipstart:ipend) > KmeansNemax(ic) : ", Membership(ipstart:ipend)
                        error stop
                    end if
#endif
                    Membership(ipstart:ipend) = Membership(ipstart:ipend) + icstart - 1_IK

                else blockSubclusterCheck

                    if (Kmeans%Size(ic) > 0_IK) then
                        KmeansNeopt(ic) = 1_IK
                        Membership(ipstart:ipend) = icstart
#if defined MAXDEN
                        if (probPoisDiff < unifrnd) then
                            LogVolNormed(icstart) = KmeansLogVolEstimate(ic)
                            scaleFactor = exp( -ndInverse * logVolRatio )
                            Kmeans%Prop%ScaleFactorSq(ic) = scaleFactor**2
                            scaleFactorSqInverse = 1._RK / Kmeans%Prop%ScaleFactorSq(ic)
                            do j = 1, nd
                                Kmeans%Prop%ChoDia(j,ic) = Kmeans%Prop%ChoDia(j,ic) * scaleFactor
                                Kmeans%Prop%InvCovMat(1:nd,j,ic) = Kmeans%Prop%InvCovMat(1:nd,j,ic) * scaleFactorSqInverse
                                do i = j + 1, nd
                                    Kmeans%Prop%ChoLowCovUpp(i,j,ic) = Kmeans%Prop%ChoLowCovUpp(i,j,ic) * scaleFactor
                                end do
                            end do
                        end if
#endif
                    else
                        KmeansNeopt(ic) = 0_IK
                    end if

                end if blockSubclusterCheck

#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
                if (any(Membership(ipstart:ipend) <= 0_IK) .or. any(Membership(ipstart:ipend) > nemax)) then
                    write(*,*) PROCEDURE_NAME
                    write(*,*) "Membership(ipstart:ipend) = 0 : ", Membership(ipstart:ipend)
                    write(*,*) "ic, icstart: ", ic, icstart
                    error stop
                end if
#endif

                blockSingleCluster: if (KmeansNeopt(ic) == 1_IK) then ! .or. KmeansNemax(ic) == 1_IK ) then ! There is only one cluster here

                    Size(icstart) = Kmeans%Size(ic)
                    Center(1:nd,icstart) = Kmeans%Center(1:nd,ic)
                    ChoDia(1:nd,icstart) = Kmeans%Prop%ChoDia(1:nd,ic)
                    ChoLowCovUpp(1:nd,1:nd,icstart) = Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic)
                    InvCovMat(1:nd,1:nd,icstart) = Kmeans%Prop%InvCovMat(1:nd,1:nd,ic)

                    ! Rescale the Cholesky factor and its diagonals such that they describe the bounded shape of the cluster

                    !scaleFactor = sqrt(Kmeans%Prop%ScaleFactorSq(ic))
                    !scaleFactorSqInverse = 1._RK / Kmeans%Prop%ScaleFactorSq(ic)
                    !do j = 1, nd
                    !    Center(j,icstart) = Kmeans%Center(j,ic)
                    !    ChoDia(j,icstart) = Kmeans%Prop%ChoDia(j,ic) * scaleFactor
                    !    InvCovMat(1:nd,j,icstart) = Kmeans%Prop%InvCovMat(1:nd,j,ic) * scaleFactorSqInverse
                    !    do i = j + 1, nd
                    !        ChoLowCovUpp(i,j,icstart) = ChoLowCovUpp(i,j,icstart) * scaleFactor
                    !    end do
                    !end do

#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
                    block
                        integer :: ip
                        logical :: isInside
                        real(RK) :: mahalSq
                        real(RK), allocatable :: NormedPoint(:)
                        do ip = ipstart, ipend
                            NormedPoint = Point(:,ip) - Center(:,icstart)
                            mahalSq = dot_product(NormedPoint,matmul(InvCovMat(:,:,icstart),NormedPoint))
                            isInside = mahalSq - 1._RK <= 1.e-6_RK
                            if (.not. isInside) then
                                ! LCOV_EXCL_START
                                write(*,"(*(g0.15,:,' '))") new_line("a"), PROCEDURE_NAME//": FATAL - POINT NOT INSIDE!, MAHAL = ", sqrt(mahalSq), new_line("a")
                                write(*,"(60(g0,:,','))") "numRecursiveCall", numRecursiveCall
                                write(*,"(60(g0,:,','))") "KmeansNemax(icstart), icstart", KmeansNemax(icstart), icstart
                                write(*,"(60(g0,:,','))") "ip, IpStart, IpEnd, size(Membership)", ip, ipstart, ipend, size(Membership(ipstart:ipend))
                                write(*,"(60(g0,:,','))") "Membership", Membership(ipstart:ipend)
                                write(*,"(60(g0,:,','))") "ChoLowCovUpp", ChoLowCovUpp(:,:,icstart)
                                write(*,"(60(g0,:,','))") "InvCovMat", InvCovMat(:,:,icstart)
                                write(*,"(60(g0,:,','))") "ChoDia", ChoDia(:,icstart)
                                write(*,"(60(g0,:,','))") "NormedPoint", NormedPoint
                                error stop
                                ! LCOV_EXCL_STOP
                            end if
                        end do
                    end block
#endif
                end if blockSingleCluster

            end do loopSubclusterSearch

            neopt = sum(KmeansNeopt)

#if defined MINVOL
        else blockSubclusterSearch  ! one cluster is better

            Membership = 1_IK
            Size(1) = np
            neopt = 1_IK
            return

        end if blockSubclusterSearch
#endif

#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
        block
            integer :: ic, ip
            logical :: isInside
            real(RK) :: mahalSq
            real(RK), allocatable :: NormedPoint(:)
            do ip = 1, np
                ic = Membership(ip)
                if (ic < 1_IK .or. ic > nemax) then
                    write(*,*) PROCEDURE_NAME
                    write(*,*) "ic, nc, ip, np, nemax, numRecursiveCall", ic, nc, ip, np, nemax, numRecursiveCall
                    write(*,*) "Kmeans%Prop%CumSumSize", Kmeans%Prop%CumSumSize
                    error stop
                end if
                NormedPoint = Point(:,ip) - Center(:,ic)
                mahalSq = dot_product(NormedPoint,matmul(InvCovMat(:,:,ic),NormedPoint))
                isInside = mahalSq - 1._RK <= 1.e-6_RK
                if (.not. isInside) then
                    ! LCOV_EXCL_START
                    write(*,"(*(g0.15,:,' '))") new_line("a"), PROCEDURE_NAME//": FATAL - POINT NOT INSIDE!, MAHALSQ = ", mahalSq, new_line("a")
                    write(*,"(60(g0,:,','))") "numRecursiveCall, ic, size(Membership)", numRecursiveCall, ic, size(Membership)
                    write(*,"(60(g0,:,','))") "Membership", Membership
                    write(*,"(60(g0,:,','))") "ChoLowCovUpp", ChoLowCovUpp(:,:,:) ! ic)
                    write(*,"(60(g0,:,','))") "InvCovMat", InvCovMat(:,:,ic)
                    write(*,"(60(g0,:,','))") "NormedPoint", NormedPoint
                    error stop
                    ! LCOV_EXCL_STOP
                end if
            end do
        end block
#endif

    end subroutine runRecursivePartition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined MAXDEN
    pure function getProbPoisDiff(count, logVolRatio) result(probPoisDiff)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbPoisDiff
#endif
        use Constants_mod, only: IK, RK, NEGINF_RK
        implicit none
        integer(IK) , intent(in)    :: count
        real(RK)    , intent(in)    :: logVolRatio
        real(RK)                    :: probPoisDiff
        probPoisDiff = count * (1._RK + logVolRatio - exp(logVolRatio))
        if (probPoisDiff < NEGINF_RK) then ! @todo: xxx NEGINF_RK could be replaced with log(epsilon) for efficiency reasons.
            probPoisDiff = 0._RK
        else
            probPoisDiff = exp(probPoisDiff)
        end if
    end function getProbPoisDiff
#endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine writePartitionMethod(self, fileUnit, nd, np, Point)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: writePartitionMethod
#endif
        use Constants_mod, only: IK, RK
        implicit none
        class(Partition_type)   , intent(in)    :: self
        integer(IK)             , intent(in)    :: fileUnit, nd, np
        real(RK)                , intent(in)    :: Point(nd,np)
        call writePartition ( fileUnit = fileUnit & ! LCOV_EXCL_LINE
                            , nd = nd & ! LCOV_EXCL_LINE
                            , np = np & ! LCOV_EXCL_LINE
                            , nc = self%neopt & ! LCOV_EXCL_LINE
                            , Size = self%Size & ! LCOV_EXCL_LINE
                            , Center = self%Center & ! LCOV_EXCL_LINE
                            , Membership = self%Membership & ! LCOV_EXCL_LINE
                            , ChoLowCovUpp = self%ChoLowCovUpp & ! LCOV_EXCL_LINE
                            , LogVolNormed = self%LogVolNormed & ! LCOV_EXCL_LINE
                            , ChoDia = self%ChoDia & ! LCOV_EXCL_LINE
                            , Point = Point & ! LCOV_EXCL_LINE
                            )
    end subroutine writePartitionMethod

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine writePartition   ( fileUnit & ! LCOV_EXCL_LINE
                                , nd,np,nc & ! LCOV_EXCL_LINE
                                , Size & ! LCOV_EXCL_LINE
                                , Center & ! LCOV_EXCL_LINE
                                , ChoDia & ! LCOV_EXCL_LINE
                                , Membership & ! LCOV_EXCL_LINE
                                , ChoLowCovUpp & ! LCOV_EXCL_LINE
                                , LogVolNormed & ! LCOV_EXCL_LINE
                                , Point & ! LCOV_EXCL_LINE
                                )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: writePartition
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK)     , intent(in)    :: fileUnit
        integer(IK)     , intent(in)    :: nd, np, nc
        integer(IK)     , intent(in)    :: Size(nc)
        real(RK)        , intent(in)    :: Center(nd,nc)
        real(RK)        , intent(in)    :: ChoDia(nd,nc)
        integer(IK)     , intent(in)    :: Membership(np)
        real(RK)        , intent(in)    :: LogVolNormed(nc)
        real(RK)        , intent(in)    :: ChoLowCovUpp(nd,nd,nc)
        real(RK)        , intent(in)    :: Point(nd,np)
        character(*)    , parameter     :: fileFormat = "(*(g0.8,:,','))"
        integer(IK)                     :: i, j, ic

        write(fileUnit,"(*(g0,:,'"//new_line("a")//"'))") "nd", nd, "np", np, "nc", nc

        write(fileUnit,"(A)") "Size"
        write(fileUnit,fileFormat) Size

        write(fileUnit,"(A)") "Center"
        write(fileUnit,fileFormat) Center

        write(fileUnit,"(A)") "LogVolume"
        write(fileUnit,fileFormat) LogVolNormed

        write(fileUnit,"(A)") "CholeskyLower"
        write(fileUnit,fileFormat) ((ChoDia(j,ic), (ChoLowCovUpp(i,j,ic), i=j+1,nd), j=1,nd), ic=1,nc)

        write(fileUnit,"(A)") "Point"
        write(fileUnit,fileFormat) Point

        write(fileUnit,"(A)") "Membership"
        write(fileUnit,fileFormat) Membership

    end subroutine writePartition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
