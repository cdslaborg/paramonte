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

    use Constants_mod, only: RK
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
        integer                     :: nsim                         !< See the interface of [constructPartition()](@ref constructpartition).
        integer                     :: nd,np                        !< See the interface of [constructPartition()](@ref constructpartition).
        integer                     :: nc,nt                        !< See the interface of [constructPartition()](@ref constructpartition).
        integer                     :: nemax                        !< See the interface of [constructPartition()](@ref constructpartition).
        integer                     :: minSize                      !< See the interface of [constructPartition()](@ref constructpartition).
#if defined MAXDEN
        integer                     :: optimizationLevel
#endif
        integer                     :: maxAllowedKmeansFailure      !< See the interface of [constructPartition()](@ref constructpartition).
        integer                     :: maxAllowedKmeansRecursion    !< See the interface of [constructPartition()](@ref constructpartition).
        integer                     :: maxAllowedKvolumeRecursion   !< See the interface of [constructPartition()](@ref constructpartition).
        real(RK)                    :: mahalSqWeightExponent        !< See the interface of [constructPartition()](@ref constructpartition).
        real(RK)                    :: parentLogVolNormed           !< See the interface of [constructPartition()](@ref constructpartition).
        real(RK)                    :: pointLogVolNormed            !< The estimated volume of a single point based on the input value of `parentLogVolNormed`.
        real(RK)                    :: inclusionFraction            !< See the interface of [constructPartition()](@ref constructpartition).
        real(RK)                    :: logExpansion                 !< See the interface of [constructPartition()](@ref constructpartition).
        real(RK)                    :: logShrinkage                 !< See the interface of [constructPartition()](@ref constructpartition).
        real(RK)                    :: kmeansRelTol                 !< See the interface of [constructPartition()](@ref constructpartition).
        logical                     :: stanEnabled                  !< See the interface of [constructPartition()](@ref constructpartition).
        integer                     :: neopt                        !< The predicted optimal number of clusters identified in the input data.
        integer                     :: numRecursiveCall             !< The total number of recursive calls to the partitioning algorithm.
        integer                     :: convergenceFailureCount      !< The number of times the partitioning algorithm has failed to converge.
        integer     , allocatable   :: Membership(:)                !< An array of size `(np)` representing the bounding-ellipsoid membership IDs of the corresponding data points.
        integer     , allocatable   :: PointIndex(:)                !< An array of size `(np)` of indices such that the input `Point(1:nd,Index(ip))` is saved at the output `Point(1:nd,ip)`.
#if defined MAXDEN
        integer     , allocatable   :: EffectiveSize(:)             !< An array of size `(nemax)` indicating the likelihood of subclustering (`0<<1` if successful, else `1<<2` if failed), enlargement (if < 0), being warranted. If it is warranted but fails, the likelihood will be negative.
        real(RK)    , allocatable   :: LogLikeFitness(:)            !< An array of size `(nemax)` indicating the likelihood of subclustering (`0<<1` if successful, else `1<<2` if failed), enlargement (if < 0), being warranted. If it is warranted but fails, the likelihood will be negative.
#endif
        real(RK)    , allocatable   :: LogVolNormed(:)              !< An array of size `(nemax)` representing the log-volumes of the corresponding bounding ellipsoids.
        real(RK)    , allocatable   :: ChoLowCovUpp(:,:,:)          !< An array of size `(nd,nd,nemax)` representing the Cholesky lower triangle, diagonal, and covariance matrices of the bounding ellipsoids.
        real(RK)    , allocatable   :: InvCovMat(:,:,:)             !< An array of size `(nd,nd,nemax)` representing the full symmetric inverse covariance matrices of the bounding ellipsoids.
        real(RK)    , allocatable   :: ChoDia(:,:)                  !< An array of size `(nd,nemax)` representing the Cholesky lower triangle, diagonal, and covariance matrices of the bounding ellipsoids.
        real(RK)    , allocatable   :: Center(:,:)                  !< An array of size `(nd,nemax)` representing the centers of the bounding ellipsoids.
        integer     , allocatable   :: Size(:)                      !< An array of size `(nemax)` representing the sizes of the corresponding bounding ellipsoids.
        type(Err_type)              :: Err                          !< An object of class [Err_type](@ref err_mod::err_type) containing information about error occurrence.
    contains
        procedure, pass             :: write => writePartition
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
    !> \param[in]   logShrinkage                :   The logarithm of the factor by which the sum of the volumes of the children bounds must be smaller than their parents to warrant further sub-partitioning (**optional**, default = `0.`).
    !> \param[in]   inclusionFraction           :   The fraction of non-member points that **are** inside the partition, to be used in estimating of the true volumes of the partitions (**optional**, default = `1.`).
    !> \param[in]   parentLogVolNormed          :   The logarithm of the volume of the input set of points normalized to the volume of the unit `nd`-ball (**optional**, default = `-infinity`).
    !> \param[in]   mahalSqWeightExponent       :   The exponent with which the density of the partitions are exponentiated to be used as the Mahalanobis distance weights (**optional**, default = `0.`).
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
                                , logShrinkage & ! LCOV_EXCL_LINE
                                , inclusionFraction & ! LCOV_EXCL_LINE
#if defined MAXDEN
                                , optimizationLevel & ! LCOV_EXCL_LINE
#endif
                                , parentLogVolNormed & ! LCOV_EXCL_LINE
                                , mahalSqWeightExponent & ! LCOV_EXCL_LINE
                                , maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                , maxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                                , maxAllowedKvolumeRecursion & ! LCOV_EXCL_LINE
                                ) result(Partition)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructPartition
#endif
        use Constants_mod, only: IK, RK, NEGINF_RK

        implicit none

        real(RK)    , intent(inout)             :: Point(:,:)
        integer     , intent(in)    , optional  :: nc,nt
        integer     , intent(in)    , optional  :: nsim
        integer     , intent(in)    , optional  :: nemax
        integer     , intent(in)    , optional  :: minSize
        logical     , intent(in)    , optional  :: stanEnabled
        logical     , intent(in)    , optional  :: trimEnabled
        real(RK)    , intent(in)    , optional  :: kmeansRelTol
        real(RK)    , intent(in)    , optional  :: logExpansion
        real(RK)    , intent(in)    , optional  :: logShrinkage
        real(RK)    , intent(in)    , optional  :: inclusionFraction
        real(RK)    , intent(in)    , optional  :: parentLogVolNormed
        real(RK)    , intent(in)    , optional  :: mahalSqWeightExponent
        integer     , intent(in)    , optional  :: maxAllowedKmeansFailure
        integer     , intent(in)    , optional  :: maxAllowedKmeansRecursion
        integer     , intent(in)    , optional  :: maxAllowedKvolumeRecursion
#if defined MAXDEN
        integer     , intent(in)    , optional  :: optimizationLevel
#endif
        type(Partition_type)                    :: Partition

        Partition%nd = size(Point(:,1))
        Partition%np = size(Point(1,:))

        Partition%nc = 2; if (present(nc)) Partition%nc = max(Partition%nc, nc)
        Partition%nt = 1; if (present(nt)) Partition%nt = max(Partition%nt, nt)
        Partition%nsim = 0; if (present(nsim)) Partition%nsim = max(Partition%nsim, nsim)
        Partition%minSize = Partition%nd + 1; if (present(minSize)) Partition%minSize = minSize
        Partition%nemax = Partition%np / max(1, Partition%minSize); if (present(nemax)) Partition%nemax = max(1, nemax)
#if defined MAXDEN
        Partition%optimizationLevel = 1; if (present(optimizationLevel)) Partition%optimizationLevel = optimizationLevel
#endif
        Partition%maxAllowedKvolumeRecursion = 3; if (present(maxAllowedKvolumeRecursion)) Partition%maxAllowedKvolumeRecursion = maxAllowedKvolumeRecursion
        Partition%maxAllowedKmeansRecursion = 300; if (present(maxAllowedKmeansRecursion)) Partition%maxAllowedKmeansRecursion = maxAllowedKmeansRecursion
        Partition%maxAllowedKmeansFailure = 10; if (present(maxAllowedKmeansFailure)) Partition%maxAllowedKmeansFailure = maxAllowedKmeansFailure
        Partition%mahalSqWeightExponent = 0._RK; if (present(mahalSqWeightExponent)) then; if (abs(mahalSqWeightExponent)>1.e-4_RK) Partition%mahalSqWeightExponent = mahalSqWeightExponent; endif
        Partition%inclusionFraction = 0._RK; if (present(inclusionFraction)) then; if (abs(inclusionFraction)>1.e-5_RK) Partition%inclusionFraction = inclusionFraction; endif
        Partition%logExpansion = 0._RK; if (present(logExpansion)) Partition%logExpansion = logExpansion
        Partition%logShrinkage = 0._RK; if (present(logShrinkage)) Partition%logShrinkage = logShrinkage
        Partition%kmeansRelTol = 1.e-4_RK; if (present(kmeansRelTol)) Partition%kmeansRelTol = kmeansRelTol
        Partition%stanEnabled = .true.; if (present(stanEnabled)) Partition%stanEnabled = stanEnabled

        if (Partition%mahalSqWeightExponent == 0._RK .and. .not. present(parentLogVolNormed)) then ! xxx could we instead approximate `parentLogVolNormed`?
            Partition%Err%msg = "If `mahalSqWeightExponent` is given as input argument, then `parentLogVolNormed` must be also given as input."
            Partition%Err%occurred = .true.
            error stop
        end if

        allocate( Partition%Size            (Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%Center          (Partition%nd,Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%ChoDia          (Partition%nd,Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%InvCovMat       (Partition%nd,Partition%nd,Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%ChoLowCovUpp    (Partition%nd,Partition%nd,Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%LogVolNormed    (Partition%nemax) & ! LCOV_EXCL_LINE
#if defined MAXDEN
                , Partition%EffectiveSize   (Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%LogLikeFitness  (Partition%nemax) & ! LCOV_EXCL_LINE
#endif
                , Partition%Membership      (Partition%np) & ! LCOV_EXCL_LINE
                , Partition%PointIndex      (Partition%np) & ! LCOV_EXCL_LINE
                )

        call Partition%run(Point, parentLogVolNormed)

        if (present(trimEnabled)) then
            if (trimEnabled) then
                Partition%Size              = Partition%Size            (1:Partition%neopt)
                Partition%Center            = Partition%Center          (1:Partition%nd,1:Partition%neopt)
                Partition%ChoDia            = Partition%ChoDia          (1:Partition%nd,1:Partition%neopt)
                Partition%InvCovMat         = Partition%InvCovMat       (1:Partition%nd,1:Partition%nd,1:Partition%neopt)
                Partition%ChoLowCovUpp      = Partition%ChoLowCovUpp    (1:Partition%nd,1:Partition%nd,1:Partition%neopt)
                Partition%LogVolNormed      = Partition%LogVolNormed    (1:Partition%neopt)
#if defined MAXDEN
                Partition%EffectiveSize     = Partition%EffectiveSize   (1:Partition%neopt)
                Partition%LogLikeFitness    = Partition%LogLikeFitness  (1:Partition%neopt)
#endif
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
        use Constants_mod, only: IK, RK, NEGINF_RK, POSINF_RK
        use Matrix_mod, only: getInvMatFromCholFac
        use Matrix_mod, only: getCholeskyFactor
        use String_mod, only: num2str
        use Sort_mod, only: indexArray!, sortAscending

        implicit none

        class(Partition_type)   , intent(inout) :: Partition
        real(RK)                , intent(inout) :: Point(Partition%nd,Partition%np) ! This must appear after Partition declaration.
#if defined MINVOL
        real(RK)    , optional  , intent(in)    :: parentLogVolNormed
#elif defined MAXDEN
        real(RK)                , intent(in)    :: parentLogVolNormed
        ! optimization
        real(RK)                                :: ndHalf
        real(RK)                                :: unifrnd
        real(RK)                                :: logVolRatio
        real(RK)                                :: logVolNormed
        real(RK)                                :: ndHalfInverse
        real(RK)                                :: logLikeFitness
        real(RK)                                :: Center(Partition%nd)
        real(RK)                                :: ChoDia(Partition%nd)
        real(RK)                                :: InvCovMat(Partition%nd,Partition%nd)
        real(RK)                                :: ChoLowCovUpp(Partition%nd,Partition%nd)
        integer                                 :: MahalSqSortedIndex(Partition%np)
        integer                                 :: effectiveSize

        real(RK)                                :: MahalSq(Partition%np,Partition%nemax)
       !real(RK)                                :: logLikeShrinkage
       !real(RK)                                :: SumPoint(Partition%nd)
       !real(RK)                                :: SinglePartitionCenter(Partition%nd)
       !real(RK)                                :: MahalSqSorted(Partition%np)
       !real(RK)    , allocatable               :: MahalSqAcrossPartition(:)
       !integer                                 :: newEffeciveSize
       !integer                                 :: ipstart, ipend
        integer                                 :: nps, npe, npc
       !integer                                 :: ipMaxMahal
       !integer                                 :: icMinMahal
       !integer                                 :: newSize
       !integer                                 :: jp
#endif
        integer                                 :: i, j, ip, ic
        integer , allocatable                   :: PointIndexSorted(:)
        real(RK)                                :: NormedPoint(Partition%nd,Partition%np)
        logical                                 :: boundedRegionIsTooLarge
        real(RK)                                :: PointAvg(Partition%nd)
        real(RK)                                :: PointStd(Partition%nd)
        real(RK)                                :: scaleFactorSqInverse
        real(RK)                                :: scaleFactorSq
        real(RK)                                :: scaleFactor
        real(RK)                                :: MahalSqSinglePartition(Partition%np)
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
            MahalSqSinglePartition(ip) = dot_product( NormedPoint(1:Partition%nd,ip) , matmul(Partition%InvCovMat(1:Partition%nd,1:Partition%nd,1), NormedPoint(1:Partition%nd,ip)) )
            if (scaleFactorSq < MahalSqSinglePartition(ip)) scaleFactorSq = MahalSqSinglePartition(ip)
        end do
        scaleFactor = sqrt(scaleFactorSq)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the largest possible LogVolume of the first cluster based on the parent cluster and the input estimated volume.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Partition%LogVolNormed(1) = sum( log( Partition%ChoDia(1:Partition%nd,1) ) ) + Partition%nd * log(scaleFactor)

        ! Rescale the scaleFactor, if needed, to enlarge the bounded cluster in the future.
        ! This must be done here, although it will be used at the very end.
        ! Otherwise, Partition%LogVolNormed(1) will be potentially overwritten.

#if defined MAXDEN
        Partition%parentLogVolNormed = Partition%logExpansion + parentLogVolNormed
        logVolRatio = Partition%LogVolNormed(1) - Partition%parentLogVolNormed
        Partition%pointLogVolNormed = Partition%parentLogVolNormed - log(real(Partition%np,RK))
        Partition%LogLikeFitness(1) = getLogProbPoisDiff(Partition%np, logVolRatio)
#if defined DEBUG_ENABLED || defined TESTING_ENABLED || defined CODECOVE_ENABLED
        if (Partition%LogLikeFitness(1) < NEGINF_RK .or. Partition%LogLikeFitness(1) > 0._RK) then
            write(*,*) "Internal error occurred: Partition%LogLikeFitness(1) < NEGINF_RK .or. Partition%LogLikeFitness(1) > 0._RK"
            write(*,*) "Partition%LogVolNormed(1), Partition%parentLogVolNormed, Partition%LogLikeFitness(1):" & ! LCOV_EXCL_LINE
                     ,  Partition%LogVolNormed(1), Partition%parentLogVolNormed, Partition%LogLikeFitness(1)
            error stop
        end if
#endif
        boundedRegionIsTooLarge = logVolRatio > 0._RK
        if (boundedRegionIsTooLarge) then
            call random_number(unifrnd)
            if (unifrnd > 0._RK) then; unifrnd = log(unifrnd); else; unifrnd = -huge(unifrnd); end if ! LCOV_EXCL_LINE or unifrnd = log(max(tiny(unifrnd),unifrnd))
            boundedRegionIsTooLarge = Partition%LogLikeFitness(1) < unifrnd
        else ! enlarge, no question asked.
            !Partition%LogLikeFitness(1) = 1._RK + Partition%LogLikeFitness(1)
            Partition%LogLikeFitness(1) = 0._RK
            scaleFactor = scaleFactor * exp( (Partition%parentLogVolNormed - Partition%LogVolNormed(1)) / Partition%nd )
            Partition%LogVolNormed(1) = Partition%parentLogVolNormed
            scaleFactorSq = scaleFactor**2
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

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Search for subclusters.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        blockSubclusterSearch: if (Partition%nemax > 1 .and. boundedRegionIsTooLarge) then

            Partition%numRecursiveCall = 0
            Partition%convergenceFailureCount = 0

            blockStandardization: if (Partition%stanEnabled) then

                PointAvg(1:Partition%nd) = Partition%Center(1:Partition%nd,1)
                call runRecursivePartition  ( nd    = Partition%nd & ! LCOV_EXCL_LINE
                                            , np    = Partition%np & ! LCOV_EXCL_LINE
                                            , nc    = min(Partition%nc, Partition%nemax) & ! LCOV_EXCL_LINE
                                            , nt    = Partition%nt & ! LCOV_EXCL_LINE
                                            , npc   = Partition%np & ! LCOV_EXCL_LINE
                                            , nps   = 1 & ! LCOV_EXCL_LINE
                                            , npe   = Partition%np & ! LCOV_EXCL_LINE
                                            , nsim  = Partition%nsim & ! LCOV_EXCL_LINE
                                            , nemax = Partition%nemax & ! LCOV_EXCL_LINE
                                            , minSize = Partition%minSize & ! LCOV_EXCL_LINE
                                            , kmeansRelTol = Partition%kmeansRelTol & ! LCOV_EXCL_LINE
                                            , logExpansion = Partition%logExpansion & ! LCOV_EXCL_LINE
                                            , logShrinkage = Partition%logShrinkage & ! LCOV_EXCL_LINE
                                            , inclusionFraction = Partition%inclusionFraction & ! LCOV_EXCL_LINE
                                            , pointLogVolNormed = Partition%pointLogVolNormed & ! LCOV_EXCL_LINE
                                            , parentLogVolNormed = Partition%parentLogVolNormed & ! LCOV_EXCL_LINE
                                            , mahalSqWeightExponent = Partition%mahalSqWeightExponent & ! LCOV_EXCL_LINE
                                            , maxAllowedKmeansFailure = Partition%maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                            , maxAllowedKmeansRecursion = Partition%maxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                                            , maxAllowedKvolumeRecursion = Partition%maxAllowedKvolumeRecursion & ! LCOV_EXCL_LINE
                                            , Point = Point & ! LCOV_EXCL_LINE
                                            , PartitionSize = Partition%Size & ! LCOV_EXCL_LINE
                                            , neopt = Partition%neopt & ! LCOV_EXCL_LINE
                                            , Center = Partition%Center & ! LCOV_EXCL_LINE
                                            , ChoDia = Partition%ChoDia & ! LCOV_EXCL_LINE
#if defined MAXDEN
                                            , LogLikeFitness = Partition%LogLikeFitness & ! LCOV_EXCL_LINE
                                            , EffectiveSize = Partition%EffectiveSize & ! LCOV_EXCL_LINE
                                            , MahalSq = MahalSq & ! LCOV_EXCL_LINE
#endif
                                            , InvCovMat = Partition%InvCovMat & ! LCOV_EXCL_LINE
                                            , Membership = Partition%Membership & ! LCOV_EXCL_LINE
                                            , PointIndex = Partition%PointIndex & ! LCOV_EXCL_LINE
                                            , ChoLowCovUpp = Partition%ChoLowCovUpp & ! LCOV_EXCL_LINE
                                            , LogVolNormed = Partition%LogVolNormed & ! LCOV_EXCL_LINE
                                            , numRecursiveCall = Partition%numRecursiveCall & ! LCOV_EXCL_LINE
                                            , convergenceFailureCount = Partition%convergenceFailureCount & ! LCOV_EXCL_LINE
                                            , PointAvg = PointAvg(1:Partition%nd) & ! LCOV_EXCL_LINE
                                            , PointStd = PointStd & ! LCOV_EXCL_LINE
                                            )
                !block
                !use Timer_mod, only: Timer_type
                !type(Timer_type) :: Timer
                !call Timer%tic()
                ! The following step takes about 1 microseconds (in low dims) to 1 ms (in high dims).
                ! Sorting PointIndex at dimensions > 20 with 500-2000 tends to be superior to sorting Point.
                !xx!if (Partition%nd < 20) then
                !xx!    NormedPoint(1:Partition%nd, 1:Partition%np) = Point(1:Partition%nd, Partition%PointIndex(1:Partition%np))
                !xx!    Point(1:Partition%nd, 1:Partition%np) = NormedPoint(1:Partition%nd, 1:Partition%np)
                !xx!else
                !xx!    if (.not. allocated(PointIndexSorted)) allocate(PointIndexSorted(Partition%np))
                !xx!    call indexArray(Partition%np, Partition%PointIndex, PointIndexSorted, Partition%Err)
                !xx!    if (Partition%Err%occurred) return ! LCOV_EXCL_LINE
                !xx!    Partition%Membership(1:Partition%np) = Partition%Membership(PointIndexSorted)
                !xx!end if
                !call Timer%toc()
                !write(*,*) "delta = ", Timer%Time%delta
                !end block

            else blockStandardization

                call runRecursivePartition  ( nd    = Partition%nd & ! LCOV_EXCL_LINE
                                            , np    = Partition%np & ! LCOV_EXCL_LINE
                                            , nc    = min(Partition%nc, Partition%nemax) & ! LCOV_EXCL_LINE
                                            , nt    = Partition%nt & ! LCOV_EXCL_LINE
                                            , npc   = Partition%np & ! LCOV_EXCL_LINE
                                            , nps   = 1 & ! LCOV_EXCL_LINE
                                            , npe   = Partition%np & ! LCOV_EXCL_LINE
                                            , nsim  = Partition%nsim & ! LCOV_EXCL_LINE
                                            , nemax = Partition%nemax & ! LCOV_EXCL_LINE
                                            , minSize = Partition%minSize & ! LCOV_EXCL_LINE
                                            , kmeansRelTol = Partition%kmeansRelTol & ! LCOV_EXCL_LINE
                                            , logExpansion = Partition%logExpansion & ! LCOV_EXCL_LINE
                                            , logShrinkage = Partition%logShrinkage & ! LCOV_EXCL_LINE
                                            , inclusionFraction = Partition%inclusionFraction & ! LCOV_EXCL_LINE
                                            , pointLogVolNormed = Partition%pointLogVolNormed & ! LCOV_EXCL_LINE
                                            , parentLogVolNormed = Partition%parentLogVolNormed & ! LCOV_EXCL_LINE
                                            , mahalSqWeightExponent = Partition%mahalSqWeightExponent & ! LCOV_EXCL_LINE
                                            , maxAllowedKmeansFailure = Partition%maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                            , maxAllowedKmeansRecursion = Partition%maxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                                            , maxAllowedKvolumeRecursion = Partition%maxAllowedKvolumeRecursion & ! LCOV_EXCL_LINE
                                            , Point = Point & ! LCOV_EXCL_LINE
                                            , PartitionSize = Partition%Size & ! LCOV_EXCL_LINE
                                            , neopt = Partition%neopt & ! LCOV_EXCL_LINE
                                            , Center = Partition%Center & ! LCOV_EXCL_LINE
                                            , ChoDia = Partition%ChoDia & ! LCOV_EXCL_LINE
#if defined MAXDEN
                                            , LogLikeFitness = Partition%LogLikeFitness & ! LCOV_EXCL_LINE
                                            , EffectiveSize = Partition%EffectiveSize & ! LCOV_EXCL_LINE
                                            , MahalSq = MahalSq & ! LCOV_EXCL_LINE
#endif
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

            Partition%neopt = 1

        end if blockSubclusterSearch

#if defined DEBUG_ENABLED || defined TESTING_ENABLED || defined CODECOVE_ENABLED
        if (Partition%neopt < 1) then
            write(*,*) PROCEDURE_NAME//": Internal error occurred : Partition%neopt < 1 : ", Partition%neopt
            error stop
        end if
#endif

        blockPartitionPostProcessing: if (Partition%neopt == 1) then

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Rescale the Cholesky lower and diagonals to describe the bounding ellipsoid of the single partition.
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            Partition%Size(1) = Partition%np
            Partition%Membership(1:Partition%np) = 1
            scaleFactorSqInverse = 1._RK / scaleFactorSq
            do j = 1, Partition%nd
                Partition%ChoDia(j,1) = Partition%ChoDia(j,1) * scaleFactor
                !Partition%ChoLowCovUpp(1:j,j,1) = Partition%ChoLowCovUpp(1:j,j,1) * scaleFactorSq
                Partition%InvCovMat(1:Partition%nd,j,1) = Partition%InvCovMat(1:Partition%nd,j,1) * scaleFactorSqInverse
                do i = j+1, Partition%nd
                    Partition%ChoLowCovUpp(i,j,1) = Partition%ChoLowCovUpp(i,j,1) * scaleFactor
                end do
            end do

#if defined MAXDEN
            MahalSq(1:Partition%np,1) = MahalSqSinglePartition(1:Partition%np) * scaleFactorSqInverse
            Partition%EffectiveSize(1) = Partition%np

        else blockPartitionPostProcessing

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! maximize and refine partitions.
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            !if (.not. allocated(MahalSqSortedIndex)) allocate(MahalSqSortedIndex(Partition%np))
            !if (.not. allocated(MahalSqAcrossPartition)) allocate(MahalSqAcrossPartition(Partition%neopt))
            !SinglePartitionCenter(1:Partition%nd) = Partition%Center(1:Partition%nd,1)

            blockOptimizationEnabled: if (Partition%optimizationLevel > 0) then

                ndHalf = 0.5_RK * Partition%nd
                ndHalfInverse = 1._RK / ndHalf

                loopOverPartition: do ic = 1, Partition%neopt

                    blockNormalPartition: if (Partition%LogLikeFitness(ic) <= 0._RK) then

                        loopCurrentPartitionOptimization: do

                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            ! Do scale optimization.
                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                            !NewCenter(1:Partition%nd) = SinglePartitionCenter(1:Partition%nd)

                            call indexArray(Partition%np, MahalSq(1:Partition%np,ic), MahalSqSortedIndex(1:Partition%np), Partition%Err)
                            !MahalSqSorted(1:Partition%np) = MahalSq(1:Partition%np,ic)
                            !call sortAscending(Partition%np, MahalSq(1:Partition%np,ic), Partition%Err)
                            if (Partition%Err%occurred) cycle loopOverPartition

                            npe = Partition%np
                            nps = Partition%EffectiveSize(ic)
                            if (npe - nps < 2) cycle loopOverPartition ! xxx is this necessary?

                            loopScaleOptimization: do

                                npc = (nps + npe) / 2

                                logVolNormed = Partition%LogVolNormed(ic) + ndHalf * log(MahalSq(MahalSqSortedIndex(npc),ic))
                                logVolRatio = logVolNormed - Partition%pointLogVolNormed - log(real(npc,RK))
                                boundedRegionIsTooLarge = logVolRatio > 0._RK
                                if (boundedRegionIsTooLarge) then
                                    logLikeFitness = getLogProbPoisDiff(npc, logVolRatio)
                                    !call random_number(unifrnd)
                                    !if (unifrnd > 0._RK) then; unifrnd = log(unifrnd); else; unifrnd = -huge(unifrnd); end if ! LCOV_EXCL_LINE
                                    !boundedRegionIsTooLarge = logLikeFitness < unifrnd .and. Partition%LogLikeFitness(ic) > logLikeFitness
                                    boundedRegionIsTooLarge = Partition%LogLikeFitness(ic) > logLikeFitness
                                    if (boundedRegionIsTooLarge) then
                                        call random_number(unifrnd)
                                        if (unifrnd > 0._RK) then; unifrnd = log(unifrnd); else; unifrnd = -huge(unifrnd); end if ! LCOV_EXCL_LINE
                                        boundedRegionIsTooLarge = logLikeFitness - Partition%LogLikeFitness(ic) < unifrnd
                                        !boundedRegionIsTooLarge = logLikeFitness < unifrnd
                                    end if
                                end if

                                if (boundedRegionIsTooLarge) then
                                    npe = npc
                                else
                                    nps = npc
                                end if

                                if (npe - nps < 2) then
                                    if (boundedRegionIsTooLarge) then
                                        npc = nps
                                    else
                                        npc = npe
                                    end if
                                    write(*,*) "ic, npc, Partition%EffectiveSize(ic)", ic, npc, Partition%EffectiveSize(ic)
                                    if (Partition%EffectiveSize(ic) == npc) exit loopCurrentPartitionOptimization
#if defined DEBUG_ENABLED || defined TESTING_ENABLED || defined CODECOVE_ENABLED
                                    if (Partition%EffectiveSize(ic) > npc) then
                                        write(*,*) "Partition%EffectiveSize(ic) > npc: ", Partition%EffectiveSize(ic), npc
                                        error stop
                                    end if
#endif
                                    Partition%EffectiveSize(ic) = npc
                                    Partition%LogVolNormed(ic) = Partition%LogVolNormed(ic) + ndHalf * log(MahalSq(MahalSqSortedIndex(npc),ic))
                                    logVolRatio = Partition%LogVolNormed(ic) - Partition%pointLogVolNormed - log(real(npc,RK))
                                    Partition%LogLikeFitness(ic) = getLogProbPoisDiff(npc, logVolRatio)
                                    if (logVolRatio < 0._RK) then
                                        Partition%LogLikeFitness(ic) = 0._RK
                                        scaleFactorSq = exp(ndHalfInverse * logVolRatio)
                                        Partition%LogVolNormed(ic) = Partition%LogVolNormed(ic) - logVolRatio
                                    else
                                        scaleFactorSq = MahalSq(MahalSqSortedIndex(Partition%EffectiveSize(ic)),ic)
                                        Partition%LogVolNormed(ic) = Partition%LogVolNormed(ic) - logVolRatio
                                    end if
                                    scaleFactor = sqrt(scaleFactorSq)
                                    scaleFactorSqInverse = 1._RK / scaleFactorSq
                                    do j = 1, Partition%nd
                                        Partition%ChoDia(j,ic) = Partition%ChoDia(j,ic) * scaleFactor
                                        !Partition%ChoLowCovUpp(1:j,j,ic) = Partition%ChoLowCovUpp(1:j,j,ic) * scaleFactorSq
                                        Partition%InvCovMat(1:Partition%nd,j,ic) = Partition%InvCovMat(1:Partition%nd,j,ic) * scaleFactorSqInverse
                                        do i = j+1, Partition%nd
                                            Partition%ChoLowCovUpp(i,j,ic) = Partition%ChoLowCovUpp(i,j,ic) * scaleFactor
                                        end do
                                    end do
                                    exit loopScaleOptimization
                                end if

                            end do loopScaleOptimization

                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            ! Do shape optimization
                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                            blockShapeOptimizationEnabled: if (Partition%optimizationLevel > 1) then

                                blockTempShapeOptimization: block

                                !real(RK) :: logVolNormed
                                !integer  :: effectiveSize
                                !real(RK) :: logLikeFitness
                                !real(RK) :: Center(Partition%nd)
                                !real(RK) :: ChoDia(Partition%nd)
                                !real(RK) :: InvCovMat(Partition%nd,Partition%nd)
                                !real(RK) :: ChoLowCovUpp(Partition%nd,Partition%nd)
                                !real(RK) :: PartitionPoint(1:Partition%nd,Partition%EffectiveSize(ic))

                                NormedPoint(1:Partition%nd,1:Partition%np) = Point(1:Partition%nd,MahalSqSortedIndex)

                                ! Initialize the center of the first cluster

                                Center(1:Partition%nd) = 0._RK
                                do ip = 1, Partition%EffectiveSize(ic)
                                    Center(1:Partition%nd) = Center(1:Partition%nd) + NormedPoint(1:Partition%nd,ip)
                                end do
                                Center(1:Partition%nd) = Center(1:Partition%nd) / Partition%EffectiveSize(ic)

                                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                ! Compute the covariance matrix and the Cholesky factorization of the first cluster.
                                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                ChoLowCovUpp(1:Partition%nd,1:Partition%nd) = 0._RK
                                do ip = 1, Partition%EffectiveSize(ic)
                                    do j = 1, Partition%nd
                                        NormedPoint(j,ip) = NormedPoint(j,ip) - Center(j)
                                        ChoLowCovUpp(1:j,j) = ChoLowCovUpp(1:j,j) + NormedPoint(1:j,ip) * NormedPoint(j,ip)
                                    end do
                                end do
                                ChoLowCovUpp(1:Partition%nd,1:Partition%nd) = ChoLowCovUpp(1:Partition%nd,1:Partition%nd) / (Partition%EffectiveSize(ic) - 1)

                                ! Compute the Cholesky factor lower.

                                call getCholeskyFactor(Partition%nd, ChoLowCovUpp(1:Partition%nd,1:Partition%nd), ChoDia(1:Partition%nd))
                                if (ChoDia(1)<0._RK) then
                                    ! LCOV_EXCL_START
                                    Partition%Err%msg = PROCEDURE_NAME//": Singular covariance matrix detected for nemax, nd, Partition%EffectiveSize(ic): "//num2str(Partition%nemax)//", "//num2str(Partition%nd)//", "//num2str(Partition%EffectiveSize(ic))
                                    write(*,*) Partition%Err%msg ! xxx This needs to be fixed.
                                    error stop
                                    ! LCOV_EXCL_STOP
                                end if

                                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                ! Compute the inverse covariance matrix.
                                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                InvCovMat(1:Partition%nd,1:Partition%nd) = getInvMatFromCholFac ( Partition%nd & ! LCOV_EXCL_LINE
                                                                                                , ChoLowCovUpp(1:Partition%nd,1:Partition%nd) & ! LCOV_EXCL_LINE
                                                                                                , ChoDia(1:Partition%nd) & ! LCOV_EXCL_LINE
                                                                                                )

                                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                ! Compute the Mahal distances squared and the scale factor which makes the covariance matrix a bounding ellipsoid.
                                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                scaleFactorSq = NEGINF_RK
                                do ip = 1, Partition%EffectiveSize(ic)
                                    MahalSqSinglePartition(ip) = dot_product( NormedPoint(1:Partition%nd,ip) , matmul(InvCovMat(1:Partition%nd,1:Partition%nd), NormedPoint(1:Partition%nd,ip)) )
                                    if (scaleFactorSq < MahalSqSinglePartition(ip)) scaleFactorSq = MahalSqSinglePartition(ip)
                                end do
                                !scaleFactor = sqrt(scaleFactorSq)

                                do concurrent(ip = Partition%EffectiveSize(ic)+1:Partition%np)
                                    NormedPoint(1:Partition%nd,ip) = NormedPoint(1:Partition%nd,ip) - Center(1:Partition%nd)
                                    MahalSqSinglePartition(ip) = dot_product( NormedPoint(1:Partition%nd,ip) , matmul(InvCovMat(1:Partition%nd,1:Partition%nd), NormedPoint(1:Partition%nd,ip)) )
                                end do

                                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                ! Compute the effective fraction of points inside the bounded region.
                                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                if (Partition%inclusionFraction > 0._RK) then
                                    effectiveSize   = Partition%EffectiveSize(ic) & ! LCOV_EXCL_LINE
                                                    + nint(Partition%inclusionFraction * & ! LCOV_EXCL_LINE
                                                    ( count(MahalSqSinglePartition(Partition%EffectiveSize(ic)+1:Partition%np)<=scaleFactorSq) & ! LCOV_EXCL_LINE
                                                    ) )
                                else
                                    effectiveSize = Partition%EffectiveSize(ic)
                                end if

                                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                ! Assess the partition size.
                                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                logVolNormed = sum(log(ChoDia)) + ndHalf * log(scaleFactorSq)
                                logVolRatio = logVolNormed - Partition%pointLogVolNormed - log(real(effectiveSize,RK))
                                boundedRegionIsTooLarge = logVolRatio > 0._RK
                                if (boundedRegionIsTooLarge) then
                                    logLikeFitness = getLogProbPoisDiff(effectiveSize, logVolRatio)
                                    !call random_number(unifrnd)
                                    !if (unifrnd > 0._RK) then; unifrnd = log(unifrnd); else; unifrnd = -huge(unifrnd); end if ! LCOV_EXCL_LINE
                                    !boundedRegionIsTooLarge = logLikeFitness < unifrnd .and. Partition%LogLikeFitness(ic) > logLikeFitness
                                    boundedRegionIsTooLarge = Partition%LogLikeFitness(ic) > logLikeFitness
                                    if (boundedRegionIsTooLarge) then
                                        call random_number(unifrnd)
                                        if (unifrnd > 0._RK) then; unifrnd = log(unifrnd); else; unifrnd = -huge(unifrnd); end if ! LCOV_EXCL_LINE
                                        boundedRegionIsTooLarge = logLikeFitness - Partition%LogLikeFitness(ic) < unifrnd
                                        !boundedRegionIsTooLarge = logLikeFitness < unifrnd
                                    end if
                                end if

                                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                ! Decide on the fate of the partition.
                                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                if (boundedRegionIsTooLarge) then
                                    exit loopCurrentPartitionOptimization
                                else
                                    Partition%LogVolNormed(ic) = logVolNormed
                                    Partition%EffectiveSize(ic) = effectiveSize
                                    Partition%LogLikeFitness(ic) = logLikeFitness
                                    Partition%Center(1:Partition%nd,ic) = Center(1:Partition%nd)
                                    if (logVolRatio < 0._RK) then
                                        Partition%LogLikeFitness(ic) = 0._RK
                                        scaleFactorSq = scaleFactorSq * exp(ndHalfInverse * logVolRatio)
                                        Partition%LogVolNormed(ic) = Partition%LogVolNormed(ic) - logVolRatio
                                    end if
                                    scaleFactor = sqrt(scaleFactorSq)
                                    scaleFactorSqInverse = 1._RK / scaleFactorSq
                                    do j = 1, Partition%nd
                                        Partition%ChoDia(j,ic) = ChoDia(j) * scaleFactor
                                        !Partition%ChoLowCovUpp(1:j,j,ic) = ChoLowCovUpp(1:j,j) * scaleFactorSq
                                        Partition%InvCovMat(1:Partition%nd,j,ic) = InvCovMat(1:Partition%nd,j) * scaleFactorSqInverse
                                        do i = j+1, Partition%nd
                                            Partition%ChoLowCovUpp(i,j,ic) = ChoLowCovUpp(i,j) * scaleFactor
                                        end do
                                    end do
                                    MahalSq(MahalSqSortedIndex,ic) = MahalSqSinglePartition * scaleFactorSqInverse
                                    cycle loopCurrentPartitionOptimization
                                end if

                                end block blockTempShapeOptimization

                            else blockShapeOptimizationEnabled

                                exit loopCurrentPartitionOptimization

                            end if blockShapeOptimizationEnabled

                            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        end do loopCurrentPartitionOptimization

                    end if blockNormalPartition

                end do loopOverPartition

            end if blockOptimizationEnabled
#endif
        end if blockPartitionPostProcessing

#if defined DEBUG_ENABLED || defined TESTING_ENABLED || defined CODECOVE_ENABLED
        block
            integer :: ipstart, ipend
            ipend = 0
            do ic = 1, Partition%neopt
                ipstart = ipend + 1
                ipend = ipstart + Partition%Size(ic) - 1
                if (any(Partition%Membership(ipstart:ipend) /= Partition%Membership(ipstart))) then
                    write(*,*) PROCEDURE_NAME//": Internal error occurred: Partition%Membership is not unique:"
                    write(*,*) "Membership:", Partition%Membership(ipstart:ipend)
                    write(*,*) "ic, ipstart, ipend:", ic, ipstart, ipend
                    error stop
                end if
            end do
        end block
#endif

    end subroutine runPartition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    ! Recursively partition the input `Point`.
    !>
    !> \warning
    !> if `parentLogVolNormed > NEGINF_RK` , then `pointLogVolNormed > NEGINF_RK` must also hold.
    !>
    !> \warning
    !> if `mahalSqWeightExponent /= 0._RK` , then `pointLogVolNormed > NEGINF_RK` and `parentLogVolNormed > NEGINF_RK` must also hold.
    recursive subroutine runRecursivePartition  ( nd & ! LCOV_EXCL_LINE
                                                , np & ! LCOV_EXCL_LINE
                                                , nc & ! LCOV_EXCL_LINE
                                                , nt & ! LCOV_EXCL_LINE
                                                , npc & ! LCOV_EXCL_LINE
                                                , nps & ! LCOV_EXCL_LINE
                                                , npe & ! LCOV_EXCL_LINE
                                                , nsim & ! LCOV_EXCL_LINE
                                                , nemax & ! LCOV_EXCL_LINE
                                                , minSize & ! LCOV_EXCL_LINE
                                                , kmeansRelTol & ! LCOV_EXCL_LINE
                                                , logExpansion & ! LCOV_EXCL_LINE
                                                , logShrinkage & ! LCOV_EXCL_LINE
                                                , inclusionFraction & ! LCOV_EXCL_LINE
                                                , pointLogVolNormed & ! LCOV_EXCL_LINE
                                                , parentLogVolNormed & ! LCOV_EXCL_LINE
                                                , mahalSqWeightExponent & ! LCOV_EXCL_LINE
                                                , maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                                , maxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                                                , maxAllowedKvolumeRecursion & ! LCOV_EXCL_LINE
                                                , Point & ! LCOV_EXCL_LINE
                                                , PartitionSize & ! LCOV_EXCL_LINE
                                                , neopt & ! LCOV_EXCL_LINE
                                                , Center & ! LCOV_EXCL_LINE
                                                , ChoDia & ! LCOV_EXCL_LINE
#if defined MAXDEN
                                                , LogLikeFitness & ! LCOV_EXCL_LINE
                                                , EffectiveSize & ! LCOV_EXCL_LINE
                                                , MahalSq & ! LCOV_EXCL_LINE
#endif
                                                , InvCovMat & ! LCOV_EXCL_LINE
                                                , Membership & ! LCOV_EXCL_LINE
                                                , PointIndex & ! LCOV_EXCL_LINE
                                                , ChoLowCovUpp & ! LCOV_EXCL_LINE
                                                , LogVolNormed & ! LCOV_EXCL_LINE
                                                , numRecursiveCall & ! LCOV_EXCL_LINE
                                                , convergenceFailureCount & ! LCOV_EXCL_LINE
                                                , PointAvg & ! LCOV_EXCL_LINE
                                                , PointStd & ! LCOV_EXCL_LINE
                                                )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: runRecursivePartition
#endif
!#if defined DEBUG_ENABLED
!        use, intrinsic :: ieee_arithmetic, only: ieee_support_underflow_control
!        use, intrinsic :: ieee_arithmetic, only: ieee_set_underflow_mode
!#endif
        use Constants_mod, only: IK, RK, NEGINF_RK !, POSINF_RK !, SPR
        use Matrix_mod, only: getCholeskyFactor, getInvMatFromCholFac, getEye
        use Kmeans_mod, only: Kmeans_type
        use Math_mod, only: getLogSumExp

        implicit none

        integer , intent(in)                    :: nd, np, nc, nt
        integer , intent(in)                    :: npc, nps, npe
        integer , intent(in)                    :: nsim, nemax, minSize
        real(RK)    , intent(in)                :: logExpansion
        real(RK)    , intent(in)                :: logShrinkage
        real(RK)    , intent(in)                :: kmeansRelTol
        real(RK)    , intent(in)                :: inclusionFraction
        real(RK)    , intent(in)                :: pointLogVolNormed
        real(RK)    , intent(in)                :: parentLogVolNormed
        real(RK)    , intent(in)                :: mahalSqWeightExponent
        integer , intent(in)                    :: maxAllowedKmeansFailure
        integer , intent(in)                    :: maxAllowedKmeansRecursion
        integer , intent(in)                    :: maxAllowedKvolumeRecursion
        integer , intent(out)                   :: neopt
        integer , intent(out)                   :: PartitionSize(nemax)
        integer , intent(inout)                 :: numRecursiveCall
        integer , intent(inout)                 :: convergenceFailureCount
        integer , intent(out)                   :: Membership(nps:npe)
        integer , intent(inout)                 :: PointIndex(nps:npe)
        real(RK)    , intent(inout)             :: Point(nd,np)
        real(RK)    , intent(inout)             :: Center(nd,nemax)
        real(RK)    , intent(inout)             :: ChoLowCovUpp(nd,nd,nemax)
        real(RK)    , intent(inout)             :: InvCovMat(nd,nd,nemax)
        real(RK)    , intent(inout)             :: LogVolNormed(nemax)
#if defined MAXDEN
        real(RK)    , intent(inout)             :: LogLikeFitness(nemax)
        integer , intent(out)                   :: EffectiveSize(nemax)
        real(RK)    , intent(out)               :: MahalSq(np,nemax)
#endif
        real(RK)    , intent(inout)             :: ChoDia(nd,nemax)
        real(RK)    , intent(inout) , optional  :: PointAvg(nd) ! must coexist with PointStd
        real(RK)    , intent(inout) , optional  :: PointStd(nd)

        type(Kmeans_type)                       :: Kmeans
        real(RK)    , allocatable               :: InvStd(:)
        real(RK)                                :: ndInverse
        real(RK)                                :: NormedPoint(nd)
        real(RK)                                :: MahalSqWeight(nc)
        real(RK)                                :: StanPoint(nd,nps:npe)
        real(RK)                                :: KmeansLogVolExpected(nc) ! New Cluster Volume estimates.
        real(RK)                                :: pointLogVolNormedDefault
        real(RK)                                :: scaleFactorSqInverse
        real(RK)                                :: scaleFactor
        integer                                 :: KmeansNemax(nc)      ! maximum allowed number of clusters in each of the two K-means clusters.
        integer                                 :: KmeansNeopt(0:nc)    ! optimal number of clusters in each of the child clusters.
        integer                                 :: KmeansMemberCounter(nc)
        integer                                 :: minSizePartition
        integer                                 :: i, j, ip, jp, ic, jc, ipstart, ipend
        integer                                 :: icstart, icend, recursionCounter
        integer                                 :: counterEffectiveSize
        integer                                 :: icmin, nemaxRemained
        logical                                 :: subclusterSearchWarranted
        logical                                 :: boundedRegionIsTooLarge
        logical                                 :: mahalSqWeightEnabled
        logical                                 :: reclusteringNeeded
        real(RK)                                :: maxLogVolNormed
#if defined MAXDEN
        real(RK)                                :: unifrnd
        real(RK)                                :: logVolRatio
       !real(RK)                                :: KmeansMahalSq(np,nc)
        real(RK)                                :: scaleFactorSqInverse4KPM ! for MahalSq component of Kmeans Prop
        real(RK)                                :: parEllLogLike    ! recovery
        real(RK)                                :: logLikeSum       ! recovery
#endif
        real(RK)                                :: ParEllCenter(nd) ! recovery
        real(RK)                                :: ParEllChoDia(1:nd) ! recovery
        real(RK)                                :: parEllLogVolNormed ! recovery
        real(RK)                                :: ParEllInvCovMat(1:nd,1:nd) ! recovery
        real(RK)                                :: ParEllChoLowCovUpp(1:nd,1:nd) ! recovery

        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME//"@runRecursivePartition()"

        ndInverse = 1._RK / nd
        minSizePartition = max(1,minSize)
        numRecursiveCall = numRecursiveCall + 1

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Perform the initial Kmeans clustering.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        blockKmeansClustering: if (present(PointStd)) then

            InvStd = 1._RK / PointStd(1:nd)
            do concurrent(ip = nps:npe)
                StanPoint(1:nd,ip) = (Point(1:nd,ip) - PointAvg(1:nd)) * InvStd(1:nd)
            end do

            Kmeans = Kmeans_type( nd = nd & ! LCOV_EXCL_LINE
                                , np = npc & ! LCOV_EXCL_LINE
                                , nc = nc & ! LCOV_EXCL_LINE
                                , nt = nt & ! LCOV_EXCL_LINE
                                , Point = StanPoint & ! LCOV_EXCL_LINE
                                , relTol = kmeansRelTol & ! LCOV_EXCL_LINE
                                , minSize = minSize & ! LCOV_EXCL_LINE
                                , nfailMax = maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                , niterMax = maxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                                )
            if (Kmeans%Err%occurred) then
                PartitionSize(1) = npc
                neopt = 1
                return
            end if

            ! Shift centers back.

            do concurrent(ic = 1:nc)
                Kmeans%Center(1:nd,ic) = Kmeans%Center(1:nd,ic) * PointStd(1:nd) + PointAvg(1:nd)
            end do

        else blockKmeansClustering

            Kmeans = Kmeans_type( nd = nd & ! LCOV_EXCL_LINE
                                , np = npc & ! LCOV_EXCL_LINE
                                , nc = nc & ! LCOV_EXCL_LINE
                                , nt = nt & ! LCOV_EXCL_LINE
                                , Point = Point(1:nd,nps:npe) & ! LCOV_EXCL_LINE
                                , relTol = kmeansRelTol & ! LCOV_EXCL_LINE
                                , minSize = minSize & ! LCOV_EXCL_LINE
                                , nfailMax = maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                , niterMax = maxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                                )
            if (Kmeans%Err%occurred) then
                PartitionSize(1) = npc
                neopt = 1
                return
            end if

        end if blockKmeansClustering

#if (defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED)
        block
            integer :: sizeDum
            real(RK) :: CenterDum(nd)
            do ic = 1, nc
                sizeDum = 0
                CenterDum = 0._RK
                do ip = nps, npe
                    jp = ip - nps + 1
                    if (Kmeans%Membership(jp) == ic) then
                        sizeDum = sizeDum + 1
                        CenterDum(1:nd) = CenterDum(1:nd) + Point(1:nd,ip)
                    end if
                end do
                CenterDum(1:nd) = CenterDum(1:nd) / sizeDum
                if (any(abs(CenterDum(1:nd) - Kmeans%Center(1:nd,ic)) > 1.e-6_RK)) then
                    if (sizeDum /= Kmeans%Size(ic)) then; write(*,*) "sizeDum, Kmeans%Size = ", sizeDum, Kmeans%Size; error stop; endif
                    write(*,*) "Kmeans%Center = ", Kmeans%Center(1:nd,ic)
                    write(*,*) "CenterDum     = ", CenterDum(1:nd)
                    write(*,*) "ic = ", ic
                endif
            end do
        end block
#endif

        ! Compute cluster properties.

        allocate(Kmeans%Prop%ChoDia         (nd,nc))
        allocate(Kmeans%Prop%MahalSq        (np,nc))
        allocate(Kmeans%Prop%InvCovMat      (nd,nd,nc))
        allocate(Kmeans%Prop%ChoLowCovUpp   (nd,nd,nc))
        allocate(Kmeans%Prop%EffectiveSize  (nc))
        allocate(Kmeans%Prop%ScaleFactorSq  (nc))
        allocate(Kmeans%Prop%LogVolNormed   (nc))
        allocate(Kmeans%Prop%CumSumSize     (0:nc))
        allocate(Kmeans%Prop%Index          (nps:npe))

#include "Partition_mod.prop.inc.f90"

!        call Kmeans%getProp ( nd = nd & ! LCOV_EXCL_LINE
!                            , np = npc & ! LCOV_EXCL_LINE
!                            , Point = Point(1:nd,nps:npe) & ! LCOV_EXCL_LINE
!                            , Index = PointIndex & ! LCOV_EXCL_LINE
!#if defined MINVOL
!                            , inclusionFraction = inclusionFraction & ! LCOV_EXCL_LINE
!#endif
!                            , pointLogVolNormed = pointLogVolNormed & ! LCOV_EXCL_LINE
!                            )
!        if (Kmeans%Err%occurred) then
!            ! LCOV_EXCL_START
!            PartitionSize(1) = npc
!            neopt = 1
!            return
!            ! LCOV_EXCL_STOP
!        end if

#if (defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED)
        if (any(abs(Kmeans%Prop%LogVolNormed(1:nc)) < 1.e-6_RK)) then ! xxx this must be removed as this condition could normally happen, even though rarely.
            write(*,*) "Internal error occurred #1: Kmeans%Prop%LogVolNormed == 0. : Kmeans%Prop%LogVolNormed :"
            write(*,*) Kmeans%Prop%LogVolNormed
            error stop
        end if
        loopComputeScaleFactor: do ic = 1, nc
                block
                    use Statistics_mod, only: getSamCholFac, getMahalSq
                    real(RK), parameter :: TOLERANCE = 1.e-2_RK
                    real(RK) :: mahalSqScalar
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

                    do ip = 1, npc
                        NormedPoint(1:nd) = Point(1:nd,ip+nps-1) - Kmeans%Center(1:nd,ic)
                        mahalSqScalar = dot_product( NormedPoint , matmul(Kmeans%Prop%InvCovMat(1:nd,1:nd,ic), NormedPoint) )
                        if (mahalSqScalar<0._RK) then
                        ! LCOV_EXCL_START
                            mahalSqScalar = -1._RK
                            write(*,*) PROCEDURE_NAME
                            write(*,*) "mahalSqScalar<0._RK", ip, mahalSqScalar
                            error stop
                        end if
                        ! LCOV_EXCL_STOP
                        dum = getMahalSq( nd, Kmeans%Center(1:nd,ic), Kmeans%Prop%InvCovMat(1:nd,1:nd,ic), Point(1:nd,ip+nps-1) )
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
                end block
                !MahalSqScalar(1:np,ic) = MahalSqScalar(1:np,ic) / Kmeans%Prop%ScaleFactorSq(ic)
            end do loopComputeScaleFactor
#endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! If optimization is enabled (i.e., not kmeans), reassign points to clusters based on minimum-volume, maximum-density, ...
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !KmeansBackup = Kmeans_type(nd, npc, nc)
        !call KmeansBackup%Prop%allocate(nd, npc, nc)
        !KmeansBackup%Prop%LogVolNormed = 2 * Kmeans%Prop%LogVolNormed

        blockRecursivePartitionOptimization: if (maxAllowedKvolumeRecursion > 0) then

#if defined MAXDEN
            mahalSqWeightEnabled = mahalSqWeightExponent /= 0._RK
#elif defined MINVOL
            mahalSqWeightEnabled = mahalSqWeightExponent /= 0._RK .and. parentLogVolNormed > NEGINF_RK
#endif

            loopRecursivePartitionOptimization: do recursionCounter = 1, maxAllowedKvolumeRecursion

                ! Ensure the stability of the minimum bounding regions by making a backup.

                !if (Kmeans%Prop%logSumVolNormed < KmeansBackup%Prop%logSumVolNormed) then ! Create a backup of the current clustering properties. This must never occur on the first try.
                !    !KmeansBackup = Kmeans
                !    BackupPoint(1:nd,1:npc)                          = Point(1:nd,1:npc)
                !    BackupPointIndex(1:npc)                          = PointIndex(1:npc)
                !    KmeansBackup%Membership(1:npc)                   = Kmeans%Membership(1:npc)
                !    KmeansBackup%Center(1:nd,1:nc)                  = Kmeans%Center(1:nd,1:nc)
                !    KmeansBackup%Prop%ChoDia(1:nd,1:nc)             = Kmeans%Prop%ChoDia(1:nd,1:nc)
                !    KmeansBackup%Prop%MahalSq(1:npc,1:nc)            = Kmeans%Prop%MahalSq(1:npc,1:nc)
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
                !    Point(1:nd,nps:npe)                             = BackupPoint(1:nd,nps:npe)
                !    PointIndex(1:npc)                                = BackupPointIndex(1:npc)
                !    Kmeans%Membership(1:npc)                         = KmeansBackup%Membership(1:npc)
                !    Kmeans%Center(1:nd,1:nc)                        = KmeansBackup%Center(1:nd,1:nc)
                !    Kmeans%Prop%ChoDia(1:nd,1:nc)                   = KmeansBackup%Prop%ChoDia(1:nd,1:nc)
                !    Kmeans%Prop%MahalSq(1:npc,1:nc)                  = KmeansBackup%Prop%MahalSq(1:npc,1:nc)
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

                ! Convert Mean(Point) to sum(Point) per cluster and compute variational weights.

                do concurrent(ic = 1:nc)
                    Kmeans%Center(1:nd,ic) = Kmeans%Center(1:nd,ic) * Kmeans%Size(ic)
                end do

                ! Determine the cluster to which the point is the closest and switch the point membership if needed.

                reclusteringNeeded = .false.
                if (mahalSqWeightEnabled) then

                    do concurrent(ic = 1:nc)
                        if (Kmeans%Size(ic) > 0) then
                            KmeansLogVolExpected(ic) = pointLogVolNormed + log(real(Kmeans%Prop%EffectiveSize(ic),RK)) ! @attn: Kmeans%Prop%EffectiveSize(ic) >= Kmeans%Size(ic) > 0
                            MahalSqWeight(ic) = exp(mahalSqWeightExponent * (Kmeans%Prop%LogVolNormed(ic) - KmeansLogVolExpected(ic)))
                        end if
                    end do

                    loopReassignWeightedPoint: do ip = nps, npe
                        jp = ip - nps + 1
                        icmin = Kmeans%Membership(jp)
                        do ic = 1, nc
                            if (Kmeans%Size(ic) > 0 .and. ic /= icmin .and. MahalSqWeight(ic) * Kmeans%Prop%MahalSq(ip,ic) < MahalSqWeight(icmin) * Kmeans%Prop%MahalSq(ip,icmin)) icmin = ic
                       end do
                        if (icmin /= Kmeans%Membership(jp)) then
#if (defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED)
                            if (Kmeans%Size(Kmeans%Membership(jp)) <= 0) then
                                write(*,*) PROCEDURE_NAME//": Internal error occurred: Kmeans%Size(Kmeans%Membership(jp)) <= 0"
                                write(*,*) "jp, Kmeans%Membership(jp), Kmeans%Size(Kmeans%Membership(jp)): ", jp, Kmeans%Membership(jp), Kmeans%Size(Kmeans%Membership(jp))
                                error stop
                            end if
#endif
                            Kmeans%Size(icmin) = Kmeans%Size(icmin) + 1
                            Kmeans%Size(Kmeans%Membership(jp)) = Kmeans%Size(Kmeans%Membership(jp)) - 1
                            Kmeans%Center(1:nd,icmin) = Kmeans%Center(1:nd,icmin) + Point(1:nd,ip)
                            Kmeans%Center(1:nd,Kmeans%Membership(jp)) = Kmeans%Center(1:nd,Kmeans%Membership(jp)) - Point(1:nd,ip)
                            Kmeans%Membership(jp) = icmin
                            reclusteringNeeded = .true.
                        end if
                    end do loopReassignWeightedPoint

                else

                    loopReassignPoint: do ip = nps, npe
                        jp = ip - nps + 1
                        icmin = Kmeans%Membership(jp)
                        do ic = 1, nc
                            if (Kmeans%Size(ic) > 0 .and. ic /= icmin .and. Kmeans%Prop%MahalSq(ip,ic) < Kmeans%Prop%MahalSq(ip,icmin)) icmin = ic
                        end do
                        if (icmin /= Kmeans%Membership(jp)) then ! .and. Kmeans%Size(Kmeans%Membership(jp)) > minSize) then
#if (defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED)
                            if (Kmeans%Size(Kmeans%Membership(jp)) <= 0) then
                                write(*,*) PROCEDURE_NAME//": Internal error occurred: Kmeans%Size(Kmeans%Membership(jp)) <= 0"
                                write(*,*) "jp, Kmeans%Membership(jp), Kmeans%Size(Kmeans%Membership(jp)): ", jp, Kmeans%Membership(jp), Kmeans%Size(Kmeans%Membership(jp))
                                error stop
                            end if
#endif
                            Kmeans%Size(icmin) = Kmeans%Size(icmin) + 1
                            Kmeans%Size(Kmeans%Membership(jp)) = Kmeans%Size(Kmeans%Membership(jp)) - 1
                            Kmeans%Center(1:nd,icmin) = Kmeans%Center(1:nd,icmin) + Point(1:nd,ip)
                            Kmeans%Center(1:nd,Kmeans%Membership(jp)) = Kmeans%Center(1:nd,Kmeans%Membership(jp)) - Point(1:nd,ip)
                            Kmeans%Membership(jp) = icmin
                            reclusteringNeeded = .true.
                        end if
                    end do loopReassignPoint

                end if

#if (defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED)
                !if (any(Kmeans%Size<minSize)) then
                !    !write(*,*) "Kmeans%Size < minSize : minSize, Kmeans%Size = ", minSize, Kmeans%Size
                !    !error stop ! @attn: This is NOT an error anymore as Size can become less than minSize.
                !end if
#endif

                ! Convert sum(Point) back to center of clusters.

                do concurrent(ic = 1:nc) ! xxx revive this line
                    !write(*,*) "Kmeans%Center(1:nd,ic), Kmeans%Size", Kmeans%Center(1:nd,ic), Kmeans%Size ! xxx
                    if (Kmeans%Size(ic) > 0) Kmeans%Center(1:nd,ic) = Kmeans%Center(1:nd,ic) / Kmeans%Size(ic)
                    !if (any(abs(Kmeans%Center(:,ic) - 0.677385375580642_RK)<1.e-6_RK)) then
                    !    write(*,*) "Kmeans%Center(nd,ic) = ", Kmeans%Center(:,ic)
                    !    write(*,*) "0.677385375580642 = ", 0.677385375580642_RK
                    !    write(*,*) "Kmeans%Size(ic) = ", Kmeans%Size(ic)
                    !    write(*,*) "reclusteringNeeded = ", reclusteringNeeded
                    !    write(*,*) "ic = ", ic
                    !    block
                    !        integer :: sizeDum
                    !        real(RK) :: CenterDum(nd)
                    !        sizeDum = 0
                    !        CenterDum = 0._RK
                    !        do ip = nps, npe
                    !            jp = ip - nps + 1
                    !            if (Kmeans%Membership(jp) == 2) then
                    !                sizeDum = sizeDum + 1
                    !                CenterDum(1:nd) = CenterDum(1:nd) + Point(1:nd,ip)
                    !            end if
                    !        end do
                    !        write(*,*) "sizeDum = ", sizeDum
                    !        write(*,*) "CenterDum(1:nd) = ", CenterDum(1:nd) / sizeDum
                    !    end block
                    !end if
                end do

                ! Restart the process if anything has changed.

                blockReclusteringNeeded: if (reclusteringNeeded) then ! perform reassignment

                    ! Reorder Point based on the identified clusters and recompute the cluster properties.

#include "Partition_mod.prop.inc.f90"

#if (defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED)
                    block
                        real(RK), allocatable :: LogVolNormed(:)
                        if (any(abs(Kmeans%Prop%LogVolNormed(1:nc)) < 1.e-6_RK)) then ! xxx : this check must be removed as the condition could normally though rarely be true.
                            write(*,*) "Internal error occurred #2 : Kmeans%Prop%LogVolNormed == 0. : Kmeans%Prop%LogVolNormed :"
                            write(*,*) Kmeans%Prop%LogVolNormed
                            error stop
                        end if
                    end block
#endif
                    cycle loopRecursivePartitionOptimization

                end if blockReclusteringNeeded

                exit loopRecursivePartitionOptimization ! mission accomplished, current volumes minimized.

            end do loopRecursivePartitionOptimization

        else blockRecursivePartitionOptimization

            mahalSqWeightEnabled = .false.

        end if blockRecursivePartitionOptimization

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Rescale the volumes if needed and for as long as needed, based on the user-input volume estimate
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#if defined MINVOL
        maxLogVolNormed = NEGINF_RK ! indicator for rescaling.
        if (parentLogVolNormed > NEGINF_RK) then

            do ic = 1, nc

                !KmeansLogVolExpected(ic) = pointLogVolNormed + log(real(Kmeans%Prop%EffectiveSize(ic),RK)) ! @attn: xxx how is Kmeans%Prop%EffectiveSize(ic) > 0 ensured?
                if (Kmeans%Size(ic) > 0) then

                    if (.not. mahalSqWeightEnabled) KmeansLogVolExpected(ic) = pointLogVolNormed + log(real(Kmeans%Prop%EffectiveSize(ic),RK)) ! @attn: Kmeans%Prop%EffectiveSize(ic) >= Kmeans%Size(ic) > 0

                    if (Kmeans%Prop%LogVolNormed(ic) < KmeansLogVolExpected(ic)) then

                        Kmeans%Prop%LogVolNormed(ic) = KmeansLogVolExpected(ic)
                        if (Kmeans%Prop%LogVolNormed(ic) > maxLogVolNormed) maxLogVolNormed = Kmeans%Prop%LogVolNormed(ic)

                        ! Rescale the Cholesky factor and its diagonals such that they describe the bounded shape of the cluster
                        ! @todo: xxx This will have to be moved to a later location for better efficiency.
                        ! pay attention to the correct value of Kmeans%Prop%ScaleFactorSq(ic).

                        scaleFactor = exp( ndInverse * (KmeansLogVolExpected(ic) - Kmeans%Prop%LogVolNormed(ic)) )
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
                elseif (Kmeans%Size(ic) == 0 .and. Kmeans%Prop%LogVolNormed(ic) /= NEGINF_RK) then
                    write(*,*) PROCEDURE_NAME//": Kmeans%Size(ic) == 0 .and. Kmeans%Prop%LogVolNormed(ic) /= NEGINF_RK", ic, Kmeans%Prop%LogVolNormed(ic)
                    error stop
#endif

                end if

            end do
            boundedRegionIsTooLarge = LogVolNormed(1) > logExpansion + parentLogVolNormed

        else

            KmeansLogVolExpected(1:nc) = NEGINF_RK
            boundedRegionIsTooLarge = .false.

        end if

        if (maxLogVolNormed > NEGINF_RK) Kmeans%Prop%logSumVolNormed = getLogSumExp(nc, Kmeans%Prop%LogVolNormed, maxLogVolNormed)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Check if further clustering is warranted.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        blockSubclusterSearch: if (LogVolNormed(1) > logShrinkage + Kmeans%Prop%logSumVolNormed .or. boundedRegionIsTooLarge) then
#endif

#if defined MAXDEN
            if (logShrinkage > NEGINF_RK) then
                parEllLogLike = LogLikeFitness(1) ! recovery
#endif
                ParEllCenter(1:nd) = Center(1:nd,1) ! recovery
                ParEllChoDia(1:nd) = ChoDia(1:nd,1) ! recovery
                parEllLogVolNormed = LogVolNormed(1) ! recovery
                ParEllInvCovMat(1:nd,1:nd) = InvCovMat(1:nd,1:nd,1) ! recovery
                ParEllChoLowCovUpp(1:nd,1:nd) = ChoLowCovUpp(1:nd,1:nd,1) ! recovery
#if defined MAXDEN
            end if
#endif

            ! At least nc sub-clusters is better than one cluster. Now search for more sub-sub-clusters.

            nemaxRemained = nemax
            KmeansNeopt(0) = 0
            icstart = 1

            ! Search for more subclusters.

            loopSubclusterSearch: do ic = 1, nc

                nemaxRemained = nemaxRemained - KmeansNeopt(ic-1)
                KmeansNemax(ic) = min(nemaxRemained - nc + ic, Kmeans%Size(ic) / minSizePartition)
                subclusterSearchWarranted = KmeansNemax(ic) > 1

#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
                if (nemaxRemained < KmeansNemax(ic)) then
                    write(*,"(*(g0,:,', '))") PROCEDURE_NAME//": nemaxRemained < KmeansNemax(ic) : numRecursiveCall, np, ic, minSize, nemaxRemained, KmeansNemax(ic), nemax" &
                                            , numRecursiveCall, npc, ic, minSize, nemaxRemained, KmeansNemax(ic), nemax
                    error stop
                end if
#endif

                ipstart = Kmeans%Prop%CumSumSize(ic-1) + 1
                ipend = Kmeans%Prop%CumSumSize(ic)

                icstart = icstart + KmeansNeopt(ic-1)
                icend = icstart + KmeansNemax(ic) - 1

#if defined MAXDEN
                !scaleFactorSqInverse = 1._RK
                if (Kmeans%size(ic) > 0) then
!                    if (inclusionFraction > 0._RK) then
!                        counterEffectiveSize = 0
!                        do ip = 1, nps - 1
!                            NormedPoint(1:nd) = Point(1:nd,ip) - Kmeans%Center(1:nd,ic)
!                            Kmeans%Prop%MahalSq(ip,ic) = dot_product( NormedPoint(1:nd) , matmul(Kmeans%Prop%InvCovMat(1:nd,1:nd,ic), NormedPoint(1:nd)) )
!                            if (KmeansMahalSq(ip,ic) <= 1._RK) counterEffectiveSize = counterEffectiveSize + 1
!                        end do
!#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
!                        !EffectiveSizeOld(ic) = Kmeans%Prop%EffectiveSize(ic)
!                        if (npc /= npe - nps + 1) then
!                            write(*,*) PROCEDURE_NAME//": Internal error occurred : npc /= npe - nps + 1"
!                            error stop
!                        end if
!                        block
!                            real(RK) :: mahalSqScalar
!                            do ip = ipstart + nps - 1, ipend + nps - 1
!                                NormedPoint(1:nd) = Point(1:nd,ip) - Kmeans%Center(1:nd,ic)
!                                mahalSqScalar = dot_product( NormedPoint(1:nd) , matmul(Kmeans%Prop%InvCovMat(1:nd,1:nd,ic), NormedPoint(1:nd)) )
!                                if (mahalSqScalar - 1._RK > 1.e-6_RK) then
!                                    write(*,*) "Internal error occurred : mahalSqScalar > 1 :", mahalSqScalar
!                                    write(*,*) "nps, ip, npe, np", nps, ip, npe, np
!                                    write(*,*) "ic, Kmeans%size(ic)", ic, Kmeans%size(ic)
!                                    write(*,*) "numRecursiveCall", numRecursiveCall
!                                    error stop
!                                endif
!                            end do
!                        end block
!#endif
!                        do ip = npe + 1, np
!                            NormedPoint(1:nd) = Point(1:nd,ip) - Kmeans%Center(1:nd,ic)
!                            KmeansMahalSq(ip,ic) = dot_product( NormedPoint(1:nd) , matmul(Kmeans%Prop%InvCovMat(1:nd,1:nd,ic), NormedPoint(1:nd)) )
!                            if (KmeansMahalSq(ip,ic) <= 1._RK) counterEffectiveSize = counterEffectiveSize + 1
!                        end do
!                        Kmeans%Prop%EffectiveSize(ic)   = Kmeans%Size(ic) & ! LCOV_EXCL_LINE
!                                                        + nint(inclusionFraction * & ! LCOV_EXCL_LINE
!                                                        ( count(Kmeans%Prop%MahalSq(1:ipstart-1,ic)<=Kmeans%Prop%ScaleFactorSq(ic)) & ! LCOV_EXCL_LINE
!                                                        + count(Kmeans%Prop%MahalSq(ipend+1:npc,ic)<=Kmeans%Prop%ScaleFactorSq(ic)) & ! LCOV_EXCL_LINE
!                                                        + counterEffectiveSize & ! LCOV_EXCL_LINE
!                                                        ), kind = IK)
!                    end if
                    KmeansLogVolExpected(ic) = pointLogVolNormed + log(real(Kmeans%Prop%EffectiveSize(ic),RK)) ! @attn: Kmeans%Prop%EffectiveSize(ic) >= Kmeans%Size(ic) > 0
                    logVolRatio = Kmeans%Prop%LogVolNormed(ic) - KmeansLogVolExpected(ic)
                    LogLikeFitness(icstart) = getLogProbPoisDiff(Kmeans%Prop%EffectiveSize(ic), logVolRatio)
#if defined DEBUG_ENABLED || defined TESTING_ENABLED || defined CODECOVE_ENABLED
                    !if (LogLikeFitness(icstart) < 0._RK .or. LogLikeFitness(icstart) > 1._RK) then
                    if (LogLikeFitness(icstart) < NEGINF_RK .or. LogLikeFitness(icstart) > 0._RK) then
                        write(*,*) "Internal error occurred: LogLikeFitness(icstart) < NEGINF_RK .or. LogLikeFitness(icstart) > 0._RK"
                        write(*,*) "Kmeans%Prop%LogVolNormed(ic), KmeansLogVolExpected(ic), LogLikeFitness(icstart):" & ! LCOV_EXCL_LINE
                                 , Kmeans%Prop%LogVolNormed(ic), KmeansLogVolExpected(ic), LogLikeFitness(icstart)
                        error stop
                    end if
#endif
                    boundedRegionIsTooLarge = logVolRatio > 0._RK
                    if (boundedRegionIsTooLarge) then
                        call random_number(unifrnd)
                        if (unifrnd > 0._RK) then; unifrnd = log(unifrnd); else; unifrnd = -huge(unifrnd); end if ! LCOV_EXCL_LINE or unifrnd = log(max(tiny(unifrnd),unifrnd))
                        boundedRegionIsTooLarge = LogLikeFitness(icstart) < unifrnd
                    end if
                    subclusterSearchWarranted = subclusterSearchWarranted .and. boundedRegionIsTooLarge
                end if
#endif
                blockSubclusterCheck: if (subclusterSearchWarranted) then ! .and. Kmeans%Prop%LogVolNormed(ic)>1.1_RK*KmeansLogVolExpected(ic)) then

                    LogVolNormed(icstart) = Kmeans%Prop%LogVolNormed(ic)

                    if (present(PointStd)) then

                        PointAvg(1:nd) = Kmeans%Center(1:nd,ic)
                        do concurrent(i = 1:nd)
                            PointStd(i) = sqrt(Kmeans%Prop%ChoLowCovUpp(i,i,ic))
                        end do

                        !do concurrent(ip = ipstart+nps-1:ipend+nps-1)
                        !    Point(1:nd,ip) = Point(1:nd,ip) - PointAvg(1:nd)
                        !end do

                    end if

                    call runRecursivePartition  ( nd = nd & ! LCOV_EXCL_LINE
                                                !, np = Kmeans%Size(ic) & ! LCOV_EXCL_LINE
                                                , np = np & ! LCOV_EXCL_LINE
                                                , nc = min(nc, KmeansNemax(ic)) & ! LCOV_EXCL_LINE
                                                , nt = nt & ! LCOV_EXCL_LINE
                                                , npc = Kmeans%Size(ic) & ! LCOV_EXCL_LINE
                                                , nps = ipstart & ! LCOV_EXCL_LINE
                                                , npe = ipend & ! LCOV_EXCL_LINE
                                                , nsim = nsim & ! LCOV_EXCL_LINE
                                                , nemax = KmeansNemax(ic) & ! LCOV_EXCL_LINE
                                                , minSize = minSize & ! LCOV_EXCL_LINE
                                                , kmeansRelTol = kmeansRelTol & ! LCOV_EXCL_LINE
                                                , logExpansion = logExpansion & ! LCOV_EXCL_LINE
                                                , logShrinkage = logShrinkage & ! LCOV_EXCL_LINE
                                                , inclusionFraction = inclusionFraction & ! LCOV_EXCL_LINE
                                                , pointLogVolNormed = pointLogVolNormed & ! LCOV_EXCL_LINE
                                                , parentLogVolNormed = KmeansLogVolExpected(ic) & ! LCOV_EXCL_LINE
                                                , mahalSqWeightExponent = mahalSqWeightExponent & ! LCOV_EXCL_LINE
                                                , maxAllowedKmeansFailure = maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                                , maxAllowedKmeansRecursion = maxAllowedKmeansRecursion & ! LCOV_EXCL_LINE
                                                , maxAllowedKvolumeRecursion = maxAllowedKvolumeRecursion & ! LCOV_EXCL_LINE
                                                , Point = Point & ! (1:nd,ipstart:ipend) & ! LCOV_EXCL_LINE
                                                , PartitionSize = PartitionSize(icstart:icend) & ! LCOV_EXCL_LINE
                                                , neopt = KmeansNeopt(ic) & ! LCOV_EXCL_LINE
                                                , Center = Center(1:nd,icstart:icend) & ! LCOV_EXCL_LINE
                                                , ChoDia = ChoDia(1:nd,icstart:icend) & ! LCOV_EXCL_LINE
#if defined MAXDEN
                                                , LogLikeFitness = LogLikeFitness(icstart:icend) & ! LCOV_EXCL_LINE
                                                , EffectiveSize = EffectiveSize(icstart:icend) & ! LCOV_EXCL_LINE
                                                , MahalSq = MahalSq(1:np,icstart:icend) & ! LCOV_EXCL_LINE
#endif
                                                , InvCovMat = InvCovMat(1:nd,1:nd,icstart:icend) & ! LCOV_EXCL_LINE
                                                , Membership = Membership(ipstart:ipend) & ! LCOV_EXCL_LINE
                                                , PointIndex = PointIndex(ipstart:ipend) & ! LCOV_EXCL_LINE
                                                , ChoLowCovUpp = ChoLowCovUpp(1:nd,1:nd,icstart:icend) & ! LCOV_EXCL_LINE
                                                , LogVolNormed = LogVolNormed(icstart:icend) & ! LCOV_EXCL_LINE
                                                , numRecursiveCall = numRecursiveCall & ! LCOV_EXCL_LINE
                                                , convergenceFailureCount = convergenceFailureCount & ! LCOV_EXCL_LINE
                                                , PointAvg = PointAvg & ! LCOV_EXCL_LINE
                                                , PointStd = PointStd & ! LCOV_EXCL_LINE
                                                )

                    !if (KmeansNeopt(ic) > 1) then ! xxx remove
                        !if (present(PointStd)) then
                        !    do concurrent(jc = icstart : icstart + KmeansNeopt(ic) - 1)
                        !        Center(1:nd,jc) = Center(1:nd,jc) + Kmeans%Center(1:nd,ic)
                        !    end do
                        !end if
#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
                    if (any(abs(LogVolNormed(icstart:icstart+KmeansNeopt(ic)-1)) < 1.e-6_RK)) then ! xxx This must be removed later as the condition can be true in normal situations.
                        write(*,*) LogVolNormed(icstart:icstart+KmeansNeopt(ic)-1)
                        error stop
                    elseif (KmeansNeopt(ic) < 1) then
                        write(*,*) "KmeansNeopt(ic) < 1 : ", KmeansNeopt(ic), ic
                        error stop
                    end if
#endif

                    if (KmeansNeopt(ic) > 1) then
#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
                        if (any(Membership(ipstart:ipend) < 1) .or. any(Membership(ipstart:ipend) > KmeansNemax(ic))) then
                            write(*,*) PROCEDURE_NAME//": Internal error occurred."
                            write(*,*) "KmeansNemax(ic) : ", KmeansNemax(ic)
                            write(*,*) "Membership(ipstart:ipend) <= 0 : ", pack(Membership(ipstart:ipend), Membership(ipstart:ipend) <= 0)
                            write(*,*) "Membership(ipstart:ipend) > KmeansNemax(ic) : ", pack(Membership(ipstart:ipend), Membership(ipstart:ipend) > KmeansNemax(ic))
                            write(*,*) "Membership(ipstart:ipend) : ", Membership(ipstart:ipend)
                            write(*,*) "ic, icstart: ", ic, icstart
                            error stop
                        end if
#endif
                        Membership(ipstart:ipend) = Membership(ipstart:ipend) + icstart - 1
                    end if

                elseif (Kmeans%Size(ic) > 0) then blockSubclusterCheck

                    KmeansNeopt(ic) = 1

                else blockSubclusterCheck ! should happen only when Kmeans%Size(ic) == 0

                    KmeansNeopt(ic) = 0
#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
                if (Kmeans%Size(ic) /= 0) then
                    write(*,*) PROCEDURE_NAME//": Internal error occurred: Kmeans%Size(ic) /= 0: ic, Kmeans%Size(ic):", ic, Kmeans%Size(ic)
                    error stop
                end if
#endif

                end if blockSubclusterCheck

                blockSingleCluster: if (KmeansNeopt(ic) == 1) then ! .or. KmeansNemax(ic) == 1 ) then ! There is only one cluster here

                    Membership(ipstart:ipend) = icstart
                    PartitionSize(icstart) = Kmeans%Size(ic)
                    Center(1:nd,icstart) = Kmeans%Center(1:nd,ic)
#if defined MAXDEN
                    !if (LogLikeFitness(icstart) < unifrnd .and. logVolRatio < 0._RK) then ! enlarge
                    blockExpansionCheck: if (logVolRatio < 0._RK) then ! enlarge, no questions asked.

                        !LogLikeFitness(icstart) = 1._RK + LogLikeFitness(icstart)
                        LogLikeFitness(icstart) = 0._RK
                        LogVolNormed(icstart) = KmeansLogVolExpected(ic)
                        scaleFactor = exp( -ndInverse * logVolRatio )
                        scaleFactorSqInverse = 1._RK / scaleFactor**2
                        scaleFactorSqInverse4KPM = scaleFactorSqInverse / Kmeans%Prop%ScaleFactorSq(ic) ! assumes no MahalSq normalization in Kmeans%getProp()
                        MahalSq(1:np,icstart) = Kmeans%Prop%MahalSq(1:np,ic) * scaleFactorSqInverse4KPM
                        !MahalSq(1:nps-1,icstart) = KmeansMahalSq(1:nps-1,ic) * scaleFactorSqInverse
                        !MahalSq(nps:npe,icstart) = Kmeans%Prop%MahalSq(1:npc,ic) * scaleFactorSqInverse4KPM
                        !MahalSq(npe+1:np,icstart) = KmeansMahalSq(npe+1:np,ic) * scaleFactorSqInverse
                        EffectiveSize(icstart) = Kmeans%Size(ic) + nint(inclusionFraction * (count(MahalSq(1:ipstart-1,icstart) <= 1._RK) + count(MahalSq(ipend+1:np,icstart) <= 1._RK)))
#if (defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED)
                        !if (count(MahalSq(ipstart+nps-1:ipend+nps-1,icstart) <= 1._RK) < Kmeans%Size(ic)) then
                        !    write(*,*) PROCEDURE_NAME//": Internal error occurred : count(MahalSq(ipstart+nps-1:ipend+nps-1,icstart) <= 1._RK) >= Kmeans%Size(ic)"
                        !    write(*,*) "count(MahalSq(nps:npe,icstart) <= 1._RK) : ", count(MahalSq(ipstart+nps-1:ipend+nps-1,icstart) <= 1._RK)
                        !    write(*,*) "Kmeans%Size(ic) : ", Kmeans%Size(ic)
                        !    write(*,*) "logVolRatio < 0._RK :", logVolRatio < 0._RK
                        !    error stop
                        !else
                        !    block
                        !        integer :: extraSize
                        !        extraSize = nint(inclusionFraction * (count(MahalSq(1:ipstart+nps-2,icstart) <= 1._RK) + count(MahalSq(1:ipend+nps,icstart) <= 1._RK)))
                        !        if (Kmeans%Size(ic) + nint(inclusionFraction * extraSize) < Kmeans%Prop%EffectiveSize(ic)) then
                        !            write(*,*) PROCEDURE_NAME//": Internal error occurred : Kmeans%Size(ic) + nint(inclusionFraction * extraSize) /= Kmeans%Prop%EffectiveSize(ic)"
                        !            write(*,*) "count(MahalSq(1:ipstart+nps-2,icstart) <= 1._RK) : ", count(MahalSq(1:ipstart+nps-2,icstart) <= 1._RK)
                        !            write(*,*) "count(MahalSq(1:ipend+nps,icstart) <= 1._RK) : ", count(MahalSq(1:ipend+nps,icstart) <= 1._RK)
                        !            write(*,*) "Kmeans%Size(ic) : ", Kmeans%Size(ic)
                        !            write(*,*) "extraSize : ", extraSize
                        !            write(*,*) "Kmeans%Size(ic) + extraSize : ", Kmeans%Size(ic) + extraSize
                        !            write(*,*) "Kmeans%Prop%EffectiveSize(ic) : ", Kmeans%Prop%EffectiveSize(ic)
                        !            write(*,*) "logVolRatio < 0._RK :", logVolRatio < 0._RK
                        !            error stop
                        !        end if
                        !    end block
                        !end if
#endif
                        do j = 1, nd
                            ChoDia(j,icstart) = Kmeans%Prop%ChoDia(j,ic) * scaleFactor
                            InvCovMat(1:nd,j,icstart) = Kmeans%Prop%InvCovMat(1:nd,j,ic) * scaleFactorSqInverse
                            do i = j + 1, nd
                                ChoLowCovUpp(i,j,icstart) = Kmeans%Prop%ChoLowCovUpp(i,j,ic) * scaleFactor
                            end do
                        end do

                    else blockExpansionCheck
#endif
                        ChoDia(1:nd,icstart) = Kmeans%Prop%ChoDia(1:nd,ic)
                        LogVolNormed(icstart) = Kmeans%Prop%LogVolNormed(ic)
                        InvCovMat(1:nd,1:nd,icstart) = Kmeans%Prop%InvCovMat(1:nd,1:nd,ic)
                        ChoLowCovUpp(1:nd,1:nd,icstart) = Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic)
#if defined MAXDEN
                        if (boundedRegionIsTooLarge) LogLikeFitness(icstart) = -LogLikeFitness(icstart)
                        EffectiveSize(icstart) = Kmeans%Prop%EffectiveSize(ic)
                        scaleFactorSqInverse4KPM = 1._RK / Kmeans%Prop%ScaleFactorSq(ic) ! assumes no MahalSq normalization in Kmeans%getProp()
                        MahalSq(1:np,icstart) = Kmeans%Prop%MahalSq(1:np,ic) * scaleFactorSqInverse4KPM
                        !MahalSq(1:nps-1,icstart) = KmeansMahalSq(1:nps-1,ic)
                        !MahalSq(nps:npe,icstart) = Kmeans%Prop%MahalSq(1:npc,ic) * scaleFactorSqInverse4KPM
                        !MahalSq(npe+1:np,icstart) = KmeansMahalSq(npe+1:np,ic)
#if (defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED)
                        !with the new approach to computing EffectiveSize, these tests should be irrelevant.
                        !if (any(MahalSq(ipstart+nps-1:ipend+nps-1,icstart) > 1._RK)) then
                        !    write(*,*) pack(MahalSq(ipstart+nps-1:ipend+nps-1,icstart), MahalSq(ipstart+nps-1:ipend+nps-1,icstart) > 1._RK)
                        !    write(*,*) PROCEDURE_NAME//": Internal error occurred : MahalSq(1:ipstart+nps-1,icstart) > 1._RK :"
                        !    write(*,*) "logVolRatio < 0._RK :", logVolRatio < 0._RK
                        !    error stop
                        !else
                        !    block
                        !        integer :: effectiveSize
                        !        !write(*,*) "1", inclusionFraction
                        !        !write(*,*) "2", count(MahalSq(ipstart+nps-1:ipend+nps-1,icstart) <= 1._RK)
                        !        !write(*,*) "3", count(MahalSq(nps:ipstart+nps-2,icstart) <= 1._RK) + count(MahalSq(ipend+nps:npe,icstart) <= 1._RK)
                        !        !write(*,*) "4", count(MahalSq(1:nps-1,icstart) <= 1._RK) + count(MahalSq(npe+1:np,icstart) <= 1._RK)
                        !        effectiveSize   = count(MahalSq(ipstart+nps-1:ipend+nps-1,icstart) <= 1._RK) & ! LCOV_EXCL_LINE ! Kmeans%Size(ic)
                        !                        + nint(inclusionFraction * & ! LCOV_EXCL_LINE
                        !                        ( count(MahalSq(nps:ipstart+nps-2,icstart) <= 1._RK) + count(MahalSq(ipend+nps:npe,icstart) <= 1._RK) & ! LCOV_EXCL_LINE
                        !                        + count(MahalSq(1:nps-1,icstart) <= 1._RK) + count(MahalSq(npe+1:np,icstart) <= 1._RK) & ! LCOV_EXCL_LINE
                        !                        ) )
                        !        if (effectiveSize /= Kmeans%Prop%EffectiveSize(ic)) then
                        !            write(*,*) PROCEDURE_NAME//": Internal error occurred : effectiveSize /= Kmeans%Prop%EffectiveSize(ic) :"
                        !            write(*,*) "Kmeans%Prop%EffectiveSize(ic) :", Kmeans%Prop%EffectiveSize(ic)
                        !            write(*,*) "effectiveSize :", effectiveSize
                        !            write(*,*) "count(MahalSq(ipstart+nps-1:ipend+nps-1,icstart) <= 1._RK) == Kmeans%Size(ic) :", count(MahalSq(ipstart+nps-1:ipend+nps-1,icstart) <= 1._RK), Kmeans%Size(ic)
                        !            !write(*,*) "EffectiveSizeOld(ic) == effSizeOld :", EffectiveSizeOld(ic), effSizeOld
                        !            error stop
                        !        end if
                        !    end block
                        !end if
#endif
                    end if blockExpansionCheck
#endif
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
                        real(RK) :: mahalSqScalar
                        real(RK), allocatable :: NormedPoint(:)
                        do ip = ipstart, ipend
                            NormedPoint = Point(:,ip) - Center(:,icstart)
                            mahalSqScalar = dot_product(NormedPoint,matmul(InvCovMat(:,:,icstart),NormedPoint))
                            isInside = mahalSqScalar - 1._RK <= 1.e-6_RK
                            if (.not. isInside) then
                                ! LCOV_EXCL_START
                                write(*,"(*(g0.15,:,' '))") new_line("a"), PROCEDURE_NAME//": FATAL - POINT NOT INSIDE!, MAHAL = ", sqrt(mahalSqScalar), new_line("a")
                                write(*,"(60(g0,:,','))") "numRecursiveCall", numRecursiveCall
                                write(*,"(60(g0,:,','))") "boundedRegionIsTooLarge", boundedRegionIsTooLarge
                                write(*,"(60(g0,:,','))") "subclusterSearchWarranted", subclusterSearchWarranted
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

#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
            if (neopt < 1) then
                write(*,*) "Internal error occurred: neopt <= 1 :: neopt = ", neopt
                error stop
            end if
#endif

#if defined MINVOL
            if (parEllLogVolNormed < sum(LogVolNormed(1:neopt)) + logShrinkage) then ! recovery
#elif defined MAXDEN
            if (logShrinkage > NEGINF_RK) then
                logLikeSum = 0._RK
                loopSumLogLikeComputation: do ic = 1, neopt
                    if (LogLikeFitness(ic) > 0._RK) then
                        logLikeSum = logLikeSum - LogLikeFitness(ic)
                        !logLikeSum = POSINF_RK
                        !exit loopSumLogLikeComputation
                    else
                        logLikeSum = logLikeSum + LogLikeFitness(ic)
                    end if
                end do loopSumLogLikeComputation
                !if (logLikeSum < POSINF_RK) then
                    !call random_number(unifrnd);
                    !if (unifrnd > 0._RK) then; unifrnd = log(unifrnd); else; unifrnd = -huge(unifrnd); end if ! LCOV_EXCL_LINE or unifrnd = log(max(tiny(unifrnd),unifrnd))
                    !boundedRegionIsTooLarge = logLikeSum - parEllLogLike - logShrinkage < unifrnd ! dummy
                    !if (boundedRegionIsTooLarge) LogLikeFitness(1) = -parEllLogLike
                !else ! ill-defined
                !    boundedRegionIsTooLarge = .true.
                !    LogLikeFitness(1) = parEllLogLike
                !end if
                !if (boundedRegionIsTooLarge) then ! recovery
                if (logLikeSum + logShrinkage < parEllLogLike) then ! recovery
write(*,*) "sub2par likelihood: ", exp(logLikeSum - parEllLogLike - logShrinkage)
                    LogLikeFitness(1) = parEllLogLike ! -parEllLogLike
#endif
                    Center(1:nd,1) = ParEllCenter(1:nd)
                    ChoDia(1:nd,1) = ParEllChoDia(1:nd)
                    LogVolNormed(1) = parEllLogVolNormed
                    InvCovMat(1:nd,1:nd,1) = ParEllInvCovMat(1:nd,1:nd)
                    ChoLowCovUpp(1:nd,1:nd,1) = ParEllChoLowCovUpp(1:nd,1:nd)
                    PartitionSize(1) = npc
                    neopt = 1
                    return
                end if
#if defined MAXDEN
            end if
#endif

#if defined MINVOL
        else blockSubclusterSearch  ! one cluster is better

            PartitionSize(1) = npc
            neopt = 1
            return

        end if blockSubclusterSearch
#endif

#if defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED
        block
            integer :: ic, ip
            logical :: isInside
            real(RK) :: mahalSqScalar
            real(RK), allocatable :: NormedPoint(:)
            do ip = nps, npe
                ic = Membership(ip)
                if (ic < 1 .or. ic > nemax) then
                    write(*,*) PROCEDURE_NAME
                    write(*,*) "ic, nc, ip, npc, nemax, numRecursiveCall", ic, nc, ip, npc, nemax, numRecursiveCall
                    write(*,*) "Kmeans%Prop%CumSumSize", Kmeans%Prop%CumSumSize
                    error stop
                end if
                NormedPoint = Point(:,ip) - Center(:,ic)
                mahalSqScalar = dot_product(NormedPoint,matmul(InvCovMat(:,:,ic),NormedPoint))
                isInside = mahalSqScalar - 1._RK <= 1.e-6_RK
                if (.not. isInside) then
                    ! LCOV_EXCL_START
                    write(*,"(*(g0.15,:,' '))") new_line("a"), PROCEDURE_NAME//": FATAL - POINT NOT INSIDE!, MAHALSQ = ", mahalSqScalar, new_line("a")
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

#if defined DEBUG_ENABLED || defined TESTING_ENABLED || defined CODECOVE_ENABLED
    function getLogProbPoisDiff(count, logVolRatio) result(logProbPoisDiff)
#else
    pure function getLogProbPoisDiff(count, logVolRatio) result(logProbPoisDiff)
#endif
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbPoisDiff
#endif
        use Constants_mod, only: IK, RK, NEGINF_RK
        implicit none
        integer , intent(in)    :: count
        real(RK)    , intent(in)    :: logVolRatio
        real(RK)                    :: logProbPoisDiff
        logProbPoisDiff = count * (1._RK + logVolRatio - exp(logVolRatio))
        !if (probPoisDiff < NEGINF_RK) then ! @todo: xxx NEGINF_RK could be replaced with log(epsilon) for efficiency reasons.
        !    probPoisDiff = 0._RK
        !else
        !    probPoisDiff = exp(probPoisDiff)
        !end if
#if defined DEBUG_ENABLED || defined TESTING_ENABLED || defined CODECOVE_ENABLED
        if (logProbPoisDiff < NEGINF_RK .or. logProbPoisDiff > 0._RK) then
            write(*,*) "Internal error occurred: logProbPoisDiff < NEGINF_RK .or. logProbPoisDiff > 0._RK"
            write(*,*) "count, logVolRatio, logProbPoisDiff:", count, logVolRatio, logProbPoisDiff
            error stop
        end if
#endif
    end function getLogProbPoisDiff

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !function optimizePartition  ( nd & ! LCOV_EXCL_LINE
    !                            , np & ! LCOV_EXCL_LINE
    !                            , nps & ! LCOV_EXCL_LINE
    !                            , npe & ! LCOV_EXCL_LINE
    !                            , Point & ! LCOV_EXCL_LINE
    !                            , psize & ! LCOV_EXCL_LINE
    !                            , Center & ! LCOV_EXCL_LINE
    !                            , ChoDia & ! LCOV_EXCL_LINE
    !                            , MahalSq & ! LCOV_EXCL_LINE
    !                            , InvCovMat & ! LCOV_EXCL_LINE
    !                            , ChoLowCovUpp & ! LCOV_EXCL_LINE
    !                            , logVolNormed & ! LCOV_EXCL_LINE
    !                            , effectiveSize & ! LCOV_EXCL_LINE
    !                            , logLikeFitness & ! LCOV_EXCL_LINE
    !                            )
    !    implicit none
    !    integer , intent(in)    :: nd, np
    !end function optimizePartition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine writePartition   ( Partition & ! LCOV_EXCL_LINE
                                , fileUnit & ! LCOV_EXCL_LINE
                                , Point & ! LCOV_EXCL_LINE
                                )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: writePartition
#endif
        use Constants_mod, only: IK, RK
        implicit none
        class(Partition_type)   , intent(in)    :: Partition
        integer                 , intent(in)    :: fileUnit
        real(RK)                , intent(in)    :: Point(Partition%nd,Partition%np)
        character(*)            , parameter     :: fileFormat = "(*(g0,:,','))"
        integer                                 :: i, j, ic

        write(fileUnit,"(*(g0,:,'"//new_line("a")//"'))") "nd", Partition%nd, "np", Partition%np, "nc", Partition%neopt

        write(fileUnit,"(A)") "numRecursiveCall"
        write(fileUnit,fileFormat) Partition%numRecursiveCall

        write(fileUnit,"(A)") "pointLogVolNormed"
        write(fileUnit,fileFormat) Partition%pointLogVolNormed

        write(fileUnit,"(A)") "parentLogVolNormed"
        write(fileUnit,fileFormat) Partition%parentLogVolNormed

        write(fileUnit,"(A)") "Size"
        write(fileUnit,fileFormat) Partition%Size(1:Partition%neopt)

        write(fileUnit,"(A)") "Center"
        write(fileUnit,fileFormat) Partition%Center(1:Partition%nd,1:Partition%neopt)

        write(fileUnit,"(A)") "LogVolume"
        write(fileUnit,fileFormat) Partition%LogVolNormed(1:Partition%neopt)

        write(fileUnit,"(A)") "CholeskyLower"
        write(fileUnit,fileFormat) ((Partition%ChoDia(j,ic), (Partition%ChoLowCovUpp(i,j,ic), i=j+1,Partition%nd), j=1,Partition%nd), ic=1,Partition%neopt)

        write(fileUnit,"(A)") "Point"
        write(fileUnit,fileFormat) Point
#if defined MAXDEN
        write(fileUnit,"(A)") "LogLikeFitness"
        write(fileUnit,fileFormat) Partition%LogLikeFitness(1:Partition%neopt)
#endif
        write(fileUnit,"(A)") "Membership"
        write(fileUnit,fileFormat) Partition%Membership(1:Partition%np)

    end subroutine writePartition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%