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

module PartitionOptDen_mod

    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type

    implicit none

    character(*), parameter :: MODULE_NAME = "@PartitionOptDen_mod"

    type :: ParLev_type
        integer                     :: nc, nemax
        integer                     :: nps, npe, np
        integer     , allocatable   :: EffectiveSize(:)
        integer     , allocatable   :: Membership   (:)
        integer     , allocatable   :: CumSumSize   (:)
        integer     , allocatable   :: Size         (:)
        real(RK)    , allocatable   :: Center       (:,:)
        real(RK)    , allocatable   :: ChoDia       (:,:)
        real(RK)    , allocatable   :: InvCovMat    (:,:,:)
        real(RK)    , allocatable   :: ChoLowCovUpp (:,:,:)
        real(RK)    , allocatable   :: ScaleFactorSq(:)
        real(RK)    , allocatable   :: LogVolNormed (:)
    end type ParLev_type

    type KmeansTry_type
        real(RK)                    :: potential
        real(RK)    , allocatable   :: Center       (:.:)
        integer     , allocatable   :: Membership   (:)
        integer     , allocatable   :: Size         (:)
        type(Err_type)              :: Err
    end type KmeansTry_type

    !> \brief
    !> The `Partition_type` class. Partitions an input array `Point(nd,np)` with `nd` attributes and `np` observations (points).
    type :: Partition_type
        integer                     :: nsim                         !< See the interface of [constructPartition()](@ref constructpartition).
        integer                     :: nd,np                        !< See the interface of [constructPartition()](@ref constructpartition).
        integer                     :: nc,nt                        !< See the interface of [constructPartition()](@ref constructpartition).
        integer                     :: nemax                        !< See the interface of [constructPartition()](@ref constructpartition).
        integer                     :: minSize                      !< See the interface of [constructPartition()](@ref constructpartition).
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
        integer     , allocatable   :: EffectiveSize(:)             !< An array of size `(nemax)` indicating the likelihood of subclustering (`0<<1` if successful, else `1<<2` if failed), enlargement (if < 0), being warranted. If it is warranted but fails, the likelihood will be negative.
        real(RK)    , allocatable   :: LogLikeFitness(:)            !< An array of size `(nemax)` indicating the likelihood of subclustering (`0<<1` if successful, else `1<<2` if failed), enlargement (if < 0), being warranted. If it is warranted but fails, the likelihood will be negative.
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
                                , parentLogVolNormed & ! LCOV_EXCL_LINE
                                , mahalSqWeightExponent & ! LCOV_EXCL_LINE
                                , maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                !, scaleOptimizationEnabled & ! LCOV_EXCL_LINE
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
       !logical     , intent(in)    , optional  :: scaleOptimizationEnabled
        type(Partition_type)                    :: Partition

        Partition%nd = size(Point(:,1))
        Partition%np = size(Point(1,:))

        Partition%nc = 2; if (present(nc)) Partition%nc = max(Partition%nc, nc)
        Partition%nt = 1; if (present(nt)) Partition%nt = max(Partition%nt, nt)
        Partition%nsim = 0; if (present(nsim)) Partition%nsim = max(Partition%nsim, nsim)
        Partition%minSize = Partition%nd + 1; if (present(minSize)) Partition%minSize = minSize
        Partition%nemax = Partition%np / max(1, Partition%minSize); if (present(nemax)) Partition%nemax = max(1, nemax)
        !Partition%scaleOptimizationEnabled = .true.; if (present(scaleOptimizationEnabled)) Partition%scaleOptimizationEnabled = scaleOptimizationEnabled
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
                , Partition%EffectiveSize   (Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%LogLikeFitness  (Partition%nemax) & ! LCOV_EXCL_LINE
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
                Partition%EffectiveSize     = Partition%EffectiveSize   (1:Partition%neopt)
                Partition%LogLikeFitness    = Partition%LogLikeFitness  (1:Partition%neopt)
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
        use Matrix_mod, only: getInvMatFromCholFac, getEye
        use Matrix_mod, only: getCholeskyFactor
        use Kmeans_mod, only: runKmeans
        use String_mod, only: num2str
        use Sort_mod, only: sortAscending !, indexArray
        use Statistics_mod, only: getMean

        implicit none

        class(Partition_type)   , intent(inout) :: Partition
        real(RK)                , intent(inout) :: Point(Partition%nd,Partition%np) ! This must appear after Partition declaration.
        real(RK)                , intent(in)    :: parentLogVolNormed

        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME//"@runPartition()"

       !real(RK)                    :: logLikeShrinkage
       !real(RK)                    :: SumPoint(Partition%nd)
       !real(RK)                    :: NewCenter(Partition%nd)
       !real(RK)                    :: NewChoDia(Partition%nd)
       !real(RK)                    :: SinglePartitionCenter(Partition%nd)
       !real(RK)                    :: NewInvCovMat(Partition%nd,Partition%nd)
       !real(RK)                    :: NewChoLowCovUpp(Partition%nd,Partition%nd)
       !real(RK)                    :: MahalSqSorted(Partition%np)
       !real(RK)    , allocatable   :: MahalSqAcrossPartition(:)
       !integer     , allocatable   :: MahalSqSortedIndex(:)
       !integer                     :: newEffeciveSize
       !integer                     :: ipstart, ipend
       !integer                     :: ipMaxMahal
       !integer                     :: icMinMahal
       !integer                     :: newSize
       !integer                     :: jp
       !real(RK)                    :: MahalSqSinglePartition(Partition%np)

        real(RK)                :: ndHalf
        real(RK)                :: unifrnd
        real(RK)                :: ndInverse
        real(RK)                :: potential
        real(RK)                :: logVolRatio
        real(RK)                :: scaleFactor
        real(RK)                :: scaleFactorSq
        real(RK)                :: ndHalfInverse
        real(RK)                :: maxLogVolNormed
        real(RK)                :: scaleFactorSqInverse
        real(RK)                :: PointAvg(Partition%nd)
        real(RK)                :: PointStd(Partition%nd)
        real(RK)                :: PointStdInv(Partition%nd)
        real(RK)                :: MinDistanceSq(Partition%np)
        real(RK)                :: expectedPartitionLogVolNormed
        real(RK)                :: MahalSq(Partition%np,Partition%nc)
        real(RK)                :: NormedPoint(Partition%nd,Partition%np)

        real(RK)                :: InvCovMatSinglePoint(Partition%nd,Partition%nd)
        real(RK)                :: ChoLowCovUppSinglePoint(Partition%nd,Partition%nd)
        real(RK)                :: scaleFactorSqInverseSinglePoint
        real(RK)                :: scaleFactorSinglePoint

        type(ParLev_type)       :: ParLev(0:Partition%nemax)
        type(KmeansTry_type)    :: KmeansTry(1:Partition%nt)
        logical                 :: boundedRegionIsTooLarge

        integer                 :: nps
       !integer                 :: nemax
        integer                 :: ndPlusOne
        integer                 :: minSizePartition
        integer                 :: knc(Partition%nemax)
        integer                 :: knps(Partition%nemax)
        integer                 :: knpe(Partition%nemax)
        integer                 :: knpc(Partition%nemax)
        integer                 :: PointIndex(Partition%np)
        integer                 :: KmeansMemberCounter(Partition%nc)
        integer                 :: knemax(Partition%nc,Partition%nemax)
        integer                 :: i,j,ip,jp,ic,ie,nd,is,it,itmin
        integer                 :: ipstart,ipend
        integer                 :: niter,nfail

        nd = Partition%nd
        ndPlusOne = nd + 1
        minSizePartition = max(1,Partition%minSize)

        ndHalf = 0.5_RK * nd
        ndInverse = 1._RK / nd
        ndHalfInverse = 1._RK / ndHalf

        Partition%Err%occurred = .false.
        Partition%numRecursiveCall = 0
        Partition%convergenceFailureCount = 0
        Partition%Membership(1:Partition%np) = 0
        Partition%parentLogVolNormed = parentLogVolNormed ! xxx Does this need to be a component of Partition anymore?
        Partition%pointLogVolNormed = Partition%parentLogVolNormed - log(real(Partition%np,RK))

        do concurrent(ip = 1:Partition%np)
            Partition%PointIndex(ip) = ip
        end do

        ! Set the singular partition properties. ! xxx this could be perhaps done only when needed.

        scaleFactorSinglePoint = exp(ndInverse * Partition%pointLogVolNormed)
        scaleFactorSq = scaleFactorSinglePoint**2
        scaleFactorSqInverseSinglePoint = 1._RK / scaleFactorSq
        ChoLowCovUppSinglePoint(1:nd,1:nd) = getEye(nd, nd, scaleFactorSq)
        InvCovMatSinglePoint(1:nd,1:nd) = getEye(nd, nd, scaleFactorSqInverseSinglePoint)

        if (Partition%nt > 1) then
            do it = 1, Partition%nt
                allocate( KmeansTry(it)%Size(Partition%nc) & ! LCOV_EXCL_LINE
                        , KmeansTry(it)%Center(nd,Partition%nc) & ! LCOV_EXCL_LINE
                        , KmeansTry(it)%Membership(Partition%np) & ! LCOV_EXCL_LINE
                        )
            end do
        end if

        ie = 1
        nps = 1

        loopRecursivePartition: do

            ParLev(0)%nc = 1
            ParLev(0)%ic = 1
            ParLev(0)%nps = nps
            ParLev(0)%npe = Partition%np
            ParLev(0)%np = Partition%np - nps + 1
            ParLev(0)%nemax = Partition%nemax - ie + 1
            ParLev(0)%Size(1) = Partition%np - nps + 1

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Compute the properties of singular clusters.
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            blockSinglePointPartition: if (ParLev(0)%np < ndPlusOne) then

                ! Assume a hyper-sphere as the bounding region of each point.

                loopOverSinglePointPartition: do ip = nps, Partition%np

                    if (ie == Partition%nemax) then ! only one partition left.
                        Partition%neopt = ie
                        Partition%Size(ie) = Partition%np - ip + 1
                        Partition%Membership(ip) = ie
                        Partition%Center(1:nd,ie) = getMean(nd, Partition%Size(ie), Point(1:nd,ip:Partition%np))
                        scaleFactorSq = NEGINF_RK
                        do jp = ip, Partition%np
                            MahalSq(jp,1) = sum( ( Point(1:nd,jp) - Partition%Center(1:nd,ie) )**2 )
                            if (MahalSq(jp,1) > scaleFactorSq) scaleFactorSq = MahalSq(jp,1)
                        end do
                        Partition%LogVolNormed(ie) = max( ndHalf * log(scaleFactorSq) , Partition%pointLogVolNormed + log(real(Partition%Size(ie),RK)) )
                        scaleFactor = exp(ndInverse * Partition%LogVolNormed(ie))
                        scaleFactorSq = scaleFactor**2
                        scaleFactorSqInverse = 1._RK / scaleFactorSq
                        Partition%ChoDia(1:nd,ie) = scaleFactor
                        Partition%ChoLowCovUpp(1:nd,1:nd,ie) = getEye(nd, nd, scaleFactorSq)
                        Partition%InvCovMat(1:nd,1:nd,ie) = getEye(nd, nd, scaleFactorSqInverse)
                        Partition%LogLikeFitness(ie) = 0._RK
                        exit loopRecursivePartition
                    end if

                    Partition%Size(ie) = 1
                    Partition%Membership(ip) = ie
                    Partition%Center(1:nd,ie) = Point(1:nd,ip)
                    Partition%ChoDia(1:nd,ie) = scaleFactorSinglePoint
                    Partition%LogVolNormed(ie) = Partition%pointLogVolNormed
                    Partition%InvCovMat(1:nd,1:nd,ie) = InvCovMatSinglePoint
                    Partition%ChoLowCovUpp(1:nd,1:nd,ie) = ChoLowCovUppSinglePoint
                    Partition%LogLikeFitness(ie) = 0._RK

                    ie = ie + 1

                end do loopOverSinglePointPartition

                Partition%neopt = ie
                exit loopRecursivePartition

            end if blockSinglePointPartition

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Initialize the center of the cluster
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ParLev(0)%Center(1:nd,1) = 0._RK
            do ip = nps, Partition%np
              ParLev(0)%Center(1:nd,1) = ParLev(0)%Center(1:nd,1) + Point(1:nd,ip)
            end do
            ParLev(0)%Center(1:nd,1) = ParLev(0)%Center(1:nd,1) / ParLev(0)%Size(1)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Compute the covariance matrix.
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ParLev(0)%ChoLowCovUpp(1:nd,1:nd,1) = 0._RK
            do ip = nps, Partition%np
                do j = 1, nd
                    NormedPoint(j,ip) = Point(j,ip) - ParLev(0)%Center(j,1)
                    ParLev(0)%ChoLowCovUpp(1:j,j,1) = ParLev(0)%ChoLowCovUpp(1:j,j,1) + NormedPoint(1:j,ip) * NormedPoint(j,ip)
                end do
            end do
            ParLev(0)%ChoLowCovUpp(1:nd,1:nd,1) = ParLev(0)%ChoLowCovUpp(1:nd,1:nd,1) / (ParLev(0)%Size(1) - 1)

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Compute the Cholesky lower factor.
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            call getCholeskyFactor(nd, ParLev(0)%ChoLowCovUpp(1:nd,1:nd,1), ParLev(0)%ChoDia(1:nd,1))
            if (ParLev(0)%ChoDia(1,1)<0._RK) then
                ! LCOV_EXCL_START
                Partition%Err%msg = PROCEDURE_NAME//": Singular covariance matrix detected for ParLev(0)nemax, nd, Partition%Size(1): "// & ! LCOV_EXCL_LINE
                num2str(ParLev(0)%nemax)//", "//num2str(nd)//", "//num2str(ParLev(0)%Size(1))
                write(*,*) Partition%Err%msg ! xxx This needs to be fixed.
                error stop
                ! LCOV_EXCL_STOP
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Compute the inverse covariance matrix.
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ParLev(0)%InvCovMat(1:nd,1:nd,1) = getInvMatFromCholFac( nd & ! LCOV_EXCL_LINE
                                                                    , ParLev(0)%ChoLowCovUpp(1:nd,1:nd,1) & ! LCOV_EXCL_LINE
                                                                    , ParLev(0)%ChoDia(1:nd,1) & ! LCOV_EXCL_LINE
                                                                    )


            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Compute the MahalSq distances and the scale factor which makes the covariance matrix a bounding ellipsoid.
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            scaleFactorSq = NEGINF_RK
            do ip = nps, Partition%np
                MahalSq(ip,1) = dot_product( NormedPoint(1:nd,ip) , matmul(ParLev(0)%InvCovMat(1:nd,1:nd,1), NormedPoint(1:nd,ip)))
                if (scaleFactorSq < MahalSq(ip,1)) scaleFactorSq = MahalSq(ip,1)
            end do

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Find the optimal effective size of the bounded region.
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if (Partition%inclusionFraction > 1.e-6_RK .and. nps > 1) then
                ParLev(0)%EffectiveSize(1) = ParLev(0)%Size(1) + nint(Partition%inclusionFraction * count(MahalSq(1:nps-1,1)<=scaleFactorSq))
            else
                ParLev(0)%EffectiveSize(1) = ParLev(0)%Size(1)
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Compute the LogVolume and occurrence likelihood of the cluster based on the input estimated volume.
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            expectedPartitionLogVolNormed = Partition%pointLogVolNormed * log(real(ParLev(0)%EffectiveSize(1),RK))
            ParLev(0)%LogVolNormed(1) = sum( log( ParLev(0)%ChoDia(1:nd,1) ) ) + ndHalf * log(scaleFactorSq)
            logVolRatio = ParLev(0)%LogVolNormed(1) - expectedPartitionLogVolNormed
            ParLev(0)%LogLikeFitness(1) = getLogProbPoisDiff(ParLev(0)%Size(1), logVolRatio)

            boundedRegionIsTooLarge = logVolRatio > 0._RK

            if (boundedRegionIsTooLarge) then
                call random_number(unifrnd)
                if (unifrnd > 0._RK) then; unifrnd = log(unifrnd); else; unifrnd = -huge(unifrnd); end if ! LCOV_EXCL_LINE
                boundedRegionIsTooLarge = ParLev(0)%LogLikeFitness(1) < unifrnd
            else ! enlarge, no question asked.
                ParLev(0)%LogLikeFitness(1) = 0._RK
                ParLev(0)%LogVolNormed(1) = expectedPartitionLogVolNormed
                scaleFactorSq = scaleFactorSq * exp(-logVolRatio * ndHalfInverse)
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Subcluster or finish.
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            blockSubclusterSearch: if (ie < Partition%nemax .and. boundedRegionIsTooLarge) then

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! Search for subclusters.
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                is = 1
                ParLev(is)%nps = nps
                ParLev(is)%npe = Partition%np
                ParLev(is)%nemax = ParLev(0)%nemax - ie + 1
                ParLev(is)%nc = min(Partition%nc,ParLev(is)%nemax)
                ParLev(is)%np = ParLev(is)%npe - ParLev(is)%nps + 1
                ParLev(is)%ic = 1

                if (Partition%stanEnabled) then
                    PointAvg(1:Partition%nd) = Partition%Center(1:Partition%nd,ie)
                    do concurrent(i = 1:nd)
                        PointStd(i) = sqrt(ParLev(0)%ChoLowCovUpp(i,i,ie))
                    end do
                end if

                loopSubclusterSearch: do

                    Partition%numRecursiveCall = Partition%numRecursiveCall + 1

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ! Perform Kmeans clustering.
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    allocate( ParLev(is)%Size(ParLev(is)%nc) & ! LCOV_EXCL_LINE
                            , ParLev(is)%Center(nd,ParLev(is)%nc)) & ! LCOV_EXCL_LINE
                            , ParLev(is)%Membership(ParLev(is)%nps:ParLev(is)%npe) & ! LCOV_EXCL_LINE
                            )

                    blockStandardization: if (Partition%stanEnabled) then

                        PointStdInv(1:nd) = 1._RK / PointStd(1:nd)

                        do concurrent(ip = ParLev(is)%nps:ParLev(is)%npe)
                            NormedPoint(1:nd,ip) = (Point(1:nd,ip) - PointAvg(1:nd)) * PointStdInv(1:nd)
                        end do

#define STAN_ENABLED
#include "PartitionOptDen_mod.kmeans.inc.f90"
#undef STAN_ENABLED

                        do concurrent(ic = 1:ParLev(is)%nc)
                            ParLev(is)%Center(1:nd,ic) = ParLev(is)%Center(1:nd,ic) * PointStd(1:nd) + PointAvg(1:nd)
                        end do

                    else blockStandardization

#include "PartitionOptDen_mod.kmeans.inc.f90"

                    end if blockStandardization

#if (defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED)
                    block
                        integer :: sizeDum
                        real(RK) :: CenterDum(nd)
                        do ic = 1, ParLev(is)%nc
                            sizeDum = 0
                            CenterDum = 0._RK
                            do ip = ParLev(is)%nps, ParLev(is)%npe
                                if (ParLev(is)%Membership(ip) == ic) then
                                    CenterDum(1:nd) = CenterDum(1:nd) + Point(1:nd,ip)
                                    sizeDum = sizeDum + 1
                                end if
                            end do
                            CenterDum(1:nd) = CenterDum(1:nd) / sizeDum
                            if (any(abs(CenterDum(1:nd) - ParLev(is)%Center(1:nd,ic)) > 1.e-6_RK)) then
                                if (sizeDum /= ParLev(is)%Size(ic)) then; write(*,*) "sizeDum, Kmeans%Size = ", sizeDum, ParLev(is)%Size; error stop; endif
                                write(*,*) "nd,minSize,ic,nc,is     = ", nd, Partition%minSize, ic, ParLev(is)%nc, is
                                write(*,*) "ParLev(is)%Center(:)    = ", ParLev(is)%Center(1:nd,ic)
                                write(*,*) "CenterDum(:)            = ", CenterDum(1:nd)
                                write(*,*) "ParLev(is)%Size(:)      = ", ParLev(is)%Size
                                write(*,*) "ParLev(is)%Err%stat     = ", ParLev(is)%Err%stat
                                write(*,*) "ParLev(is)%Err%occurred = ", ParLev(is)%Err%occurred
                            endif
                        end do
                    end block
#endif

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ! Allocate and compute the partition properties.
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    allocate(ParLev(is)%ChoDia          (nd,ParLev(is)%nc))
                    allocate(ParLev(is)%InvCovMat       (nd,nd,ParLev(is)%nc))
                    allocate(ParLev(is)%ChoLowCovUpp    (nd,nd,ParLev(is)%nc))
                    allocate(ParLev(is)%EffectiveSize   (1:ParLev(is)%nc))
                    allocate(ParLev(is)%ScaleFactorSq   (1:ParLev(is)%nc))
                    allocate(ParLev(is)%LogVolNormed    (1:ParLev(is)%nc))
                    allocate(ParLev(is)%CumSumSize      (0:ParLev(is)%nc))

#include "PartitionOptDen_mod.prop.inc.f90"




                end do loopSubclusterSearch

            end if blockSubclusterSearch

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Rescale the Cholesky lower and diagonals to describe the bounding ellipsoid of the single partition.
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            Partition%neopt = ie
            scaleFactor = sqrt(scaleFactorSq)
            scaleFactorSqInverse = 1._RK / scaleFactorSq
            Partition%Membership(nps:Partition%np) = ie
            do j = 1, Partition%nd
                Partition%ChoDia(j,1) = Partition%ChoDia(j,1) * scaleFactor
                !Partition%ChoLowCovUpp(1:j,j,1) = Partition%ChoLowCovUpp(1:j,j,1) * scaleFactorSq
                Partition%InvCovMat(1:Partition%nd,j,1) = Partition%InvCovMat(1:Partition%nd,j,1) * scaleFactorSqInverse
                do i = j+1, Partition%nd
                    Partition%ChoLowCovUpp(i,j,1) = Partition%ChoLowCovUpp(i,j,1) * scaleFactor
                end do
            end do
            !Partition%EffectiveSize(ne) = Partition%np ! this must be fixed.

            exit loopRecursivePartition

        end do loopRecursivePartition

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
        integer             , intent(in)    :: fileUnit
        real(RK)                , intent(in)    :: Point(Partition%nd,Partition%np)
        character(*)            , parameter     :: fileFormat = "(*(g0,:,','))"
        integer                             :: i, j, ic

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

        write(fileUnit,"(A)") "LogLikeFitness"
        write(fileUnit,fileFormat) Partition%LogLikeFitness(1:Partition%neopt)

        write(fileUnit,"(A)") "Membership"
        write(fileUnit,fileFormat) Partition%Membership(1:Partition%np)

    end subroutine writePartition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module PartitionOptDen_mod ! LCOV_EXCL_LINE