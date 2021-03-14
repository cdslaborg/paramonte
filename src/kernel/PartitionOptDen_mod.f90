!Point has to become intent(in).
!remapping must be done for this.
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

    type :: ParTry_type
        logical                     :: deallocated = .true.
        integer                     :: ic
        integer                     :: nc
        integer                     :: np
        integer                     :: nps
        integer                     :: npe
        integer                     :: itp  ! parent partition level
        integer                     :: its  ! subpartition level
        integer     , allocatable   :: EffectiveSize(:)
       !integer     , allocatable   :: Membership(:)
        integer     , allocatable   :: CumSumSize(:)
        integer     , allocatable   :: Size(:)
        real(RK)    , allocatable   :: Center(:,:)
        real(RK)    , allocatable   :: ChoDia(:,:)
        real(RK)    , allocatable   :: MahalSq(:,:)
        real(RK)    , allocatable   :: InvCovMat(:,:,:)
        real(RK)    , allocatable   :: ChoLowCovUpp(:,:,:)
        real(RK)    , allocatable   :: LogVolNormed(:)
        real(RK)    , allocatable   :: ScaleFactorSq(:)
        real(RK)    , allocatable   :: LogLikeFitness(:)
        real(RK)    , allocatable   :: LogVolRatio(:)
        logical     , allocatable   :: IsTooLarge(:)
    end type ParTry_type

    type KmeansTry_type
        real(RK)                    :: potential
        real(RK)    , allocatable   :: Center(:.:)
        integer     , allocatable   :: Membership(:)
        integer     , allocatable   :: Size(:)
        type(Err_type)              :: Err
    end type KmeansTry_type

    !> \brief
    !> The `Partition_type` class. Partitions an input array `Point(nd,np)` with `nd` attributes and `np` observations (points).
    type :: Partition_type
        integer                         :: nsim                         !< See the interface of [constructPartition()](@ref constructpartition).
        integer                         :: ntry                         !< The number of partition tries (to be determined after partitioning).
        integer                         :: nd,np                        !< See the interface of [constructPartition()](@ref constructpartition).
        integer                         :: nc,nt                        !< See the interface of [constructPartition()](@ref constructpartition).
        integer                         :: nemax                        !< See the interface of [constructPartition()](@ref constructpartition).
        integer                         :: minSize                      !< See the interface of [constructPartition()](@ref constructpartition).
        integer                         :: maxAllowedKmeansFailure      !< See the interface of [constructPartition()](@ref constructpartition).
        integer                         :: maxAllowedKmeansRecursion    !< See the interface of [constructPartition()](@ref constructpartition).
        integer                         :: maxAllowedKvolumeRecursion   !< See the interface of [constructPartition()](@ref constructpartition).
        real(RK)                        :: mahalSqWeightExponent        !< See the interface of [constructPartition()](@ref constructpartition).
        real(RK)                        :: parentLogVolNormed           !< See the interface of [constructPartition()](@ref constructpartition).
        real(RK)                        :: pointLogVolNormed            !< The estimated volume of a single point based on the input value of `parentLogVolNormed`.
        real(RK)                        :: inclusionFraction            !< See the interface of [constructPartition()](@ref constructpartition).
        real(RK)                        :: logExpansion                 !< See the interface of [constructPartition()](@ref constructpartition).
        real(RK)                        :: logShrinkage                 !< See the interface of [constructPartition()](@ref constructpartition).
        real(RK)                        :: kmeansRelTol                 !< See the interface of [constructPartition()](@ref constructpartition).
        real(RK)                        :: nplog                        !< The natural logarithm of `np`.
        integer                         :: neopt                        !< The predicted optimal number of clusters identified in the input data.
        integer                         :: numRecursiveCall             !< The total number of recursive calls to the partitioning algorithm.
        integer                         :: convergenceFailureCount      !< The number of times the partitioning algorithm has failed to converge.
        logical                         :: stanEnabled                  !< See the interface of [constructPartition()](@ref constructpartition).
        type(ParTry_type), allocatable  :: Try(:)                       !< An array of size `(0:nemax)` containing the partition tries.
        integer         , allocatable   :: Size(:)                      !< An array of size `(nemax)` representing the sizes of the corresponding bounding ellipsoids.
        integer         , allocatable   :: CumSumSize(:)                !< An array of size `(0:nemax)` containing the Cumulative sum of the array `Size()`.
        integer         , allocatable   :: Membership(:)                !< An array of size `(np)` representing the bounding-ellipsoid membership IDs of the corresponding data points.
        integer         , allocatable   :: PointIndex(:)                !< An array of size `(np)` of indices such that the input `Point(1:nd,Index(ip))` is saved at the output `Point(1:nd,ip)`.
        integer         , allocatable   :: EffectiveSize(:)             !< An array of size `(nemax)` indicating the likelihood of subclustering (`0<<1` if successful, else `1<<2` if failed), enlargement (if < 0), being warranted. If further partitioning is warranted but fails, the likelihood will be negative.
        real(RK)        , allocatable   :: LogLikeFitness(:)            !< An array of size `(nemax)` indicating the likelihood of subclustering (`0<<1` if successful, else `1<<2` if failed), enlargement (if < 0), being warranted. If further partitioning is warranted but fails, the likelihood will be negative.
        real(RK)        , allocatable   :: LogVolNormed(:)              !< An array of size `(nemax)` representing the log-volumes of the corresponding bounding ellipsoids.
        real(RK)        , allocatable   :: ChoLowCovUpp(:,:,:)          !< An array of size `(nd,nd,nemax)` representing the Cholesky lower triangle, diagonal, and covariance matrices of the bounding ellipsoids.
        real(RK)        , allocatable   :: InvCovMat(:,:,:)             !< An array of size `(nd,nd,nemax)` representing the full symmetric inverse covariance matrices of the bounding ellipsoids.
       !real(RK)        , allocatable   :: MahalSq(:,:)                 !< An array of size `(np,nemax)` containing the MahalSq distances of all points from all partition centers.
        real(RK)        , allocatable   :: ChoDia(:,:)                  !< An array of size `(nd,nemax)` representing the Cholesky lower triangle, diagonal, and covariance matrices of the bounding ellipsoids.
        real(RK)        , allocatable   :: Center(:,:)                  !< An array of size `(nd,nemax)` representing the centers of the bounding ellipsoids.
        type(Err_type)                  :: Err                          !< An object of class [Err_type](@ref err_mod::err_type) containing information about error occurrence.
    contains
        procedure, pass                 :: write => writePartition
        procedure, pass                 :: run => runPartition
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
        Partition%nplog = log(real(Partition%np,RK))

        if (Partition%mahalSqWeightExponent == 0._RK .and. .not. present(parentLogVolNormed)) then ! xxx could we instead approximate `parentLogVolNormed`?
            Partition%Err%msg = "If `mahalSqWeightExponent` is given as input argument, then `parentLogVolNormed` must be also given as input."
            Partition%Err%occurred = .true.
            error stop
        end if

        allocate( Partition%Try             (0:Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%Center          (Partition%nd,Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%ChoDia          (Partition%nd,Partition%nemax) & ! LCOV_EXCL_LINE
               !, Partition%MahalSq         (Partition%np,Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%InvCovMat       (Partition%nd,Partition%nd,Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%ChoLowCovUpp    (Partition%nd,Partition%nd,Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%LogVolNormed    (1:Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%LogLikeFitness  (1:Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%EffectiveSize   (1:Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%CumSumSize      (0:Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%Size            (Partition%nemax) & ! LCOV_EXCL_LINE
                , Partition%Membership      (Partition%np) & ! LCOV_EXCL_LINE
                , Partition%PointIndex      (Partition%np) & ! LCOV_EXCL_LINE
                )

        Partition%CumSumSize(0) = 0 ! this must never change.
        call Partition%run(Point, parentLogVolNormed)

        if (present(trimEnabled)) then
            if (trimEnabled) then
                Partition%Try               = Partition%Try             (0:Partition%ntry)
                Partition%Center            = Partition%Center          (1:Partition%nd,1:Partition%neopt)
                Partition%ChoDia            = Partition%ChoDia          (1:Partition%nd,1:Partition%neopt)
                Partition%InvCovMat         = Partition%InvCovMat       (1:Partition%nd,1:Partition%nd,1:Partition%neopt)
                Partition%ChoLowCovUpp      = Partition%ChoLowCovUpp    (1:Partition%nd,1:Partition%nd,1:Partition%neopt)
                Partition%LogVolNormed      = Partition%LogVolNormed    (1:Partition%neopt)
                Partition%LogLikeFitness    = Partition%LogLikeFitness  (1:Partition%neopt)
                Partition%EffectiveSize     = Partition%EffectiveSize   (1:Partition%neopt)
                Partition%CumSumSize        = Partition%CumSumSize      (0:Partition%neopt)
                Partition%Size              = Partition%Size            (1:Partition%neopt)
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
        use Math_mod, only: getLogSumExp
        use Sort_mod, only: sortAscending !, indexArray
        use Statistics_mod, only: getMean

        implicit none

        class(Partition_type)   , intent(inout)                         :: Partition
        real(RK)                , intent(inout) , contiguous            :: Point(:,:) ! Partition%nd, Partition%np - must be after Partition declaration.
        real(RK)                , intent(in)                , optional  :: parentLogVolNormed

        character(*), parameter     :: PROCEDURE_NAME = MODULE_NAME//"@runPartition()"

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
       !real(RK)                :: PointAvg(Partition%nd)
        real(RK)                :: PointStd(Partition%nd)
        real(RK)                :: PointStdInv(Partition%nd)
        real(RK)                :: MinDistanceSq(Partition%np)
        real(RK)                :: expectedPartitionLogVolNormed
        real(RK)                :: MahalSq(Partition%np,Partition%nc)
        real(RK)                :: NormedPoint(Partition%nd,Partition%np)

        type(KmeansTry_type)    :: KmeansTry(1:Partition%nt)
        logical                 :: kvolumeEnabled

        integer                 :: itp
        integer                 :: nemax
        integer                 :: itmax
        integer                 :: ndPlusOne
        integer                 :: minSizePartition
        integer                 :: knc(Partition%nemax)
        integer                 :: knps(Partition%nemax)
        integer                 :: knpe(Partition%nemax)
        integer                 :: knpc(Partition%nemax)
        integer                 :: Membership(Partition%np)
        integer                 :: PointIndex(Partition%np)
        integer                 :: KmeansMemberCounter(Partition%nc)
        integer                 :: knemax(Partition%nc,Partition%nemax)
        integer                 :: i,j,ip,jp,ic,nd,kt,it,itmin
        integer                 :: ipstart,ipend
        integer                 :: niter,nfail

        nd = Partition%nd
        ndPlusOne = nd + 1
        ndHalf = 0.5_RK * nd
        ndInverse = 1._RK / nd
        ndHalfInverse = 1._RK / ndHalf

        minSizePartition = max(1,Partition%minSize)
        kvolumeEnabled = Partition%maxAllowedKvolumeRecursion > 0

        Partition%Err%occurred = .false.
        Partition%numRecursiveCall = 0
        Partition%convergenceFailureCount = 0
        Partition%Membership(1:Partition%np) = 0
        if (present(parentLogVolNormed)) Partition%parentLogVolNormed = parentLogVolNormed
        Partition%pointLogVolNormed = Partition%parentLogVolNormed - Partition%nplog

        do concurrent(ip = 1:Partition%np)
            Partition%PointIndex(ip) = ip
        end do

        if (Partition%nt > 1) then
            do kt = 1, Partition%nt
                allocate( KmeansTry(kt)%Size(Partition%nc) & ! LCOV_EXCL_LINE
                        , KmeansTry(kt)%Center(nd,Partition%nc) & ! LCOV_EXCL_LINE
                        , KmeansTry(kt)%Membership(Partition%np) & ! LCOV_EXCL_LINE
                        )
            end do
        end if

        it = 0
        itmax = 1
        nemax = Partition%nemax
        Partition%Try(0)%nc = 1
        Partition%Try(0)%ic = 0
        Partition%neopt = 0

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        loopRecursivePartition: do

            nemax = nemax - Partition%neopt
            Partition%numRecursiveCall = Partition%numRecursiveCall + 1

            blockForwardStream: if (Partition%Try(it)%ic == 0) then

                if (it > 0) Partition%Try(it)%nc = min(Partition%nc, nemax)

#if (defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED)
                if (nemax < 0) then; write(*,*) "nemax < 0"; error stop; end if
#endif

                if (Partition%Try(it)%deallocated) then
                    allocate( Partition%Try(it)%Size(Partition%nc) & ! LCOV_EXCL_LINE
                            , Partition%Try(it)%CumSumSize(Partition%nc)) & ! LCOV_EXCL_LINE
                            , Partition%Try(it)%Center(nd,Partition%nc)) & ! LCOV_EXCL_LINE
                            , Partition%Try(it)%ChoDia(nd,Partition%nc)) & ! LCOV_EXCL_LINE
                            , Partition%Try(it)%LogVolNormed(Partition%nc)) & ! LCOV_EXCL_LINE
                            , Partition%Try(it)%ScaleFactorSq(Partition%nc)) & ! LCOV_EXCL_LINE
                            , Partition%Try(it)%LogLikeFitness(Partition%nc)) & ! LCOV_EXCL_LINE
                            , Partition%Try(it)%InvCovMat(nd,nd,Partition%nc)) & ! LCOV_EXCL_LINE
                            , Partition%Try(it)%ChoLowCovUpp(nd,nd,Partition%nc)) & ! LCOV_EXCL_LINE
                           !, Partition%Try(it)%Membership(1:Partition%np) & ! LCOV_EXCL_LINE
                            , Partition%Try(it)%EffectiveSize(Partition%nc)) & ! LCOV_EXCL_LINE
                            , Partition%Try(it)%MahalSq(1:Partition%np) & ! LCOV_EXCL_LINE
                            )
                    Partition%Try(it)%deallocated = .false.
                end if

                blockKmeans: if (Partition%Try(it)%nc > 1) then

                    Partition%Try(it)%nps = Partition%Try(Partition%Try(it)%itp)%CumSumSize(Partition%Try(Partition%Try(it)%itp)%ic-1) + 1
                    Partition%Try(it)%npe = Partition%Try(Partition%Try(it)%itp)%CumSumSize(Partition%Try(Partition%Try(it)%itp)%ic)
                    Partition%Try(it)%np = Partition%Try(it)%npe - Partition%Try(it)%nps + 1

                    if (Partition%stanEnabled) then
                        do concurrent(i = 1:nd)
                            PointStd(i) = sqrt(Partition%Try(Partition%Try(it)%itp)%ChoLowCovUpp(i,i,Partition%Try(Partition%Try(it)%itp)%ic))
                        end do
                    end if

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ! Perform Kmeans clustering.
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    blockStandardization: if (Partition%stanEnabled) then
                        PointStdInv(1:nd) = 1._RK / PointStd(1:nd)
                        do concurrent(ip = Partition%Try(it)%nps:Partition%Try(it)%npe)
                            NormedPoint(1:nd,ip) = (Point(1:nd,ip) - Partition%Try(Partition%Try(it)%itp)%Center(1:nd)) * PointStdInv(1:nd)
                        end do
#define STAN_ENABLED
#include "PartitionOptDen_mod.kmeans.inc.f90"
#undef STAN_ENABLED
                        do concurrent(ic = 1:Partition%Try(it)%nc)
                            Partition%Try(it)%Center(1:nd,ic) = Partition%Try(it)%Center(1:nd,ic) * PointStd(1:nd) + Partition%Try(Partition%Try(it)%itp)%Center(1:nd)
                        end do
                    else blockStandardization
#include "PartitionOptDen_mod.kmeans.inc.f90"
                    end if blockStandardization

                else blockKmeans

                    if (it /= 0) error stop ! xxx this must be fixed.
                    Partition%Try(it)%npe = Partition%np
                    Partition%Try(it)%nps = Partition%CumSumSize(Partition%neopt) + 1
                    Partition%Try(it)%np = Partition%Try(it)%npe - Partition%Try(it)%nps + 1
                    Partition%Try(it)%Size(1) = Partition%Try(it)%np
                    !Partition%Try(it)%Membership = Partition%neopt + 1 ! xxx is this necessary?

                end if blockKmeans

#if (defined DEBUG_ENABLED || TESTING_ENABLED || CODECOVE_ENABLED)
                block
                    integer :: sizeDum
                    real(RK) :: CenterDum(nd)
                    do ic = 1, Partition%Try(it)%nc
                        sizeDum = 0
                        CenterDum = 0._RK
                        do ip = Partition%Try(it)%nps, Partition%Try(it)%npe
                            if (Partition%Try(it)%nc == 1 .or. Partition%Try(it)%Membership(ip) == ic) then
                                CenterDum(1:nd) = CenterDum(1:nd) + Point(1:nd,ip)
                                sizeDum = sizeDum + 1
                            end if
                        end do
                        CenterDum(1:nd) = CenterDum(1:nd) / sizeDum
                        if (any(abs(CenterDum(1:nd) - Partition%Try(it)%Center(1:nd,ic)) > 1.e-6_RK)) then
                            if (sizeDum /= Partition%Try(it)%Size(ic)) then; write(*,*) "sizeDum, Kmeans%Size = ", sizeDum, Partition%Try(it)%Size; error stop; endif
                            write(*,*) "nd,minSize,ic,nc,it     = ", nd, Partition%minSize, ic, Partition%Try(it)%nc, it
                            write(*,*) "Partition%Try(it)%Center(:)    = ", Partition%Try(it)%Center(1:nd,ic)
                            write(*,*) "CenterDum(:)            = ", CenterDum(1:nd)
                            write(*,*) "Partition%Try(it)%Size(:)      = ", Partition%Try(it)%Size
                            write(*,*) "Partition%Try(it)%Err%stat     = ", Partition%Try(it)%Err%stat
                            write(*,*) "Partition%Try(it)%Err%occurred = ", Partition%Try(it)%Err%occurred
                        endif
                    end do
                end block
#endif

#include "PartitionOptDen_mod.prop.inc.f90"

            end if blockForwardStream

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Wrap up the current search
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            blockCurrentSubclusterSearch: if (Partition%Try(it)%ic > 0) then
#include "PartitionOptDen_mod.terminal.inc.f90"
                if (it == 0) exit loopRecursivePartition
                it = Partition%Try(it)%itp
                cycle loopRecursivePartition
            end if blockCurrentSubclusterSearch

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Compute the likelihood and initiate the next search.
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            blockNextSubClusterSearch: if (Partition%Try(it)%ic < Partition%Try(it)%nc) then

                Partition%Try(it)%ic = Partition%Try(it)%ic + 1

                blockComputeLikelihood: if (Partition%Try(it)%Size(Partition%Try(it)%ic) > 0) then

                    ! Compute the effective fraction of points inside the bounded region.

                    if (Partition%inclusionFraction > 1.e-6_RK) then
                        Partition%Try(it)%EffectiveSize(Partition%Try(it)%ic) = count(Partition%Try(it)%MahalSq(1:Partition%np,Partition%Try(it)%ic) <= Partition%Try(it)%ScaleFactorSq(Partition%Try(it)%ic))
                    else
                        Partition%Try(it)%EffectiveSize(Partition%Try(it)%ic) = Partition%Try(it)%Size(Partition%Try(it)%ic)
                    end if

                    ! Compute the LogVolume and occurrence likelihood of the cluster based on the input estimated volume.

                    expectedPartitionLogVolNormed = Partition%pointLogVolNormed * log(real(Partition%Try(it)%EffectiveSize(Partition%Try(it)%ic),RK))
                    Partition%Try(it)%LogVolRatio(Partition%Try(it)%ic) = Partition%Try(it)%LogVolNormed(Partition%Try(it)%ic) - expectedPartitionLogVolNormed
                    Partition%Try(it)%LogLikeFitness(Partition%Try(it)%ic) = getLogProbPoisDiff(Partition%Try(it)%EffectiveSize(Partition%Try(it)%ic), Partition%Try(it)%LogVolRatio(Partition%Try(it)%ic))

                    Partition%Try(it)%IsTooLarge(Partition%Try(it)%ic) = Partition%Try(it)%LogVolRatio(Partition%Try(it)%ic) > 0._RK

                    if (Partition%Try(it)%IsTooLarge(Partition%Try(it)%ic)) then
                        call random_number(unifrnd)
                        if (unifrnd > 0._RK) then; unifrnd = log(unifrnd); else; unifrnd = -huge(unifrnd); end if ! LCOV_EXCL_LINE
                        Partition%Try(it)%IsTooLarge(Partition%Try(it)%ic) = Partition%Try(it)%LogLikeFitness(Partition%Try(it)%ic) < unifrnd
                    !else ! enlarge, no question asked.
                    !    Partition%Try(it)%LogLikeFitness(Partition%Try(it)%ic) = 0._RK
                    !    Partition%Try(it)%LogVolNormed(Partition%Try(it)%ic) = expectedPartitionLogVolNormed
                    !    Partition%Try(it)%ScaleFactorSq(Partition%Try(it)%ic) = Partition%Try(it)%ScaleFactorSq(Partition%Try(it)%ic) * exp(-Partition%Try(it)%LogVolRatio(Partition%Try(it)%ic) * ndHalfInverse)
                    end if

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    ! Search for subclusters or do not.
                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    if (Partition%Try(it)%IsTooLarge(Partition%Try(it)%ic)) then
                        itp = it
                        it = itmax
                        itmax = itmax + 1
                        Partition%Try(it)%ic = 0
                        Partition%Try(it)%itp = itp
                        Partition%Try(itp)%its = it
                        cycle loopRecursivePartition
                    else
#include "PartitionOptDen_mod.terminal.inc.f90"
                        if (it == 0) exit loopRecursivePartition
                        it = Partition%Try(it)%itp
                        cycle loopRecursivePartition
                    end if

                end if blockComputeLikelihood

            end if blockNextSubClusterSearch

        end do loopRecursivePartition

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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