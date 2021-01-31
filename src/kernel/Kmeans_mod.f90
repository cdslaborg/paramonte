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
!> This module contains procedures and routines for the Kmeans the statistical classification and clustering methods.
!> \author
!> Amir Shahmoradi

module Kmeans_mod

    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type

    implicit none

    character(*), parameter :: MODULE_NAME = "@Kmeans_mod"

    real(RK) :: ndKuniform, btFactor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> The `Kmeans_type` class.
    !> The inclusion of the component `NormedPoint` adds ~50% to the computational cost of performing Kmeans.
    !> In this latest implementation, avoiding `NormedPoint` as a component leads to 15-50% performance improvement.
    !> The moral is, simple repeating of addition operation is much significantly cheaper than storing the operation results for later usage.
    !> In the following declaration descriptions,
    !>     +    `nc` represents the number of clusters to be found.
    !>     +    `np` represents the number of data points to be clustered.
    !>     +    `nd` represents the number of dimensions (features) of the points to be clustered.
    type :: Kmeans_type
        real(RK)                    :: potential            !< The sum of the distances-squared of all points from their corresponding group centers.
        integer(IK)                 :: nfail                !< The number of Zero-Sized Clusters Iterations.
        integer(IK)                 :: niter                !< The number of iterations to reach convergence.
        real(RK)    , allocatable   :: LogVol(:)            !< An array of size `(nc)` each element of which represents the volume of the covariance matrix of the corresponding cluster.
        real(RK)    , allocatable   :: ChoDia(:,:)          !< An array of size `(nd,nc)` representing the diagonal elements of the Cholesky factorization of the covariance matrix.
        real(RK)    , allocatable   :: Center(:,:)          !< An array of size `(nd,nc)` containing the centers of the individual clusters identified.
        real(RK)    , allocatable   :: MahalSq(:,:)         !< An array of size `(np,nc)` each element of which represents the Mahalanobis distance squared of point `ip` from the cluster `ic`.
        real(RK)    , allocatable   :: InvCovMat(:,:,:)     !< An array of size `(nd,nd,nc)` containing the inverse of the covariance matrix of the corresponding cluster.
       !real(RK)    , allocatable   :: NormedPoint(:,:,:)   !< An array of size `(nd,np,nc)` containing `nc` copies of the input array `Point(nd,np)` normalized to the cluster centers.
        real(RK)    , allocatable   :: ChoLowCovUpp(:,:,:)  !< An array of size `(nd,nd,nc)` whose upper triangle and diagonal is the covariance matrix and the lower is the Cholesky Lower.
        real(RK)    , allocatable   :: MinDistanceSq(:)     !< An array of size `(np)` of the minimum distances of the points from the cluster centers.
        real(RK)    , allocatable   :: ScaleFactorSq(:)     !< An array of size `(nc)` representing the factors by which the cluster covariance matrices must be enlarged to enclose their corresponding members.
        integer(IK) , allocatable   :: EffectiveSize(:)     !< An array of size `(nc)` representing the factors by which the cluster covariance matrices must be enlarged to enclose their corresponding members.
        integer(IK) , allocatable   :: Membership(:)        !< A vector of size `np` each element of which represents the cluster ID to which the point belongs.
        integer(IK) , allocatable   :: CumSumSize(:)        !< A vector of size `0:nc` representing the cumulative sum of all cluster sizes from cluster 1 to the last.
        integer(IK) , allocatable   :: Size(:)              !< A vector of size `nc` (the number of clusters) containing the sizes of the individual clusters identified.
        type(Err_type)              :: Err                  !< An object of class [Err_type](@ref err_mod::err_type) containing information about error occurrence.
    contains
        procedure, pass :: getProp
    end type Kmeans_type

    interface Kmeans_type
        module procedure :: getKmeans, allocateKmeans
    end interface Kmeans_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function allocateKmeans(nd, np, nc) result(Kmeans)
        implicit none
        integer(IK), intent(in) :: nd, np, nc
        type(Kmeans_type)       :: Kmeans
        allocate( Kmeans%Size(nc) & ! LCOV_EXCL_LINE
                , Kmeans%Center(nd,nc) & ! LCOV_EXCL_LINE
               !, Kmeans%NormedPoint(nd,np,nc) & ! LCOV_EXCL_LINE
                , Kmeans%MinDistanceSq(np) & ! LCOV_EXCL_LINE
                , Kmeans%Membership(np) & ! LCOV_EXCL_LINE
                )
    end function allocateKmeans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Perform the Kmeans clustering for `nt` tries on the input data set represented by the array `Point(nd, np)`
    !> and return the clustering that yields the least potential.
    !>
    !> \param[in]       nd                  :   See the description of the [runKmeans()](@ref runkmeans).
    !> \param[in]       np                  :   See the description of the [runKmeans()](@ref runkmeans).
    !> \param[in]       nc                  :   See the description of the [runKmeans()](@ref runkmeans).
    !> \param[in]       nt                  :   The number of independent instances of Kmeans to run.
    !> \param[inout]    Point               :   See the description of the [runKmeans()](@ref runkmeans).
    !> \param[inout]    Index               :   See the description of the [runKmeans()](@ref runkmeans).
    !> \param[in]       InitCenter          :   See the description of the [runKmeans()](@ref runkmeans).
    !> \param[in]       niterMax            :   See the description of the [runKmeans()](@ref runkmeans).
    !> \param[in]       nfailMax            :   See the description of the [runKmeans()](@ref runkmeans).
    !> \param[in]       minSize             :   The minimum acceptable size for any Kmeans cluster.
    !> \param[in]       relTol              :   See the description of the [runKmeans()](@ref runkmeans).
    !> \param[in]       propEnabled         :   If `.true.`, a call to `getProp()` method will be made (**optional**, default = `.false.`).
    !> \param[in]       inclusionFraction   :   See the description of the [getProp()](@ref getprop).
    !>
    !> \return
    !> `Kmeans` :   An object of type [Kmeans_type](@ref Kmeans_type) containing the properties of the identified clusters.
    !>              The `Err%stat` component of the output `Kmeans` will be set to `1` if the number of iterations reaches the input
    !>              `niterMax` or its default value before convergence occurs.
    !>              The `Err%stat` component of the output `Kmeans` will be set to `2` if more than `nfailMax` clustering tries
    !>              result in clusters with zero members.
    !>              In both cases in the above, `Err%occurred = .true.` upon exiting the procedure.
    !>
    !> \warning
    !> The input value for `nt` must be at least 1, otherwise it will lead to a runtime error or at best, the output will be undefined.
    !>
    !> \warning
    !> The input value for `niterMax` must be larger than 0 since at least one iteration will be performed by the procedure.
    !>
    !> \warning
    !> This algorithm does not check for the consistency of the input values for `nc`, `nd`, and `np`.
    !> It is currently the user's responsibility to ensure that these input values are logical and sensible.
    !>
    !> \author
    ! Amir Shahmoradi, April 03, 2017, 2:16 PM, ICES, UTEXAS
    function getKmeans(nd, np, nc, nt, Point, Index, InitCenter, niterMax, nfailMax, minSize, relTol, propEnabled, inclusionFraction) result(Kmeans)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKmeans
#endif
        use Constants_mod, only: IK, RK, HUGE_RK

        implicit none

        integer(IK) , intent(in)                :: nd, np, nc, nt
        real(RK)    , intent(inout)             :: Point(nd,np)
        integer(IK) , intent(inout) , optional  :: Index(np)
        real(RK)    , intent(in)    , optional  :: InitCenter(nd,nc)
        real(RK)    , intent(in)    , optional  :: relTol
        logical     , intent(in)    , optional  :: propEnabled
        integer(IK) , intent(in)    , optional  :: nfailMax
        integer(IK) , intent(in)    , optional  :: niterMax
        integer(IK) , intent(in)    , optional  :: minSize
        real(RK)    , intent(in)    , optional  :: inclusionFraction

        type(Kmeans_type), allocatable          :: KmeansTry(:)
        type(Kmeans_type)                       :: Kmeans
        real(RK)                                :: potential
       !real(RK)                                :: NormedPoint(nd,np)
        logical                                 :: propEnabledDefault
        integer(IK)                             :: itry, itryLeastPotential
        integer(IK)                             :: minSizeDefault

        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME//"@getKmeans()"

        minSizeDefault = 1_IK; if (present(minSize)) minSizeDefault = minSize

        if (present(propEnabled)) then
            propEnabledDefault = propEnabled
            if (propEnabledDefault) then
                minSizeDefault = max(minSizeDefault, nd + 1_IK)
            end if
        else
            propEnabledDefault = .false.
        end if

        if (allocated(KmeansTry)) deallocate(KmeansTry) ! LCOV_EXCL_LINE
        allocate(KmeansTry(nt))
        potential = HUGE_RK

        ! Compute the optimal Kmeans

        if (nt > 1_IK) then

            itryLeastPotential = 0_IK
            do itry = 1, nt
                KmeansTry(itry) = Kmeans_type(nd, np, nc)
                call runKmeans  ( nd = nd & ! LCOV_EXCL_LINE
                                , np = np & ! LCOV_EXCL_LINE
                                , nc = nc & ! LCOV_EXCL_LINE
                                , Point = Point & ! LCOV_EXCL_LINE
                                , Size = KmeansTry(itry)%Size & ! LCOV_EXCL_LINE
                                , Center = KmeansTry(itry)%Center & ! LCOV_EXCL_LINE
                                , Membership = KmeansTry(itry)%Membership & ! LCOV_EXCL_LINE
                               !, NormedPoint = KmeansTry(itry)%NormedPoint & ! LCOV_EXCL_LINE
                                , MinDistanceSq = KmeansTry(itry)%MinDistanceSq & ! LCOV_EXCL_LINE
                                , potential = KmeansTry(itry)%potential & ! LCOV_EXCL_LINE
                                , niter = KmeansTry(itry)%niter & ! LCOV_EXCL_LINE
                                , nfail = KmeansTry(itry)%nfail & ! LCOV_EXCL_LINE
                                , Err = KmeansTry(itry)%Err & ! LCOV_EXCL_LINE
                                , relTol = relTol & ! LCOV_EXCL_LINE
                                , InitCenter = InitCenter & ! LCOV_EXCL_LINE
                                , niterMax = niterMax & ! LCOV_EXCL_LINE
                                , nfailMax = nfailMax & ! LCOV_EXCL_LINE
                                )
                if ( KmeansTry(itry)%potential < potential .and. .not. any(KmeansTry(itry)%Size < minSizeDefault) ) then
                    potential = KmeansTry(itry)%potential
                    itryLeastPotential = itry
                end if
            end do

            if (itryLeastPotential > 0_IK) then
                Kmeans = KmeansTry(itryLeastPotential)
            else
                ! LCOV_EXCL_START
                Kmeans%Err%msg = "All Kmeans clustering tries failed."
                Kmeans%Err%occurred = .true.
                Kmeans%Err%stat = 1_IK
                return
                ! LCOV_EXCL_STOP
            end if

        else

            Kmeans = Kmeans_type(nd, np, nc)
            call runKmeans  ( nd = nd & ! LCOV_EXCL_LINE
                            , np = np & ! LCOV_EXCL_LINE
                            , nc = nc & ! LCOV_EXCL_LINE
                            , Point = Point & ! LCOV_EXCL_LINE
                            , Size = Kmeans%Size & ! LCOV_EXCL_LINE
                            , Center = Kmeans%Center & ! LCOV_EXCL_LINE
                            , Membership = Kmeans%Membership & ! LCOV_EXCL_LINE
                           !, NormedPoint = Kmeans%NormedPoint & ! LCOV_EXCL_LINE
                            , MinDistanceSq = Kmeans%MinDistanceSq & ! LCOV_EXCL_LINE
                            , potential = Kmeans%potential & ! LCOV_EXCL_LINE
                            , niter = Kmeans%niter & ! LCOV_EXCL_LINE
                            , nfail = Kmeans%nfail & ! LCOV_EXCL_LINE
                            , Err = Kmeans%Err & ! LCOV_EXCL_LINE
                            , relTol = relTol & ! LCOV_EXCL_LINE
                            , InitCenter = InitCenter & ! LCOV_EXCL_LINE
                            , niterMax = niterMax & ! LCOV_EXCL_LINE
                            , nfailMax = nfailMax & ! LCOV_EXCL_LINE
                            )

            if (any(Kmeans%Size < minSizeDefault)) then
                ! LCOV_EXCL_START
                Kmeans%Err%msg = "The Kmeans clustering sizes are smaller than the requested input minimum cluster size."
                Kmeans%Err%occurred = .true.
                Kmeans%Err%stat = 1_IK
                return
                ! LCOV_EXCL_STOP
            end if

        end if

        if (propEnabledDefault) call Kmeans%getProp(nd, np, Point, Index, inclusionFraction)

    end function getKmeans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> This is a method of class [Kmeans_type](@ref kmeans_type).
    !> Compute the following properties (components) of the input object of class [Kmeans_type](@ref kmeans_type):
    !>
    !> +   `ChoLowCovUpp(1:nd,1:nd,1:nc)`   :   The covariance matrix upper triangle and Cholesky factor lower triangle,
    !> +   `InvCovMat(1:nd,1:nd,1:nc)`      :   The full inverse covariance matrix,
    !> +   `EffectiveSize(1:nc)`            :   The effective sizes of the ellipsoidal bounded regions of the clusters.
    !> +   `ScaleFactorSq(1:nc)`            :   The factor by which the covariance matrix elements must be multiplied
    !>                                          to convert the matrix to a bounding ellipsoid of the cluster members.
    !> +   `ChoDia(1:nd,1:nc)`              :   The diagonal elements of the Cholesky lower triangle,
    !> +   `CumSumSize(0:nc)`               :   The cumulative sum of `Kmeans%Size(1:nc)` with `CumSumSize(0) = 0`,
    !> +   `MahalSq(1:np)`                  :   The Mahalanobis distances squared of all points from all cluster centers,
    !> +   `LogVol(1:nc)`                   :   The vector of natural-log volumes of all clusters,
    !>
    !> \param[in]       nd                  :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[in]       np                  :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[in]       nc                  :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[inout]    Point               :   The Array of size `(nd,np)` representing the original points used in Kmeans clustering.
    !> \param[inout]    Index               :   A vector of size `np` such that the output `Point(:,Index(ip),:)` points to the input `Point(:,Index(ip),:)` (**optional**).
    !> \param[in]       inclusionFraction   :   The fraction of points inside each cluster's bounding region to be considered in computing the effective sizes of 
    !>                                          the bounded regions of each cluster (**optional**, default = `0.`).
    !>
    !> \warning
    !> Upon entry to this method, the `Center`, `Size`, and `Membership` components of the Kmeans object
    !> must have been already properly set up by its class constructor or a call to `getKmeans()`.
    !>
    !> \warning
    !> On output, the arrays `Point` and `Index` will be rearranged on exit from the procedure such that the members of the same
    !> clusters appear next to each other contiguously in the array. Consequently, the output array component `Index(1:np)` will
    !> be populated such that the new `Point(1:nd,Index(ip))` points to the original input `Point(1:nd,ip)` or the original .
    !>
    !> \warning
    !> If the Cholesky factorization in the procedure fails, the control is returned with `Kmeans%Err%occurred = .true.`.
    pure subroutine getProp(Kmeans, nd, np, Point, Index, inclusionFraction)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProp
#endif
        use Constants_mod, only: IK, RK
        use Matrix_mod, only: getCholeskyFactor
        use Matrix_mod, only: getInvMatFromCholFac

        implicit none

        class(Kmeans_type) , intent(inout)      :: Kmeans
        integer(IK) , intent(in)                :: nd, np
        real(RK)    , intent(inout)             :: Point(nd,np)
        integer(IK) , intent(inout) , optional  :: Index(np)
        real(RK)    , intent(in)    , optional  :: inclusionFraction

        integer(IK) , allocatable               :: KmeansMemberCounter(:)
        integer(IK)                             :: KmeansIndex(np)
        integer(IK)                             :: j, ic, ip, nc
        integer(IK)                             :: ipstart, ipend
        real(RK)                                :: NormedPoint(nd,np)
        real(RK)                                :: ndHalf

        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME//"@getProp()"

        nc = size(Kmeans%Size)
        ndHalf = 0.5_RK * nd

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! If requested, compute the bounded region properties of the optimal Kmeans clusters.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (.not. allocated(Kmeans%ChoDia       )) allocate(Kmeans%ChoDia          (nd,nc))
        if (.not. allocated(Kmeans%MahalSq      )) allocate(Kmeans%MahalSq         (np,nc))
        if (.not. allocated(Kmeans%InvCovMat    )) allocate(Kmeans%InvCovMat       (nd,nd,nc))
        if (.not. allocated(Kmeans%ChoLowCovUpp )) allocate(Kmeans%ChoLowCovUpp    (nd,nd,nc))
        if (.not. allocated(Kmeans%EffectiveSize)) allocate(Kmeans%EffectiveSize   (nc))
        if (.not. allocated(Kmeans%ScaleFactorSq)) allocate(Kmeans%ScaleFactorSq   (nc))
        if (.not. allocated(Kmeans%CumSumSize   )) allocate(Kmeans%CumSumSize      (0:nc))
        if (.not. allocated(Kmeans%LogVol       )) allocate(Kmeans%LogVol          (nc))

        if (allocated(KmeansMemberCounter)) deallocate(KmeansMemberCounter) ! LCOV_EXCL_LINE
        allocate(KmeansMemberCounter(nc))

        ! Compute the CDF of the cluster Size array.

        Kmeans%CumSumSize(0) = 0_IK
        loopDetermineClusterBoundaries: do ic = 1, nc
            KmeansMemberCounter(ic) = Kmeans%CumSumSize(ic-1)
            Kmeans%CumSumSize(ic) = KmeansMemberCounter(ic) + Kmeans%Size(ic)
        end do loopDetermineClusterBoundaries

        ! Reorder Point based on the identified clusters.

        if (present(Index)) then
            do ip = 1, np
                KmeansMemberCounter(Kmeans%Membership(ip)) = KmeansMemberCounter(Kmeans%Membership(ip)) + 1_IK
                NormedPoint(1:nd,KmeansMemberCounter(Kmeans%Membership(ip))) = Point(1:nd,ip)
                KmeansIndex(KmeansMemberCounter(Kmeans%Membership(ip))) = Index(ip)
            end do
            Index = KmeansIndex
        else
            do ip = 1, np
                KmeansMemberCounter(Kmeans%Membership(ip)) = KmeansMemberCounter(Kmeans%Membership(ip)) + 1_IK
                NormedPoint(1:nd,KmeansMemberCounter(Kmeans%Membership(ip))) = Point(1:nd,ip)
                !KmeansIndex(KmeansMemberCounter(Kmeans%Membership(ip))) = ip
            end do
        end if
        Point = NormedPoint

        ! Compute the properties of the clusters.

        loopComputeClusterProperties: do ic = 1, nc

            ipstart = Kmeans%CumSumSize(ic-1) + 1_IK
            ipend = Kmeans%CumSumSize(ic)

            ! Correct the cluster memberships.

            Kmeans%Membership(ipstart:ipend) = ic

            ! Normalize points.

            do concurrent(ip = 1:np)
                NormedPoint(1:nd,ip) = Point(1:nd,ip) - Kmeans%Center(1:nd,ic)
            end do

            ! Compute the upper covariance matrix of the cluster covariance matrices.

            Kmeans%ChoLowCovUpp(1:nd,1:nd,ic) = 0._RK
            do ip = ipstart, ipend
                do j = 1, nd
                    Kmeans%ChoLowCovUpp(1:j,j,ic) = Kmeans%ChoLowCovUpp(1:j,j,ic) + NormedPoint(1:j,ip) * NormedPoint(j,ip)
                end do
            end do
            Kmeans%ChoLowCovUpp(1:nd,1:nd,ic) = Kmeans%ChoLowCovUpp(1:nd,1:nd,ic) / real(Kmeans%Size(ic)-1_IK, kind = RK)

            ! Compute the Cholesky Factor of the cluster covariance matrices.

            call getCholeskyFactor(nd, Kmeans%ChoLowCovUpp(1:nd,1:nd,ic), Kmeans%ChoDia(1:nd,ic))
            if (Kmeans%ChoDia(1,ic) < 0._RK) then
                ! LCOV_EXCL_START
                Kmeans%Err%msg = PROCEDURE_NAME//"Cholesky factorization failed."
                Kmeans%Err%occurred = .true.
                return
                ! LCOV_EXCL_STOP
            end if

            ! Compute the inverse of the cluster covariance matrices.

            Kmeans%InvCovMat(1:nd,1:nd,ic) = getInvMatFromCholFac( nd, Kmeans%ChoLowCovUpp(1:nd,1:nd,ic), Kmeans%ChoDia(1:nd,ic) )

            ! Compute the MahalSq of as many points as needed.

            do concurrent(ip = 1:np)
                Kmeans%MahalSq(ip,ic) = dot_product( NormedPoint(1:nd,ip) , matmul(Kmeans%InvCovMat(1:nd,1:nd,ic), NormedPoint(1:nd,ip)) )
            end do

            ! Compute the scaleFcator of the bounding region and scale the volumes to the bounded region.

            Kmeans%ScaleFactorSq(ic) = maxval(Kmeans%MahalSq(ipstart:ipend,ic))
            Kmeans%LogVol(ic) = sum( log(Kmeans%ChoDia(1:nd,ic)) ) + ndHalf * log(Kmeans%ScaleFactorSq(ic))
            !Kmeans%ScaleFactor(ic) = sqrt(Kmeans%ScaleFactorSq(ic))

        end do loopComputeClusterProperties

        ! Compute the effective fraction of points inside the bounded region.

        if (present(inclusionFraction)) then
            !ic = 1_IK; Kmeans%EffectiveSize(ic) = Kmeans%Size(ic) + nint(inclusionFraction*count(Kmeans%MahalSq(Kmeans%CumSumSize(ic)+1:np,ic)<Kmeans%ScaleFactorSq(ic)), kind=IK)
            !ic = nc  ; Kmeans%EffectiveSize(ic) = Kmeans%Size(ic) + nint(inclusionFraction*count(Kmeans%MahalSq(1:Kmeans%CumSumSize(ic-1),ic)<Kmeans%ScaleFactorSq(ic)), kind=IK)
            !do concurrent(ic = 2:nc-1)
            do concurrent(ic = 1:nc)
                Kmeans%EffectiveSize(ic) = Kmeans%Size(ic) & ! LCOV_EXCL_LINE
                                        + nint(inclusionFraction*count(Kmeans%MahalSq(1:Kmeans%CumSumSize(ic-1),ic)<Kmeans%ScaleFactorSq(ic)), kind=IK) & ! LCOV_EXCL_LINE
                                        + nint(inclusionFraction*count(Kmeans%MahalSq(Kmeans%CumSumSize(ic)+1:np,ic)<Kmeans%ScaleFactorSq(ic)), kind=IK)
                !write(*,*) ic, Kmeans%Size(ic), Kmeans%EffectiveSize(ic), abs(Kmeans%EffectiveSize(ic) - Kmeans%Size(ic))
            end do
        else
            Kmeans%EffectiveSize = Kmeans%Size
        end if

    end subroutine getProp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Perform the Kmeans clustering on the input data set represented by the array `Point(nd, np)`.
    !>
    !> \param[in]       nd              :   The number of dimensions (features) of the input data array `Point`.
    !> \param[in]       np              :   The number of observations in the input data array `Point`.
    !> \param[in]       nc              :   The number of clusters to be found via the Kmeans algorithm.
    !> \param[in]       Point           :   The input array of size `(nd,np)` of the points to be clustered via the Kmeans algorithm.
    !> \param[out]      Size            :   A vector of size `nc` (the number of clusters) containing the sizes of the individual clusters identified.
    !> \param[out]      Center          :   An array of size `(nd,nc)` containing the centers of the individual clusters identified.
    !> \param[out]      Membership      :   A vector of size `np` each element of which represents the cluster ID to which the point belongs.
    !> \param[out]      MinDistanceSq   :   An array of size `(np)` containing the minimum distances of the points from the cluster centers.
    !> \param[out]      potential       :   The sum of the distances-squared of all points from their corresponding group centers.
    !> \param[out]      niter           :   The number of iterations to reach convergence.
    !> \param[out]      nfail           :   The number of Zero-Sized Clusters Iterations.
    !> \param[out]      Err             :   An object of class [Err_type](@ref err_mod::err_type) containing information about error occurrence.
    !> \param[in]       relTol          :   If the **relative** difference between two subsequent potentials is below `relTol`, convergence is assumed (**optional**, default = 1.e-4).
    !> \param[in]       InitCenter      :   An array of size `(nd,nc)` representing the initial centers of the clusters (**optional**. If missing, Kmeans++ initializes `InitCenter`).
    !> \param[in]       niterMax    :   The maximum number of iterations beyond which it a lack of converge is assumed (**optional**, default = 1000).
    !> \param[in]       nfailMax   :   The maximum number of zero-sized cluster iterations beyond which it a lack of converge is assumed (**optional**, default = 10).
    !>
    !> \warning
    !> The input value for `niterMax` must be larger than 0 since at least one iteration will be performed by the procedure.
    !>
    !> \warning
    !> This algorithm does not check for the consistency of the input values for `nc`, `nd`, and `np`.
    !> It is currently the user's responsibility to ensure that these input values are logical and sensible.
    !>
    !> \author
    ! Amir Shahmoradi, April 03, 2017, 2:16 PM, ICES, UTEXAS
    subroutine runKmeans( nd, np, nc & ! LCOV_EXCL_LINE
                        , Point & ! LCOV_EXCL_LINE
                        , Size & ! LCOV_EXCL_LINE
                        , Center & ! LCOV_EXCL_LINE
                        , Membership & ! LCOV_EXCL_LINE
                       !, NormedPoint & ! LCOV_EXCL_LINE
                        , MinDistanceSq & ! LCOV_EXCL_LINE
                        , potential & ! LCOV_EXCL_LINE
                        , niter & ! LCOV_EXCL_LINE
                        , nfail & ! LCOV_EXCL_LINE
                        , Err & ! LCOV_EXCL_LINE
                        , relTol & ! LCOV_EXCL_LINE
                        , InitCenter & ! LCOV_EXCL_LINE
                        , nfailMax & ! LCOV_EXCL_LINE
                        , niterMax & ! LCOV_EXCL_LINE
                        )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: runKmeans
#endif
        use Constants_mod, only: IK, RK, HUGE_RK

        implicit none

        integer(IK)     , intent(in)            :: nc, nd, np
        real(RK)        , intent(in)            :: Point(nd,np)
        integer(IK)     , intent(out)           :: nfail
        integer(IK)     , intent(out)           :: niter
        real(RK)        , intent(out)           :: potential
        real(RK)        , intent(out)           :: Center(nd,nc)
       !real(RK)        , intent(out)           :: NormedPoint(nd,np,nc)
        real(RK)        , intent(out)           :: MinDistanceSq(np)
        integer(IK)     , intent(out)           :: Membership(np)
        integer(IK)     , intent(out)           :: Size(nc)
        type(Err_type)  , intent(out)           :: Err
        real(RK)        , intent(in), optional  :: relTol
        real(RK)        , intent(in), optional  :: InitCenter(nd,nc)
        integer(IK)     , intent(in), optional  :: nfailMax
        integer(IK)     , intent(in), optional  :: niterMax

        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME//"@runKmeans()"

        real(RK)                :: distanceSq               ! distance of each point ip from a given cluster center.
        real(RK)                :: potentialOld             ! The sum of all minimum distances squared.
        real(RK)                :: SumPoint(nd,nc)          ! The sum of the points in each cluster.
        integer(IK)             :: MembershipOld(np)        ! current (old) Memberships of the clusters.

        logical                 :: convergenceOccurred
        integer(IK)             :: nfailMaxDefault
        integer(IK)             :: niterMaxDefault
        real(RK)                :: relTolDefault
        integer                 :: ip, ic

        Err%occurred = .false.

        relTolDefault = 1.e-4_RK; if (present(relTol)) relTolDefault = relTol
        niterMaxDefault = 1e3_IK; if (present(niterMax)) niterMaxDefault = niterMax
        nfailMaxDefault = 1e1_IK; if (present(nfailMax)) nfailMaxDefault = nfailMax

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Find the distances of all points to all cluster centers. If not given, initialize the centers via the Kmeans++ algorithm.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if(present(InitCenter)) then
            !%%%%%%%%%%%%%%%%%%%
            ic = 1
            Size(ic) = 0_IK
            do ip = 1, np
                MinDistanceSq(ip) = sum( (InitCenter(1:nd,ic) - Point(1:nd,ip))**2 )
                MembershipOld(ip) = ic
            end do
            !%%%%%%%%%%%%%%%%%%%
            do ic = 2, nc - 1
                Size(ic) = 0_IK
                do ip = 1, np
                    distanceSq = sum( (InitCenter(1:nd,ic) - Point(1:nd,ip))**2 )
                    if (distanceSq < MinDistanceSq(ip)) then
                        MinDistanceSq(ip) = distanceSq
                        MembershipOld(ip) = ic
                    end if
                end do
            end do
            !%%%%%%%%%%%%%%%%%%%
            ic = nc
            Size(ic) = 0_IK
            SumPoint = 0._RK
            potentialOld = 0._RK
            do ip = 1, np
                distanceSq = sum( (InitCenter(1:nd,ic) - Point(1:nd,ip))**2 )
                if (distanceSq < MinDistanceSq(ip)) then
                    MinDistanceSq(ip) = distanceSq
                    MembershipOld(ip) = ic
                end if
                SumPoint(1:nd,MembershipOld(ip)) = SumPoint(1:nd,MembershipOld(ip)) + Point(1:nd,ip)
                Size(MembershipOld(ip)) = Size(MembershipOld(ip)) + 1_IK
                potentialOld = potentialOld + MinDistanceSq(ip)
            end do
            !%%%%%%%%%%%%%%%%%%%
        else
            call runKPP(nc, nd, np, Point, SumPoint, MembershipOld, Size, potentialOld)
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Refine all centers until convergence.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        nfail = 0_IK
        niter = 1_IK

        loopCenterRefinement: do

            ! compute the new cluster centers and create the new memberships.

            !%%%%%%%%%%%%%%%%%%%
            ic = 1
            Center(1:nd,ic) = SumPoint(1:nd,ic) / real(Size(ic), kind = RK)
            do ip = 1, np
                !NormedPoint(1:nd,ip,ic) = Center(1:nd,ic) - Point(1:nd,ip)
                !MinDistanceSq(ip) = sum( NormedPoint(1:nd,ip,ic)**2 )
                MinDistanceSq(ip) = sum( (Center(1:nd,ic) - Point(1:nd,ip))**2 )
                Membership(ip) = ic
            end do
            !%%%%%%%%%%%%%%%%%%%
            do ic = 2, nc - 1
                Center(1:nd,ic) = SumPoint(1:nd,ic) / real(Size(ic), kind = RK)
                do ip = 1, np
                   !NormedPoint(1:nd,ip,ic) = Center(1:nd,ic) - Point(1:nd,ip)
                   !distanceSq = sum( NormedPoint(1:nd,ip,ic)**2 )
                    distanceSq = sum( (Center(1:nd,ic) - Point(1:nd,ip))**2 )
                    if (distanceSq < MinDistanceSq(ip)) then
                        MinDistanceSq(ip) = distanceSq
                        Membership(ip) = ic
                    end if
                end do
            end do
            !%%%%%%%%%%%%%%%%%%%
            ic = nc
            convergenceOccurred = .true.
            potential = 0._RK
            Center(1:nd,ic) = SumPoint(1:nd,ic) / real(Size(ic), kind = RK)
            do ip = 1, np
               !NormedPoint(1:nd,ip,ic) = Center(1:nd,ic) - Point(1:nd,ip)
               !distanceSq = sum( NormedPoint(1:nd,ip,ic)**2 )
                distanceSq = sum( (Center(1:nd,ic) - Point(1:nd,ip))**2 )
                if (distanceSq < MinDistanceSq(ip)) then
                    MinDistanceSq(ip) = distanceSq
                    Membership(ip) = ic
                end if
                if (Membership(ip)/=MembershipOld(ip)) then
                    convergenceOccurred = .false.
                    Size(MembershipOld(ip)) = Size(MembershipOld(ip)) - 1_IK
                    Size(Membership(ip)) = Size(Membership(ip)) + 1_IK
                end if
                potential = potential + MinDistanceSq(ip)
            end do
            !%%%%%%%%%%%%%%%%%%%

            if (any(Size==0_IK)) then
                ! LCOV_EXCL_START

                niter = 1_IK
                nfail = nfail + 1_IK

                if (nfail > nfailMaxDefault) then
                    !Err%msg = "nfail > nfailMaxDefault"
                    Err%occurred = .true.
                    Err%stat = 2_IK
                    return
                end if

                ! find new random Centers and start over.

                call runKPP(nc, nd, np, Point, SumPoint, MembershipOld, Size, potentialOld)

                cycle loopCenterRefinement

#if defined DEBUG_ENABLED || defined TESTING_ENABLED
            elseif (any(Size<0_IK)) then

                block
                    use Err_mod, only: abort
                    Err%msg = PROCEDURE_NAME//": Internal error occurred - Size < 0. Please report this issue along with circumstances to the developers."
                    call abort(Err)
                    return
                end block
#endif

                ! LCOV_EXCL_STOP
            end if

            ! if `niterMax` has reached, quit.

            if (niter < niterMaxDefault) then

                ! If convergence has occurred, return.

                if (convergenceOccurred .or. abs(potential-potentialOld)/potential < relTolDefault) return

                niter = niter + 1_IK

                do ip = 1, np
                    SumPoint(1:nd,MembershipOld(ip)) = SumPoint(1:nd,MembershipOld(ip)) - Point(1:nd,ip)
                    SumPoint(1:nd,Membership(ip)) = SumPoint(1:nd,Membership(ip)) + Point(1:nd,ip)
                    MembershipOld(ip) = Membership(ip)
                end do
                potentialOld = potential

            else

                !Err%msg = "niter > niterMaxDefault"
                Err%occurred = .true.
                Err%stat = 1_IK
                return

            end if

        end do loopCenterRefinement

    end subroutine runKmeans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Find nc initial cluster centers for the Kmeans algorithm, using Kmeans++ recipe.
    !>
    !> \param[in]       nc              :   The number of clusters to be found via the Kmeans algorithm (must be > 1, though it is not verified).
    !> \param[in]       nd              :   The number of dimensions (features) of the input data array `Point`.
    !> \param[in]       np              :   The number of observations in the input data array `Point`.
    !> \param[in]       Point           :   The input array of size `(nd,np)` of the points to be clustered via the Kmeans algorithm.
    !> \param[out]      SumPoint        :   An output array of size `(nd,nc)` containing the sum of the points in each cluster.
    !>                                      When divided by the corresponding output cluster size, it will yield the final updated cluster centers.
    !> \param[out]      Membership      :   An output array of size `(np)` representing the IDs of the clusters to which individual points belong.
    !> \param[out]      Size            :   An output array of size `(nc)` representing the sizes of the identified clusters.
    !> \param[out]      potential       :   The sum of the minimum-distances-squared of all points from their corresponding group centers.
    !>
    !> \warning
    !> There is no benefit in calling this routine for `nc = 1`.
    !> As such, calling this routine with `nc = 1` will lead to undefined behavior.
    !>
    !> \todo
    !> All `minloc` queries in this procedure can be replaced with more efficient bin search algorithms.
    !>
    !> \todo
    !> The variable `CumSumDistSq` does not have to be defined as `real64`.
    !>
    !> \author
    ! Amir Shahmoradi, April 03, 2017, 2:16 PM, ICES, UTEXAS
    subroutine runKPP(nc, nd, np, Point, SumPoint, Membership, Size, potential)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: runKPP
#endif
        use Statistics_mod, only: getRandInt

        implicit none

        integer(IK) , intent(in)    :: np, nd, nc
        real(RK)    , intent(in)    :: Point(nd,np)
        real(RK)    , intent(out)   :: SumPoint(nd,nc)
        integer(IK) , intent(out)   :: Membership(np)
        integer(IK) , intent(out)   :: Size(nc)
        real(RK)    , intent(out)   :: potential

        real(RK)                    :: dummy
        real(RK)                    :: distanceSq
        real(RK)                    :: Center(nd,nc)
        real(RK)                    :: MinDistanceSq(np)    ! An array representing the minimum distances of individual points from all clusters.
        real(RK)                    :: CumSumDistSq(0:np)   ! sum of distance squared of points up to the given index of the vector.
        integer(IK)                 :: ic, ip

        CumSumDistSq(0) = 0._RK ! This element must always remain zero.

        ic = 1_IK
        Center(1:nd,ic) = Point( 1:nd , getRandInt(1,np) ) ! perhaps this random assignment is unnecessary, after all, Point is supposed to be of random order.
        do ip = 1, np
            MinDistanceSq(ip) = sum( (Center(1:nd,ic) - Point(1:nd,ip))**2 ) ! find the distance-squared of each point to the # ic cluster Center.
            CumSumDistSq(ip) = CumSumDistSq(ip-1) + MinDistanceSq(ip)
            Membership(ip) = ic
        end do

        ic = 2_IK
        call random_number(dummy)
        dummy = dummy * CumSumDistSq(np)
        Center(1:nd,ic) = Point( 1:nd, minloc(CumSumDistSq, dim = 1, mask = CumSumDistSq > dummy) - 1 ) ! -1 is essential to compensate for the 0 starting index of CumSumDistSq

        ! find the rest of cluster centers

        do ic = 2, nc - 1

            do ip = 1, np
                distanceSq = sum( (Center(1:nd,ic) - Point(1:nd,ip))**2 )
                if (distanceSq < MinDistanceSq(ip)) then
                    MinDistanceSq(ip) = distanceSq
                    Membership(ip) = ic
                end if
                CumSumDistSq(ip) = CumSumDistSq(ip-1) + MinDistanceSq(ip)
            end do

            call random_number(dummy)
            dummy = dummy * CumSumDistSq(np)
            Center(1:nd,ic+1) = Point( 1:nd , minloc(CumSumDistSq, dim = 1, mask = CumSumDistSq > dummy) - 1 )

        end do

        Size = 0_IK
        SumPoint = 0._RK
        potential = 0._RK
        do ip = 1, np
            distanceSq = sum( (Center(1:nd,nc) - Point(1:nd,ip))**2 )
            if (distanceSq < MinDistanceSq(ip)) then
                MinDistanceSq(ip) = distanceSq
                Membership(ip) = nc
            end if
            SumPoint(1:nd,Membership(ip)) = SumPoint(1:nd,Membership(ip)) + Point(1:nd,ip)
            Size(Membership(ip)) = Size(Membership(ip)) + 1_IK
            potential = potential + MinDistanceSq(ip)
        end do

    end subroutine runKPP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine writeKmeans(Kmeans, Point, nd, np, fileUnit)
        use Constants_mod, only: IK, RK, NLC
        implicit none
        type(Kmeans_type)   , intent(in)    :: Kmeans
        integer(IK)         , intent(in)    :: nd, np
        integer(IK)         , intent(in)    :: fileUnit
        real(RK)            , intent(in)    :: Point(nd,np)
        integer(IK)                         :: i, j, ic, nc
        character(*)        , parameter     :: fileFormat = "(*(g0.8,:,','))"

        nc = size(Kmeans%Size)

        write(fileUnit,"(*(g0.8,:,'"//NLC//"'))") "nd", nd, "np", np, "nc", nc, "niter", Kmeans%niter, "nfail", Kmeans%nfail

        if (allocated(Kmeans%Size)) then
            write(fileUnit,"(A)") "Size"
            write(fileUnit,fileFormat) Kmeans%Size
        end if

        if (allocated(Kmeans%Center)) then
            write(fileUnit,"(A)") "Center"
            write(fileUnit,fileFormat) Kmeans%Center
        end if

        if (allocated(Kmeans%LogVol)) then
            write(fileUnit,"(A)") "LogVol"
            write(fileUnit,fileFormat) Kmeans%LogVol
        end if

        if (allocated(Kmeans%ScaleFactorSq)) then
            write(fileUnit,"(A)") "ScaleFactorSq"
            write(fileUnit,fileFormat) Kmeans%ScaleFactorSq
        end if

        if (allocated(Kmeans%ChoLowCovUpp)) then
            write(fileUnit,"(A)") "CholeskyLower"
            write(fileUnit,fileFormat) ((Kmeans%ChoDia(j,ic), (Kmeans%ChoLowCovUpp(i,j,ic), i=j+1,nd), j=1,nd), ic=1,nc)
        end if

        write(fileUnit,"(A)") "Point"
        write(fileUnit,fileFormat) Point

        if (allocated(Kmeans%MahalSq)) then
            write(fileUnit,"(A)") "MahalSq"
            write(fileUnit,fileFormat) Kmeans%MahalSq
        end if

        if (allocated(Kmeans%Membership)) then
            write(fileUnit,"(A)") "Membership"
            write(fileUnit,fileFormat) Kmeans%Membership
        end if

        if (allocated(Kmeans%MinDistanceSq)) then
            write(fileUnit,"(A)") "MinDistanceSq"
            write(fileUnit,fileFormat) Kmeans%MinDistanceSq
        end if

    end subroutine writeKmeans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Kmeans_mod ! LCOV_EXCL_LINE