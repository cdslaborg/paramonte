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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> The `Prop_type` class, containing the properties of Kmeans clusters.
    type :: Prop_type
        real(RK)                    :: logSumVolNormed      !< The logarithm of the sum of the bounding volumes of the clusters *minus their overlaps*.
       !real(RK)                    :: logAvgDenNormed      !< The logarithm of the average of the densities of points in the bounding volumes of the clusters.
        real(RK)    , allocatable   :: ChoDia(:,:)          !< An array of size `(nd,nc)` representing the diagonal elements of the Cholesky factorization of the covariance matrix.
        real(RK)    , allocatable   :: MahalSq(:,:)         !< An array of size `(np,nc)` each element of which represents the Mahalanobis distance squared of point `ip` from the cluster `ic`.
        real(RK)    , allocatable   :: InvCovMat(:,:,:)     !< An array of size `(nd,nd,nc)` containing the inverse of the covariance matrix of the corresponding cluster.
        real(RK)    , allocatable   :: LogVolNormed(:)      !< An array of size `(nc)` each element of which represents the log(volume) of the bounding covariance matrix of the corresponding cluster.
       !real(RK)    , allocatable   :: LogDenNormed(:)      !< An array of size `(nc)` each element of which represents the log effective (mean) density of points in the corresponding cluster bounding volume.
        real(RK)    , allocatable   :: ChoLowCovUpp(:,:,:)  !< An array of size `(nd,nd,nc)` whose upper triangle and diagonal is the covariance matrix and the lower is the Cholesky Lower.
        real(RK)    , allocatable   :: ScaleFactorSq(:)     !< An array of size `(nc)` representing the factors by which the cluster covariance matrices must be enlarged to enclose their corresponding members.
        integer(IK) , allocatable   :: EffectiveSize(:)     !< An array of size `(nc)` representing the factors by which the cluster covariance matrices must be enlarged to enclose their corresponding members.
        integer(IK) , allocatable   :: CumSumSize(:)        !< A vector of size `0:nc` representing the cumulative sum of all cluster sizes from cluster 1 to the last.
        integer(IK) , allocatable   :: Index(:)             !< A vector of size `(np)` such that the output `Point(:,Index(ip),:)` points to the input `Point(:,Index(ip),:)` in [getProp](@ref getprop).
                                                            !< The `Index` component is populated only if `Index` is missing from the list of input arguments to [getProp](@ref getprop).
    contains
        procedure, pass :: allocate => allocateKmeansProp
    end type Prop_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> The `Kmeans_type` class.
    !> The inclusion of the component `NormedPoint` adds ~50% to the computational cost of performing Kmeans.
    !> In this latest implementation, avoiding `NormedPoint` as a component leads to 15-50% performance improvement.
    !> The moral is, simple repeating of addition operation is much significantly cheaper than storing the operation results for later usage.
    !> In the following declaration descriptions,
    !>     +    `nc` represents the number of clusters to be found.
    !>     +    `np` represents the number of data points to be clustered.
    !>     +    `nd` represents the number of dimensions (features) of the points to be clustered.
    type :: Kmeans_type
        integer(IK)                 :: nfail                !< The number of Zero-Sized Clusters Iterations.
        integer(IK)                 :: niter                !< The number of iterations to reach convergence.
        real(RK)                    :: potential            !< The sum of the distances-squared of all points from their corresponding group centers.
        real(RK)    , allocatable   :: Center(:,:)          !< An array of size `(nd,nc)` containing the centers of the individual clusters identified.
       !real(RK)    , allocatable   :: NormedPoint(:,:,:)   !< An array of size `(nd,np,nc)` containing `nc` copies of the input array `Point(nd,np)` normalized to the cluster centers.
        real(RK)    , allocatable   :: MinDistanceSq(:)     !< An array of size `(np)` of the minimum distances of the points from the cluster centers.
        integer(IK) , allocatable   :: Membership(:)        !< A vector of size `np` each element of which represents the cluster ID to which the point belongs.
        integer(IK) , allocatable   :: Size(:)              !< A vector of size `nc` (the number of clusters) containing the sizes of the individual clusters identified.
        type(Prop_type)             :: Prop                 !< An object of class [Prop_type](@ref prop_type) containing the Kmeans cluster properties.
        type(Err_type)              :: Err                  !< An object of class [Err_type](@ref err_mod::err_type) containing information about error occurrence.
    contains
        procedure, pass :: getProp
        !procedure, pass :: pack => packKmeans
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

    pure subroutine allocateKmeansProp(Prop, nd, np, nc)
        implicit none
        class(Prop_type), intent(inout) :: Prop
        integer(IK), intent(in)         :: nd, np, nc
        if (.not. allocated(Prop%ChoDia         )) allocate(Prop%ChoDia         (nd,nc))
        if (.not. allocated(Prop%MahalSq        )) allocate(Prop%MahalSq        (np,nc))
        if (.not. allocated(Prop%InvCovMat      )) allocate(Prop%InvCovMat      (nd,nd,nc))
        if (.not. allocated(Prop%ChoLowCovUpp   )) allocate(Prop%ChoLowCovUpp   (nd,nd,nc))
        if (.not. allocated(Prop%EffectiveSize  )) allocate(Prop%EffectiveSize  (nc))
        if (.not. allocated(Prop%ScaleFactorSq  )) allocate(Prop%ScaleFactorSq  (nc))
        if (.not. allocated(Prop%LogVolNormed   )) allocate(Prop%LogVolNormed   (nc))
       !if (.not. allocated(Prop%LogDenNormed   )) allocate(Prop%LogDenNormed   (nc))
        if (.not. allocated(Prop%CumSumSize     )) allocate(Prop%CumSumSize     (0:nc))
        if (.not. allocated(Prop%Index          )) allocate(Prop%Index          (np))
    end subroutine allocateKmeansProp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Perform the Kmeans clustering for `nt` tries on the input data set represented by the array `Point(nd, np)`
    !> and return the clustering that yields the least potential.
    !>
    !> \param[in]       nd                  :   The number of dimensions (attributes) in the input `Point`.
    !> \param[in]       np                  :   The number of observations in the input `Point`.
    !> \param[in]       nc                  :   The number of clusters to identify.
    !> \param[in]       nt                  :   The number of independent instances of Kmeans to run.
    !> \param[in]       Point               :   An array of size `(nd,np)` containing the input dataset to cluster.
    !> \param[in]       InitCenter          :   An array of size `(nd,nc)` containing the best-guess cluster centers (**optional**, default = set via Kpp).
    !> \param[in]       niterMax            :   The maximum allowed number of iterations in the algorithm of Lloyd (**optional**, default = `300`).
    !> \param[in]       nfailMax            :   The maximum allowed number of failures in the algorithm of Lloyd (**optional**, default = `2`).
    !> \param[in]       minSize             :   The minimum allowed size of a cluster (**optional**, default = `0`).
    !> \param[in]       relTol              :   See the description of the [runKmeans()](@ref runkmeans).
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
    !> It is currently the responsibility of the user to ensure that these input values are logical and sensible.
    !>
    !> \author
    ! Amir Shahmoradi, April 03, 2017, 2:16 PM, ICES, UTEXAS
    function getKmeans(nd, np, nc, nt, Point, InitCenter, relTol, niterMax, nfailMax, minSize) result(Kmeans)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getKmeans
#endif
        use Constants_mod, only: IK, RK, HUGE_RK

        implicit none

        integer(IK) , intent(in)                :: nd, np, nc, nt
        real(RK)    , intent(in)                :: Point(nd,np)
        real(RK)    , intent(in)    , optional  :: InitCenter(nd,nc)
        real(RK)    , intent(in)    , optional  :: relTol
        integer(IK) , intent(in)    , optional  :: minSize
        integer(IK) , intent(in)    , optional  :: nfailMax
        integer(IK) , intent(in)    , optional  :: niterMax

        type(Kmeans_type)                       :: Kmeans
        type(Kmeans_type), allocatable          :: KmeansTry(:)

        integer(IK)                             :: itry
        integer(IK)                             :: itryLeastPotential
        integer(IK)                             :: nfailMaxDefault
        integer(IK)                             :: niterMaxDefault
        integer(IK)                             :: minSizeDefault
        real(RK)                                :: relTolDefault
        real(RK)                                :: potential

        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME//"@getKmeans()"

        niterMaxDefault = 1e3_IK    ; if (present(niterMax)) niterMaxDefault = niterMax
        nfailMaxDefault = 1e1_IK    ; if (present(nfailMax)) nfailMaxDefault = nfailMax
        minSizeDefault  = 1_IK      ; if (present(minSize)) minSizeDefault = minSize
        relTolDefault   = 1.e-4_RK  ; if (present(relTol)) relTolDefault = relTol


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
                                , nfailMax = nfailMaxDefault & ! LCOV_EXCL_LINE
                                , niterMax = niterMaxDefault & ! LCOV_EXCL_LINE
                                , minSize = nfailMaxDefault & ! LCOV_EXCL_LINE
                                , relTol = relTolDefault & ! LCOV_EXCL_LINE
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
                                , InitCenter = InitCenter & ! LCOV_EXCL_LINE
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
                            , niterMax = niterMaxDefault & ! LCOV_EXCL_LINE
                            , nfailMax = nfailMaxDefault & ! LCOV_EXCL_LINE
                            , minSize = minSizeDefault & ! LCOV_EXCL_LINE
                            , relTol = relTolDefault & ! LCOV_EXCL_LINE
                            , Point = Point & ! LCOV_EXCL_LINE
                            , Size = Kmeans%Size & ! LCOV_EXCL_LINE
                            , Center = Kmeans%Center & ! LCOV_EXCL_LINE
                            , Membership = Kmeans%Membership & ! LCOV_EXCL_LINE
                            , MinDistanceSq = Kmeans%MinDistanceSq & ! LCOV_EXCL_LINE
                            , potential = Kmeans%potential & ! LCOV_EXCL_LINE
                            , niter = Kmeans%niter & ! LCOV_EXCL_LINE
                            , nfail = Kmeans%nfail & ! LCOV_EXCL_LINE
                            , Err = Kmeans%Err & ! LCOV_EXCL_LINE
                            , InitCenter = InitCenter & ! LCOV_EXCL_LINE
                            )
#if defined DEBUG_ENABLED || defined TESTING_ENABLED
            if (any(Kmeans%Size < minSizeDefault)) then ! @todo xxx : This need further enhancement. Is it even necessary?
                ! LCOV_EXCL_START
                if (.not. Kmeans%Err%occurred .and. Kmeans%Err%stat /= 2_IK) then
                    Kmeans%Err%msg = "The Kmeans clustering sizes are smaller than the requested input minimum cluster size."
                    write(*,*) Kmeans%Err%msg, minSizeDefault, Kmeans%Size
                    Kmeans%Err%occurred = .true.
                    Kmeans%Err%stat = 1_IK
                    error stop
                    !return
                end if
                ! LCOV_EXCL_STOP
            else
                block
                    integer(IK) :: ic
                    do ic = 1, nc
                        if (Kmeans%Size(ic) == 0_IK) then
                            if (any(Kmeans%Membership(1:np) == ic)) then
                                write(*,*) "Internal error occurred: Zero-sized cluster has associated members. ic, ip = ", ic
                                error stop
                            end if
                        end if
                    end do
                end block
            end if
#endif

        end if

    end function getKmeans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Perform the Kmeans clustering on the input data set represented by the array `Point(nd, np)`.
    !>
    !> \param[in]       nd              :   The number of dimensions (features) of the input data array `Point`.
    !> \param[in]       np              :   The number of observations in the input data array `Point`.
    !> \param[in]       nc              :   The number of clusters to be found via the Kmeans algorithm.
    !> \param[in]       nfailMax        :   A non-negative integer, representing the maximum allowed number of minSize-cluster iterations beyond which it a lack of converge is assumed.
    !> \param[in]       niterMax        :   A non-negative integer, representing the maximum allowed number of iterations beyond which it a lack of converge is assumed.
    !> \param[in]       minSize         :   A non-negative integer, representing the minimum allowed number of points in each cluster.
    !> \param[in]       relTol          :   If the **relative** difference between two subsequent potentials is below `relTol`, convergence is assumed.
    !> \param[in]       Point           :   The input array of size `(nd,np)` of the points to be clustered via the Kmeans algorithm.
    !> \param[out]      Size            :   A vector of size `nc` (the number of clusters) containing the sizes of the individual clusters identified.
    !> \param[out]      Center          :   An array of size `(nd,nc)` containing the centers of the individual clusters identified.
    !> \param[out]      Membership      :   A vector of size `np` each element of which represents the cluster ID to which the point belongs.
    !> \param[out]      MinDistanceSq   :   An array of size `(np)` containing the minimum distances of the points from the cluster centers.
    !> \param[out]      potential       :   The sum of the distances-squared of all points from their corresponding group centers.
    !> \param[out]      niter           :   The number of iterations to reach convergence.
    !> \param[out]      nfail           :   The number of iterations for which clusters with size smaller than `minSize` occurs.
    !> \param[out]      Err             :   An object of class [Err_type](@ref err_mod::err_type) containing information about error occurrence.
    !> \param[in]       InitCenter      :   An array of size `(nd,nc)` representing the initial centers of the clusters (**optional**. If missing, Kmeans++ initializes `InitCenter`).
    !>
    !> \warning
    !> The input value for `niterMax` must be larger than 0 since at least one iteration will be performed by the procedure.
    !>
    !> \warning
    !> This algorithm does not check for the consistency of the input values for `nc`, `nd`, and `np`.
    !> It is currently the responsibility of the user to ensure that these input values are logical and sensible.
    !>
    !> \remark
    !> If `niterMax = 0`, then no refinement of the centers based on the algorithm of Lloyd will be performed.
    !> In such a case, the cluster memberships will be determined based on the initial input `InitCenter` or if not provided as input,
    !> base on the randomly-guessed `Center` using Kpp algorithm and the procedure will return before further refining the `Membership` and `Center`.
    !> It is currently the responsibility of the user to ensure that these input values are logical and sensible.
    !>
    !> \author
    ! Amir Shahmoradi, April 03, 2017, 2:16 PM, ICES, UTEXAS
    subroutine runKmeans( nd, np, nc & ! LCOV_EXCL_LINE
                        , nfailMax & ! LCOV_EXCL_LINE
                        , niterMax & ! LCOV_EXCL_LINE
                        , minSize & ! LCOV_EXCL_LINE
                        , relTol & ! LCOV_EXCL_LINE
                        , Point & ! LCOV_EXCL_LINE
                        , Size & ! LCOV_EXCL_LINE
                        , niter & ! LCOV_EXCL_LINE
                        , nfail & ! LCOV_EXCL_LINE
                        , Center & ! LCOV_EXCL_LINE
                        , potential & ! LCOV_EXCL_LINE
                        , Membership & ! LCOV_EXCL_LINE
                        , MinDistanceSq & ! LCOV_EXCL_LINE
                        , Err & ! LCOV_EXCL_LINE
                        , InitCenter & ! LCOV_EXCL_LINE
                        )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: runKmeans
#endif
        use Constants_mod, only: IK, RK, HUGE_RK

        implicit none

        integer(IK)     , intent(in)            :: nc, nd, np, minSize, nfailMax, niterMax
        real(RK)        , intent(in)            :: Point(nd,np)
        real(RK)        , intent(in)            :: relTol
        integer(IK)     , intent(out)           :: nfail
        integer(IK)     , intent(out)           :: niter
        real(RK)        , intent(out)           :: potential
        real(RK)        , intent(out)           :: Center(nd,nc)
       !real(RK)        , intent(out)           :: NormedPoint(nd,np,nc)
        real(RK)        , intent(out)           :: MinDistanceSq(np)
        integer(IK)     , intent(out)           :: Membership(np)
        integer(IK)     , intent(out)           :: Size(nc)
        type(Err_type)  , intent(out)           :: Err
        real(RK)        , intent(in), optional  :: InitCenter(nd,nc)

        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME//"@runKmeans()"

        integer(IK)                             :: MembershipOld(np)        ! current (old) Memberships of the clusters.
        logical                                 :: convergenceOccurred
        integer                                 :: ip, ic
        real(RK)                                :: distanceSq               ! distance of each point ip from a given cluster center.
        real(RK)                                :: potentialOld             ! The sum of all minimum distances squared.
        real(RK)                                :: SumPoint(nd,nc)          ! The sum of the points in each cluster.

        Err%occurred = .false.

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
        niter = 0_IK

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

#if defined DEBUG_ENABLED || defined TESTING_ENABLED
            if (any(Size<0_IK)) then
                ! LCOV_EXCL_START
                block
                    use Err_mod, only: abort
                    Err%msg = PROCEDURE_NAME//": Internal error occurred - Size < 0. Please report this issue along with circumstances to the developers at: https://github.com/cdslaborg/paramonte/issues"
                    call abort(Err)
                    return
                end block
                ! LCOV_EXCL_STOP
            end if
#endif

            if (any(Size<minSize)) then
                ! LCOV_EXCL_START
                niter = 0_IK
                nfail = nfail + 1_IK
                if (nfail > nfailMax) then
                    !Err%msg = "nfail > nfailMax"
                    Err%occurred = .true.
                    Err%stat = 2_IK
                    return
                end if
                call runKPP(nc, nd, np, Point, SumPoint, MembershipOld, Size, potentialOld) ! find new random Centers and start over.
                cycle loopCenterRefinement
                ! LCOV_EXCL_STOP
            end if

            ! if `niterMax` has reached, quit.

            if (niter < niterMax) then

                ! If convergence has occurred, return.

                if (convergenceOccurred .or. abs(potential-potentialOld)/potential < relTol) return

                niter = niter + 1_IK

                do ip = 1, np
                    SumPoint(1:nd,MembershipOld(ip)) = SumPoint(1:nd,MembershipOld(ip)) - Point(1:nd,ip)
                    SumPoint(1:nd,Membership(ip)) = SumPoint(1:nd,Membership(ip)) + Point(1:nd,ip)
                    MembershipOld(ip) = Membership(ip)
                end do
                potentialOld = potential

            else

                !Err%msg = "niter > niterMax"
                Err%occurred = .true.
                Err%stat = 1_IK
                return

            end if

        end do loopCenterRefinement

    end subroutine runKmeans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> This is a method of class [Kmeans_type](@ref kmeans_type).
    !> Compute the following properties (components) of the input object of class [Kmeans_type](@ref kmeans_type):
    !>
    !> +   `ChoLowCovUpp(1:nd,1:nd,1:nc)`   :   The covariance matrix upper triangle and Cholesky factor lower triangle,
    !> +   `InvCovMat(1:nd,1:nd,1:nc)`      :   The full inverse covariance matrix,
    !> +   `EffectiveSize(1:nc)`            :   The effective sizes of the ellipsoidal bounded regions of the clusters.
    !> +   `ScaleFactorSq(1:nc)`            :   The factor by which the covariance matrix elements must be multiplied
    !>                                          to convert the matrix to a bounding ellipsoid of the cluster members.
    !> +   `MahalSq(1:np,1:nc)`             :   The Mahalanobis distances squared of all points from all cluster centers,
    !> +   `ChoDia(1:nd,1:nc)`              :   The diagonal elements of the Cholesky lower triangle,
    !> +   `CumSumSize(0:nc)`               :   The cumulative sum of `Kmeans%Size(1:nc)` with `CumSumSize(0) = 0`,
    !> +   `LogVolNormed(1:nc)`             :   The vector of natural-log volumes of all clusters normalized by the volume of unit nd-ball,
    !>
    !> \param[in]       nd                  :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[in]       np                  :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[inout]    Point               :   The Array of size `(nd,np)` representing the original points used in Kmeans clustering.
    !> \param[inout]    Index               :   A vector of size `np` such that the output `Point(:,Index(ip),:)` points to the input `Point(:,ip)`
    !>                                          (**optional**, if not provided, it will be stored in `Prop%Index` component of the input `Kmeans` object).
    !> \param[in]       inclusionFraction   :   The fraction of points inside the bounding region of each cluster to be considered in computing the effective sizes of
    !>                                          the bounded regions of each cluster (**optional**, default = `0.`).
    !> \param[in]       pointLogVolNormed   :   The logarithm of the volume of a single point, to be used for assigning the properties of singular clusters.
    !>                                          (**optional**, if not provided or if set to `NEGINF_RK`, the volume will be estimated via the other non-singular clusters).
    !>
    !> \warning
    !> If any cluster has less than `nd + 1` members, its properties will not be computed, unless `pointLogVol` is provided as input.
    !> When `pointLogVol` is provided, the properties of any non-zero-size cluster will be computed as hyper-sphere whose volume is
    !> the volume of a single point times the number of cluster members. Obviously, the properties of cluster of zero size will not be set.
    !>
    !> \warning
    !> If `inclusionFraction = 0.`, then all clusters are effective assumed to be non-overlapping.
    !> This issue becomes particularly relevant when computing the effective sum of volumes and
    !> average density across all clusters.
    !>
    !> \warning
    !> Upon entry to this method, the `Center`, `Size`, and `Membership` components of the Kmeans object
    !> must have been already properly set up by its class constructor or a call to `getKmeans()`.
    !>
    !> \warning
    !> On output, the arrays `Point` and `Index` will be rearranged on exit from the procedure such that the members of the same
    !> clusters appear next to each other contiguously in the array. Consequently, the output array component `Index(1:np)` will
    !> be populated such that the input `Point(1:nd,Index(ip))` is saved at the output `Point(1:nd,ip)`.
    !>
    !> \warning
    !> If the Cholesky factorization in the procedure fails, the control is returned with `Kmeans%Err%occurred = .true.`.
    !>
    !> \warning
    !> If `pointLogVolNormed` is provided on input, it will be solely used to define the properties of singular clusters.
    !> Consequently, the densities of singular clusters will be the same while the densities of non-singular clusters will
    !> be different from each other and from those of singular clusters.
    subroutine getProp(Kmeans, nd, np, Point, Index, inclusionFraction, pointLogVolNormed)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProp
#endif
        use Constants_mod, only: IK, RK, HUGE_RK, NEGINF_RK
        use Matrix_mod, only: getInvMatFromCholFac, getEye
        use Matrix_mod, only: getCholeskyFactor
        use Math_mod, only: getLogSumExp

        implicit none

        class(Kmeans_type) , intent(inout)      :: Kmeans
        integer(IK) , intent(in)                :: nd, np
        real(RK)    , intent(inout)             :: Point(nd,np)
        real(RK)    , intent(in)    , optional  :: pointLogVolNormed
        real(RK)    , intent(in)    , optional  :: inclusionFraction
        integer(IK) , intent(inout) , optional  :: Index(np)

        integer(IK) , allocatable               :: KmeansMemberCounter(:)
        integer(IK)                             :: i, j, ic, ip, nc
        integer(IK)                             :: ipstart, ipend

        real(RK)                                :: NormedPoint(nd,np)
        real(RK)                                :: scaleFactorSqInverse
        real(RK)                                :: pointLogVolNormedDefault
        real(RK)                                :: maxLogVolNormed
        real(RK)                                :: scaleFactor

        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME//"@getProp()"

        nc = size(Kmeans%Size)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Reorder Point based on the identified clusters.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        call Kmeans%Prop%allocate(nd,np,nc)

        if (allocated(KmeansMemberCounter)) deallocate(KmeansMemberCounter) ! LCOV_EXCL_LINE
        allocate(KmeansMemberCounter(nc))

        ! Compute the cumulative sum of the cluster Size array.

        Kmeans%Prop%CumSumSize(0) = 0_IK
        loopDetermineClusterBoundaries: do ic = 1, nc
            KmeansMemberCounter(ic) = Kmeans%Prop%CumSumSize(ic-1)
            Kmeans%Prop%CumSumSize(ic) = KmeansMemberCounter(ic) + Kmeans%Size(ic)
        end do loopDetermineClusterBoundaries

        if (present(Index)) then
            do ip = 1, np
                KmeansMemberCounter(Kmeans%Membership(ip)) = KmeansMemberCounter(Kmeans%Membership(ip)) + 1_IK
                NormedPoint(1:nd,KmeansMemberCounter(Kmeans%Membership(ip))) = Point(1:nd,ip)
                Kmeans%Prop%Index(KmeansMemberCounter(Kmeans%Membership(ip))) = Index(ip)
            end do
            Index = Kmeans%Prop%Index
        else
            do ip = 1, np
                KmeansMemberCounter(Kmeans%Membership(ip)) = KmeansMemberCounter(Kmeans%Membership(ip)) + 1_IK
                NormedPoint(1:nd,KmeansMemberCounter(Kmeans%Membership(ip))) = Point(1:nd,ip)
                Kmeans%Prop%Index(KmeansMemberCounter(Kmeans%Membership(ip))) = ip
            end do
        end if
        Point = NormedPoint

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the properties of non-singular clusters.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        maxLogVolNormed = NEGINF_RK
        pointLogVolNormedDefault = -1._RK ! indicator
        Kmeans%Prop%logSumVolNormed = NEGINF_RK
        loopComputeClusterProperties: do ic = 1, nc

            blockMinimumClusterSize: if (Kmeans%Size(ic) > nd) then

                ipstart = Kmeans%Prop%CumSumSize(ic-1) + 1_IK
                ipend = Kmeans%Prop%CumSumSize(ic)

                ! Correct the cluster memberships.

                Kmeans%Membership(ipstart:ipend) = ic

                ! Normalize points.

                do concurrent(ip = 1:np)
                    NormedPoint(1:nd,ip) = Point(1:nd,ip) - Kmeans%Center(1:nd,ic)
                end do

                ! Compute the upper covariance matrix of the cluster covariance matrices.

                Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic) = 0._RK
                do ip = ipstart, ipend
                    do j = 1, nd
                        Kmeans%Prop%ChoLowCovUpp(1:j,j,ic) = Kmeans%Prop%ChoLowCovUpp(1:j,j,ic) + NormedPoint(1:j,ip) * NormedPoint(j,ip)
                    end do
                end do
                Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic) = Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic) / real(Kmeans%Size(ic)-1_IK, kind = RK)

                ! Compute the Cholesky Factor of the cluster covariance matrices.

                call getCholeskyFactor(nd, Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic), Kmeans%Prop%ChoDia(1:nd,ic))
                if (Kmeans%Prop%ChoDia(1,ic) < 0._RK) then
                    ! LCOV_EXCL_START
                    Kmeans%Err%msg = PROCEDURE_NAME//"Cholesky factorization failed."
                    Kmeans%Err%occurred = .true.
                    return
                    ! LCOV_EXCL_STOP
                end if

                ! Compute the inverse of the cluster covariance matrices.

                Kmeans%Prop%InvCovMat(1:nd,1:nd,ic) = getInvMatFromCholFac( nd, Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic), Kmeans%Prop%ChoDia(1:nd,ic) )

                ! Compute the MahalSq of as many points as needed.

                do concurrent(ip = 1:np)
                    Kmeans%Prop%MahalSq(ip,ic) = dot_product( NormedPoint(1:nd,ip) , matmul(Kmeans%Prop%InvCovMat(1:nd,1:nd,ic), NormedPoint(1:nd,ip)) )
                end do

                ! Compute the scaleFcator of the bounding region and scale the properties to form the bounded regions.

                Kmeans%Prop%ScaleFactorSq(ic) = maxval(Kmeans%Prop%MahalSq(ipstart:ipend,ic))

                scaleFactorSqInverse = 1._RK / Kmeans%Prop%ScaleFactorSq(ic)
                scaleFactor = sqrt(Kmeans%Prop%ScaleFactorSq(ic))
                Kmeans%Prop%ChoDia(1:nd,ic) = Kmeans%Prop%ChoDia(1:nd,ic) * scaleFactor
                Kmeans%Prop%InvCovMat(1:nd,1:nd,ic) = Kmeans%Prop%InvCovMat(1:nd,1:nd,ic) * scaleFactorSqInverse
                do j = 1, nd
                    do i = j + 1, nd
                        Kmeans%Prop%ChoLowCovUpp(i,j,ic) = Kmeans%Prop%ChoLowCovUpp(i,j,ic) * scaleFactor
                    end do
                end do

                ! Rescale the MahalSq to the bounding regions

                !Kmeans%Prop%MahalSq(1:np,ic) = Kmeans%Prop%MahalSq(1:np,ic) * scaleFactorSqInverse

                ! Compute the effective fraction of points inside the bounded region.

                if (present(inclusionFraction)) then
                    Kmeans%Prop%EffectiveSize(ic)   = Kmeans%Size(ic) & ! LCOV_EXCL_LINE
                                                    + nint(inclusionFraction * &
                                                    ( count(Kmeans%Prop%MahalSq(1:Kmeans%Prop%CumSumSize(ic-1),ic)<=Kmeans%Prop%ScaleFactorSq(ic)) & ! LCOV_EXCL_LINE
                                                    + count(Kmeans%Prop%MahalSq(Kmeans%Prop%CumSumSize(ic)+1:np,ic)<=Kmeans%Prop%ScaleFactorSq(ic)) & ! LCOV_EXCL_LINE
                                                    ), kind=IK)
!write(*,*) 
!write(*,*) "Kmeans%ic, CumSumSize(ic-1), CumSumSize(ic), np", ic, Kmeans%Prop%CumSumSize(ic-1), Kmeans%Prop%CumSumSize(ic), np
!write(*,*) Kmeans%Size(ic)
!write(*,*) count(Kmeans%Prop%MahalSq(1:Kmeans%Prop%CumSumSize(ic-1),ic)<Kmeans%Prop%ScaleFactorSq(ic))
!write(*,*) count(Kmeans%Prop%MahalSq(Kmeans%Prop%CumSumSize(ic)+1:np,ic)<Kmeans%Prop%ScaleFactorSq(ic))
!write(*,*) "Kmeans%Prop%ScaleFactorSq(ic)", Kmeans%Prop%ScaleFactorSq(ic)
!write(*,*) "Kmeans%Prop%MahalSq(Kmeans%Prop%CumSumSize(ic)+1:np,ic): ", Kmeans%Prop%MahalSq(Kmeans%Prop%CumSumSize(ic)+1:np,ic)
!write(*,*) 
                else
                    Kmeans%Prop%EffectiveSize(ic) = Kmeans%Size(ic)
                end if

                Kmeans%Prop%LogVolNormed(ic) = sum( log(Kmeans%Prop%ChoDia(1:nd,ic)) )
                maxLogVolNormed = max( maxLogVolNormed, Kmeans%Prop%LogVolNormed(ic) )
               !Kmeans%Prop%LogDenNormed(ic) = log(real(Kmeans%Prop%EffectiveSize(ic),RK)) - Kmeans%Prop%LogVolNormed(ic)

            else blockMinimumClusterSize

                Kmeans%Prop%LogVolNormed(ic) = NEGINF_RK ! This initialization is essential in dependent procedures.
                Kmeans%Prop%EffectiveSize(ic) = 0_IK
                pointLogVolNormedDefault = 1._RK ! This is an indicator

            end if blockMinimumClusterSize

        end do loopComputeClusterProperties

        if (maxLogVolNormed /= NEGINF_RK) Kmeans%Prop%logSumVolNormed = getLogSumExp(nc, Kmeans%Prop%LogVolNormed, maxLogVolNormed)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Compute the properties of singular clusters.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (pointLogVolNormedDefault > 0._RK) then

            if (present(pointLogVolNormed)) then
                if (pointLogVolNormed > NEGINF_RK) pointLogVolNormedDefault = -1._RK ! flag to indicate the presence of `pointLogVolNormed`.
            end if

            if (pointLogVolNormedDefault < 0._RK) then
                pointLogVolNormedDefault = pointLogVolNormed
           !elseif (any(Kmeans%Size(1:nc) > nd)) then
            elseif (maxLogVolNormed /= NEGINF_RK) then
                pointLogVolNormedDefault = 0._RK
                do ic = 1, nc ! add all points.
                    if (Kmeans%Size(ic) > nd) pointLogVolNormedDefault = pointLogVolNormedDefault + real(Kmeans%Size(ic),RK)
                end do
                pointLogVolNormedDefault = Kmeans%Prop%logSumVolNormed - log(pointLogVolNormedDefault)
            else
                Kmeans%Err%msg = "There is not any non-singular cluster for which to compute the properties."
                Kmeans%Err%occurred = .true.
                return
            end if

            if (maxLogVolNormed /= NEGINF_RK) then
                Kmeans%Prop%logSumVolNormed = exp(Kmeans%Prop%logSumVolNormed - maxLogVolNormed)
            else
                Kmeans%Prop%logSumVolNormed = 0._RK
                maxLogVolNormed = pointLogVolNormedDefault
            end if

            ! Assume a hyper-spherical ellipsoid as the bounding region of the points.

            loopComputeSingularClusterProperties: do ic = 1, nc

                if (Kmeans%Size(ic) <= nd .and. Kmeans%Size(ic) > 0_IK) then

                    ipstart = Kmeans%Prop%CumSumSize(ic-1) + 1_IK
                    ipend = Kmeans%Prop%CumSumSize(ic)

                    ! Correct the cluster memberships.

                    Kmeans%Membership(ipstart:ipend) = ic

                    ! Compute the scale factor and other properties.

                    do concurrent(ip = 1:np)
                        Kmeans%Prop%MahalSq(ip,ic) = sum( ( Point(1:nd,ip) - Kmeans%Center(1:nd,ic) )**2 )
                    end do

                    Kmeans%Prop%ScaleFactorSq(ic) = max(exp(2*pointLogVolNormedDefault/nd), maxval(Kmeans%Prop%MahalSq(ipstart:ipend,ic)))
                    scaleFactorSqInverse = 1._RK / Kmeans%Prop%ScaleFactorSq(ic)
                    scaleFactor = sqrt(Kmeans%Prop%ScaleFactorSq(ic))

                    !Kmeans%Prop%MahalSq(ip,ic) = Kmeans%Prop%MahalSq(ip,ic) * scaleFactorSqInverse

                    Kmeans%Prop%ChoDia(1:nd,ic) = scaleFactor
                    Kmeans%Prop%LogVolNormed(ic) = nd * log(scaleFactor)
                   !Kmeans%Prop%LogDenNormed(ic) = log(real(Kmeans%Size(ic),RK)) - Kmeans%Prop%LogVolNormed(ic)
                    Kmeans%Prop%ChoLowCovUpp(1:nd,1:nd,ic) = getEye(nd, nd, Kmeans%Prop%ScaleFactorSq(ic))
                    Kmeans%Prop%InvCovMat(1:nd,1:nd,ic) = getEye(nd, nd, scaleFactorSqInverse)

                    ! Compute the effective fraction of points inside the bounded region.

                    if (present(inclusionFraction)) then
                        Kmeans%Prop%EffectiveSize(ic)   = Kmeans%Size(ic) & ! LCOV_EXCL_LINE
                                                        + nint(inclusionFraction * &
                                                        ( count(Kmeans%Prop%MahalSq(1:Kmeans%Prop%CumSumSize(ic-1),ic)<=Kmeans%Prop%ScaleFactorSq(ic)) & ! LCOV_EXCL_LINE
                                                        + count(Kmeans%Prop%MahalSq(Kmeans%Prop%CumSumSize(ic)+1:np,ic)<=Kmeans%Prop%ScaleFactorSq(ic)) & ! LCOV_EXCL_LINE
                                                        ), kind=IK)
!write(*,*) 
!write(*,*) "Kmeans%ic, CumSumSize(ic-1), CumSumSize(ic), np", ic, Kmeans%Prop%CumSumSize(ic-1), Kmeans%Prop%CumSumSize(ic), np
!write(*,*) Kmeans%Size(ic)
!write(*,*) count(Kmeans%Prop%MahalSq(1:Kmeans%Prop%CumSumSize(ic-1),ic)<Kmeans%Prop%ScaleFactorSq(ic))
!write(*,*) count(Kmeans%Prop%MahalSq(Kmeans%Prop%CumSumSize(ic)+1:np,ic)<Kmeans%Prop%ScaleFactorSq(ic))
!write(*,*) 
                    else
                        Kmeans%Prop%EffectiveSize(ic)   = Kmeans%Size(ic)
                    end if
                    !write(*,*) "Kmeans%Prop%LogVolNormed(ic), maxLogVolNormed", Kmeans%Prop%LogVolNormed(ic), maxLogVolNormed
                    Kmeans%Prop%logSumVolNormed = Kmeans%Prop%logSumVolNormed + exp(Kmeans%Prop%LogVolNormed(ic) - maxLogVolNormed)

                end if

            end do loopComputeSingularClusterProperties

            Kmeans%Prop%logSumVolNormed = maxLogVolNormed + log(Kmeans%Prop%logSumVolNormed)

        end if

    end subroutine getProp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    !> \brief
!    !> Remove zero-sized clusters.
!    !> This procedure is a method of the class [Kmeans_type](@ref kmeans_type).
!    !>
!    !> \param[inout]    Kmeans  :   An object of class [Kmeans_type](@ref kmeans_type).
!    !>
!    !> \author
!    ! Amir Shahmoradi, April 03, 2017, 2:16 PM, ICES, UTEXAS
!    pure subroutine packKmeans(Kmeans, minSize)
!#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
!        !DEC$ ATTRIBUTES DLLEXPORT :: packKmeans
!#endif
!        use Constants_mod, only: IK, RK
!        implicit none
!        class(Kmeans_type)  , intent(inout) :: Kmeans
!        integer(IK)                         :: ic, nc
!        logical                             :: Mask(size(Kmeans%Size))
!        Mask = Kmeans%Size < minSize
!        Kmeans%Size = pack(Kmeans%Size, Mask)
!        Kmeans%Center = pack(Kmeans%Center, Mask)
!        Kmeans%Prop%ChoDia = pack(Kmeans%Prop%ChoDia, Mask)
!        real(RK)                    :: logSumVolNormed      !< The logarithm of the sum of the bounding volumes of the clusters *minus their overlaps*.
!       !real(RK)                    :: logAvgDenNormed      !< The logarithm of the average of the densities of points in the bounding volumes of the clusters.
!        real(RK)    , allocatable   :: ChoDia(:,:)          !< An array of size `(nd,nc)` representing the diagonal elements of the Cholesky factorization of the covariance matrix.
!        real(RK)    , allocatable   :: MahalSq(:,:)         !< An array of size `(np,nc)` each element of which represents the Mahalanobis distance squared of point `ip` from the cluster `ic`.
!        real(RK)    , allocatable   :: InvCovMat(:,:,:)     !< An array of size `(nd,nd,nc)` containing the inverse of the covariance matrix of the corresponding cluster.
!        real(RK)    , allocatable   :: LogVolNormed(:)      !< An array of size `(nc)` each element of which represents the log(volume) of the bounding covariance matrix of the corresponding cluster.
!        real(RK)    , allocatable   :: LogDenNormed(:)      !< An array of size `(nc)` each element of which represents the log effective (mean) density of points in the corresponding cluster bounding volume.
!        real(RK)    , allocatable   :: ChoLowCovUpp(:,:,:)  !< An array of size `(nd,nd,nc)` whose upper triangle and diagonal is the covariance matrix and the lower is the Cholesky Lower.
!        real(RK)    , allocatable   :: ScaleFactorSq(:)     !< An array of size `(nc)` representing the factors by which the cluster covariance matrices must be enlarged to enclose their corresponding members.
!        integer(IK) , allocatable   :: EffectiveSize(:)     !< An array of size `(nc)` representing the factors by which the cluster covariance matrices must be enlarged to enclose their corresponding members.
!        integer(IK) , allocatable   :: CumSumSize(:)        !< A vector of size `0:nc` representing the cumulative sum of all cluster sizes from cluster 1 to the last.
!
!    end subroutine packKmeans

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

        if (allocated(Kmeans%Prop%LogVolNormed)) then
            write(fileUnit,"(A)") "LogVolNormed"
            write(fileUnit,fileFormat) Kmeans%Prop%LogVolNormed
        end if

        if (allocated(Kmeans%Prop%ScaleFactorSq)) then
            write(fileUnit,"(A)") "ScaleFactorSq"
            write(fileUnit,fileFormat) Kmeans%Prop%ScaleFactorSq
        end if

        if (allocated(Kmeans%Prop%ChoLowCovUpp)) then
            write(fileUnit,"(A)") "CholeskyLower"
            write(fileUnit,fileFormat) ((Kmeans%Prop%ChoDia(j,ic), (Kmeans%Prop%ChoLowCovUpp(i,j,ic), i=j+1,nd), j=1,nd), ic=1,nc)
        end if

        write(fileUnit,"(A)") "Point"
        write(fileUnit,fileFormat) Point

        if (allocated(Kmeans%Prop%MahalSq)) then
            write(fileUnit,"(A)") "MahalSq"
            write(fileUnit,fileFormat) Kmeans%Prop%MahalSq
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