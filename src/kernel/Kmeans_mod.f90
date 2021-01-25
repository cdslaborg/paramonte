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
    !> The inclusion of the component `NormedPoint` adds ~50% to the computational cost of performing Kmeans
    type :: Kmeans_type
        real(RK)                    :: potential            !< The sum of the distances-squared of all points from their corresponding group centers.
        integer(IK)                 :: nzsci                !< The number of Zero-Sized Clusters Iterations.
        integer(IK)                 :: niter                !< The number of iterations to reach convergence.
        integer(IK) , allocatable   :: Size(:)              !< A vector of size `nc` (the number of clusters) containing the sizes of the individual clusters identified.
        real(RK)    , allocatable   :: Center(:,:)          !< An array of size `(nd,nc)` containing the centers of the individual clusters identified.
                                                            !< Here `nd` represents the number of attributes in the data points that have been clustered.
        real(RK)    , allocatable   :: NormedPoint(:,:,:)   !< An array of size `(nd,np,nc)` containing `nc` copies of the input array `Point(nd,np)`
                                                            !< each copy of which is normalized with respect to the corresponding cluster center.
        real(RK)    , allocatable   :: MinDistanceSq(:)     !< An array of size `(np)` containing the minimum distances of the points from the cluster centers.
                                                            !< Here `np` represents the number of data points to be clustered.
        integer(IK) , allocatable   :: Membership(:)        !< A vector of size `np` each element of which represents the cluster ID to which the point belongs.
                                                            !< Here `np` represents the number of data points to be clustered.
        type(Err_type)              :: Err                  !< An object of class [Err_type](@ref err_mod::err_type) containing information about error occurrence.
    end type Kmeans_type

    interface Kmeans_type
        module procedure :: getKmeans, allocateKmeans, optimizeKmeans
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
                , Kmeans%NormedPoint(nd,np,nc) & ! LCOV_EXCL_LINE
                , Kmeans%MinDistanceSq(np) & ! LCOV_EXCL_LINE
                , Kmeans%Membership(np) & ! LCOV_EXCL_LINE
                )
    end function allocateKmeans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Perform the Kmeans clustering on the input data set represented by the array `Point(nd, np)`
    !> and return the resulting cluster characteristics as an object of type [Kmeans_type](@ref kmeans_type).
    !>
    !> \param[in]       nd          :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[in]       np          :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[in]       nc          :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[in]       Point       :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[in]       InitCenter  :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[in]       niterMax    :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[in]       nszciMax    :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[in]       relTol      :   See the description of the [runKmeans](@ref runkmeans).
    !>
    !> \return
    !> `Kmeans` :   An object of type [Kmeans_type](@ref Kmeans_type) containing the properties of the identified clusters.
    !>              The `Err%stat` component of the output `Kmeans` will be set to `1` if the number of iterations reaches the input
    !>              `niterMax` or its default value before convergence occurs.
    !>              The `Err%stat` component of the output `Kmeans` will be set to `2` if more than `nzsciMax` clustering tries
    !>              result in clusters with zero members.
    !>              In both cases in the above, `Err%occurred = .true.` upon exiting the procedure.
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
    function getKmeans(nd, np, nc, Point, InitCenter, niterMax, nzsciMax, relTol) result(Kmeans)
        implicit none
        integer(IK) , intent(in)            :: nd, np, nc
        real(RK)    , intent(in)            :: Point(nd,np)
        real(RK)    , intent(in), optional  :: InitCenter(nd,nc)
        integer(IK) , intent(in), optional  :: niterMax, nzsciMax
        real(RK)    , intent(in), optional  :: relTol
        type(Kmeans_type)                   :: Kmeans
        Kmeans = allocateKmeans(nd, np, nc)
        call runKmeans  ( nd = nd & ! LCOV_EXCL_LINE
                        , np = np & ! LCOV_EXCL_LINE
                        , nc = nc & ! LCOV_EXCL_LINE
                        , Point = Point & ! LCOV_EXCL_LINE
                        , Size = Kmeans%Size & ! LCOV_EXCL_LINE
                        , Center = Kmeans%Center & ! LCOV_EXCL_LINE
                        , Membership = Kmeans%Membership & ! LCOV_EXCL_LINE
                        , NormedPoint = Kmeans%NormedPoint & ! LCOV_EXCL_LINE
                        , MinDistanceSq = Kmeans%MinDistanceSq & ! LCOV_EXCL_LINE
                        , potential = Kmeans%potential & ! LCOV_EXCL_LINE
                        , niter = Kmeans%niter & ! LCOV_EXCL_LINE
                        , nzsci = Kmeans%nzsci & ! LCOV_EXCL_LINE
                        , Err = Kmeans%Err & ! LCOV_EXCL_LINE
                        , relTol = relTol & ! LCOV_EXCL_LINE
                        , niterMax = niterMax & ! LCOV_EXCL_LINE
                        , nzsciMax = nzsciMax & ! LCOV_EXCL_LINE
                        , InitCenter = InitCenter & ! LCOV_EXCL_LINE
                        )
    end function getKmeans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Perform the Kmeans clustering for `nt` tries on the input data set represented by the array `Point(nd, np)`
    !> and return the clustering that yields the least potential.
    !>
    !> \param[in]       nd          :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[in]       np          :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[in]       nc          :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[in]       nt          :   The number of independent instances of Kmeans to run.
    !> \param[in]       Point       :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[in]       InitCenter  :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[in]       niterMax    :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[in]       nszciMax    :   See the description of the [runKmeans](@ref runkmeans).
    !> \param[in]       relTol      :   See the description of the [runKmeans](@ref runkmeans).
    !>
    !> \return
    !> `Kmeans` :   An object of type [Kmeans_type](@ref Kmeans_type) containing the properties of the identified clusters.
    !>              The `Err%stat` component of the output `Kmeans` will be set to `1` if the number of iterations reaches the input
    !>              `niterMax` or its default value before convergence occurs.
    !>              The `Err%stat` component of the output `Kmeans` will be set to `2` if more than `nzsciMax` clustering tries
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
    function optimizeKmeans(nd, np, nc, nt, Point, InitCenter, niterMax, nzsciMax, relTol) result(Kmeans)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: optimizeKmeans
#endif
        use Constants_mod, only: IK, RK, HUGE_RK

        implicit none

        integer(IK) , intent(in)            :: nt, nc, nd, np
        real(RK)    , intent(in)            :: Point(nd,np)
        real(RK)    , intent(in), optional  :: InitCenter(nd,nc)
        integer(IK) , intent(in), optional  :: niterMax, nzsciMax
        real(RK)    , intent(in), optional  :: relTol
        type(Kmeans_type), allocatable      :: KmeansTry(:)
        type(Kmeans_type)                   :: Kmeans
        real(RK)                            :: potential
        integer(IK)                         :: itry, itryLeastPotential

        if (allocated(KmeansTry)) deallocate(KmeansTry) ! LCOV_EXCL_LINE
        allocate(KmeansTry(nt))
        potential = HUGE_RK

        do itry = 1, nt
            KmeansTry(itry) = getKmeans ( nc = nc & ! LCOV_EXCL_LINE
                                        , nd = nd & ! LCOV_EXCL_LINE
                                        , np = np & ! LCOV_EXCL_LINE
                                        , Point = Point & ! LCOV_EXCL_LINE
                                        , InitCenter = InitCenter & ! LCOV_EXCL_LINE
                                        , niterMax = niterMax & ! LCOV_EXCL_LINE
                                        , nzsciMax = nzsciMax & ! LCOV_EXCL_LINE
                                        , relTol = relTol & ! LCOV_EXCL_LINE
                                        )
            if (KmeansTry(itry)%potential < potential) then
                potential = KmeansTry(itry)%potential
                itryLeastPotential = itry
            end if
        end do

        Kmeans = KmeansTry(itryLeastPotential)

    end function optimizeKmeans

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
    !> \param[out]      NormedPoint     :   An array of size `(nd,np,nc)` containing the normalized input array `Point(nd,np)` w.r.t. each cluster center.
    !> \param[out]      MinDistanceSq   :   An array of size `(np)` containing the minimum distances of the points from the cluster centers.
    !> \param[out]      potential       :   The sum of the distances-squared of all points from their corresponding group centers.
    !> \param[out]      niter           :   The number of iterations to reach convergence.
    !> \param[out]      nzsci           :   The number of Zero-Sized Clusters Iterations.
    !> \param[out]      Err             :   An object of class [Err_type](@ref err_mod::err_type) containing information about error occurrence.
    !> \param[in]       InitCenter      :   An input array of size `(nd,nc)` representing the initial centers of the Kmeans clusters (**optional**).
    !>                                      If not provided, the cluster centers will be initialized via the Kmeans++ algorithm.
    !> \param[in]       niterMax        :   The maximum number of iterations beyond which it a lack of converge is assumed (**optional**, default = 1000).
    !> \param[in]       nszciMax        :   The maximum number of zero-sized cluster iterations beyond which it a lack of converge is assumed (**optional**, default = 10).
    !> \param[in]       relTol          :   If the **relative** difference between two subsequent potentials is below `relTol`, convergence is assumed (**optional**, default = 1.e-4).
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
                        , NormedPoint & ! LCOV_EXCL_LINE
                        , MinDistanceSq & ! LCOV_EXCL_LINE
                        , potential & ! LCOV_EXCL_LINE
                        , niter & ! LCOV_EXCL_LINE
                        , nzsci & ! LCOV_EXCL_LINE
                        , Err & ! LCOV_EXCL_LINE
                        , InitCenter & ! LCOV_EXCL_LINE
                        , niterMax & ! LCOV_EXCL_LINE
                        , nzsciMax & ! LCOV_EXCL_LINE
                        , relTol & ! LCOV_EXCL_LINE
                        )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: runKmeans
#endif
        use Constants_mod, only: IK, RK, HUGE_RK

        implicit none

        integer(IK)     , intent(in)            :: nc, nd, np
        real(RK)        , intent(in)            :: Point(nd,np)
        integer(IK)     , intent(out)           :: nzsci
        integer(IK)     , intent(out)           :: niter
        real(RK)        , intent(out)           :: potential
        real(RK)        , intent(out)           :: Center(nd,nc)
        real(RK)        , intent(out)           :: NormedPoint(nd,np,nc)
        real(RK)        , intent(out)           :: MinDistanceSq(np)
        integer(IK)     , intent(out)           :: Membership(np)
        integer(IK)     , intent(out)           :: Size(nc)
        type(Err_type)  , intent(out)           :: Err
        real(RK)        , intent(in), optional  :: InitCenter(nd,nc)
        integer(IK)     , intent(in), optional  :: niterMax, nzsciMax
        real(RK)        , intent(in), optional  :: relTol

        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME//"@runKmeans()"
        integer(IK) , parameter :: MAX_NITER = 1e3_IK       ! The maximum possible number of iterations for kmeans to converge.
        integer(IK) , parameter :: MAX_NZSCI = 1e1_IK       ! The maximum possible number of iterations for which clusters of zero-size are identified.
        real(RK)    , parameter :: TOL_DEFAULT = 1.e-4_RK   ! if the relative change in all Centers is less than this value, the algorithm has converged.

        real(RK)                :: distanceSq               ! distance of each point ip from a given cluster center.
        real(RK)                :: potentialOld             ! The sum of all minimum distances squared.
        real(RK)                :: SumPoint(nd,nc)          ! The sum of the points in each cluster.
        integer(IK)             :: MembershipOld(np)        ! current (old) Memberships of the clusters.

        logical                 :: convergenceOccurred
        integer(IK)             :: niterMaxDefault
        integer(IK)             :: nzsciMaxDefault
        real(RK)                :: relTolDefault
        integer                 :: ip, ic

        Err%occurred = .false.

        relTolDefault = TOL_DEFAULT; if (present(relTol)) relTolDefault = relTol
        niterMaxDefault = MAX_NITER; if (present(niterMax)) niterMaxDefault = niterMax
        nzsciMaxDefault = MAX_NZSCI; if (present(nzsciMax)) nzsciMaxDefault = nzsciMax

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

        nzsci = 0_IK
        niter = 1_IK

        loopCenterRefinement: do

            ! compute the new cluster centers and create the new memberships.

            !%%%%%%%%%%%%%%%%%%%
            ic = 1
            Center(1:nd,ic) = SumPoint(1:nd,ic) / real(Size(ic), kind = RK)
            do ip = 1, np
                NormedPoint(1:nd,ip,ic) = Center(1:nd,ic) - Point(1:nd,ip)
                MinDistanceSq(ip) = sum( NormedPoint(1:nd,ip,ic)**2 )
                Membership(ip) = ic
            end do
            !%%%%%%%%%%%%%%%%%%%
            do ic = 2, nc - 1
                Center(1:nd,ic) = SumPoint(1:nd,ic) / real(Size(ic), kind = RK)
                do ip = 1, np
                    NormedPoint(1:nd,ip,ic) = Center(1:nd,ic) - Point(1:nd,ip)
                    distanceSq = sum( NormedPoint(1:nd,ip,ic)**2 )
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
                NormedPoint(1:nd,ip,ic) = Center(1:nd,ic) - Point(1:nd,ip)
                distanceSq = sum( NormedPoint(1:nd,ip,ic)**2 )
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
                nzsci = nzsci + 1_IK

                if (nzsci > nzsciMaxDefault) then
                    !Err%msg = "nzsci > nzsciMaxDefault"
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

end module Kmeans_mod ! LCOV_EXCL_LINE