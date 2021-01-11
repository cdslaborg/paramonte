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

module KmeansOOP_mod

    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type

    implicit none

    character(*), parameter :: MODULE_NAME = "@KmeansOOP_mod"

    real(RK) :: ndKuniform, btFactor

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> The [KmeansTry_type](@ref kmeanstry_type) type containing information about one instance of Kmeans clustering.
    type :: KmeansTry_type
        integer(IK) , allocatable   :: nzsci                !< The number of Zero-Sized Clusters Iterations at each Kmeans try.
        integer(IK) , allocatable   :: niter                !< The number of iterations to reach convergence at each Kmeans try.
        integer(IK) , allocatable   :: Size(:)              !< An array of size `(ncMax)` containing the sizes of the individual clusters identified.
        real(RK)    , allocatable   :: potential            !< The sum of the distances-squared of all points from their corresponding group centers.
        real(RK)    , allocatable   :: Center(:,:)          !< An array of size `(ndMax,ncMax)` containing the centers of the individual clusters identified.
        real(RK)    , allocatable   :: NormedPoint(:,:,:)   !< An array of size `(ndMax,npMax,ncMax)` containing the normalized `Point(nd,np)` w.r.t. each cluster.
        real(RK)    , allocatable   :: MinDistanceSq(:)     !< An array of size `(npMax)` containing the minimum distances of the points from the cluster centers.
        integer(IK) , allocatable   :: Membership(:)        !< An array of size `(npMax)` each element of which represents the cluster ID to which the point belongs.
        type(Err_type)              :: Err                  !< An object of class [Err_type](@ref err_mod::err_type) containing information about error occurrence.
    end type KmeansTry_type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> The `KmeansOOP_type` class.
    type :: KmeansOOP_type
        logical                             :: failed               !< A logical value. If `.true.`, then none of the Kmeans clustering tries have succeeded.
        real(RK)                            :: relTol = 1.e-4_RK    !< The current requested relative tolerance below which should the potential fall, convergence is assumed.
        integer(IK)                         :: nd = 0_IK            !< The current number of data attributes to be used in clustering. The condition `nd<=ndMax` must hold at all times.
        integer(IK)                         :: np = 0_IK            !< The current number of input data points to be clustered. The condition `np<=npMax` must hold at all times.
        integer(IK)                         :: nc = 0_IK            !< The current number of Kmeans clusters to be requested. The condition `nc<=ncMax` must hold at all times.
        integer(IK)                         :: nt = 0_IK            !< The current number of Kmeans clustering trials. The condition `nt<=ntMax` must hold at all times.
        integer(IK)                         :: itopt = 0_IK         !< The Kmeans clustering try which has yielded the least potential.
        integer(IK)                         :: ndMax = 0_IK         !< The maximum possible number of data attributes to be used in clustering. Used for array allocation upon object instantiation.
        integer(IK)                         :: npMax = 0_IK         !< The maximum possible number of input data points to be clustered. Used for array allocation upon object instantiation.
        integer(IK)                         :: ncMax = 0_IK         !< The maximum possible number of Kmeans clusters to be requested. Used for array allocation upon object instantiation.
        integer(IK)                         :: ntMax = 0_IK         !< The maximum possible number of Kmeans clustering trials. Used for array allocation upon object instantiation.
        integer(IK)                         :: niterMax = 1000_IK   !< The maximum number of iterations beyond which it a lack of converge is assumed.
        integer(IK)                         :: nzsciMax = 10_IK     !< The maximum number of zero-sized cluster iterations beyond which it a lack of converge is assumed.
        type(KmeansTry_type), allocatable   :: Try(:)               !< An array of size `(ntMax)` of object of type [KmeansTry_type](@ref kmeanstry_type).
    contains
        procedure, pass :: run => runKmeans
        procedure, pass, private :: allocate => allocateKmeansComponents
    end type KmeansOOP_type

    interface KmeansOOP_type
        module procedure :: constructKmeans
    end interface KmeansOOP_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Construct and return an object of [KmeansOOP_type](@ref KmeansOOP_type) class.
    !> Note that this constructor solely allocates the components of the output
    !> `Kmeans` object to be later used for Kmeans clustering.
    !>
    !> \param[in]       ndMax       :   The maximum possible number of dimensions (features) of the input data array `Point`.
    !> \param[in]       npMax       :   The maximum possible number of observations in the input data array `Point`.
    !> \param[in]       ncMax       :   The maximum possible number of clusters to be found via the Kmeans algorithm.
    !> \param[in]       ntMax       :   The maximum possible number of observations in the input data array `Point`.
    !> \param[in]       niterMax    :   The maximum number of iterations beyond which it a lack of converge is assumed (**optional**, default = 1000).
    !> \param[in]       nszciMax    :   The maximum number of zero-sized cluster iterations beyond which it a lack of converge is assumed (**optional**, default = 10).
    !> \param[in]       relTol      :   The current requested relative tolerance below which should the potential fall, convergence is assumed (**optional**, default = 1.e-8).
    !>
    !> \return
    !> `Kmeans` :   An object of type [KmeansCluster_type](@ref kmeanscluster_type) containing the properties of the identified clusters.
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
    function constructKmeans(ndMax, npMax, ncMax, ntMax, niterMax, nzsciMax, relTol) result(Kmeans)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructKmeans
#endif
        use Constants_mod, only: IK, RK
        implicit none
        type(KmeansOOP_type)                :: Kmeans
        integer(IK) , intent(in)            :: ndMax, npMax, ncMax, ntMax
        integer(IK) , intent(in), optional  :: niterMax, nzsciMax
        real(RK)    , intent(in), optional  :: relTol
        integer(IK)                         :: it
        if (present(relTol)) Kmeans%relTol = relTol
        if (present(niterMax)) Kmeans%niterMax = niterMax
        if (present(nzsciMax)) Kmeans%nzsciMax = nzsciMax
        call Kmeans%allocate(ndMax, npMax, ncMax, ntMax)
    end function constructKmeans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine allocateKmeansComponents(Kmeans, ndMax, npMax, ncMax, ntMax)
        use Constants_mod, only: IK, RK
        implicit none
        class(KmeansOOP_type), intent(inout)    :: Kmeans
        integer(IK) , intent(in)                :: ndMax, npMax, ncMax, ntMax
        integer(IK)                             :: it
        Kmeans%nd = ndMax
        Kmeans%np = npMax
        Kmeans%nc = ncMax
        Kmeans%nt = ntMax
        Kmeans%ndMax = ndMax
        Kmeans%npMax = npMax
        Kmeans%ncMax = ncMax
        Kmeans%ntMax = ntMax
        if (allocated(Kmeans%Try)) deallocate(Kmeans%Try)
        allocate(Kmeans%Try(ntMax))
        do it = 1, Kmeans%ntMax
            allocate(Kmeans%Try(it)%Size(ncMax), source = 0_IK)
            allocate(Kmeans%Try(it)%Membership(npMax), source = 0_IK)
            allocate(Kmeans%Try(it)%Center(ndMax,ncMax), source = 0._RK)
            allocate(Kmeans%Try(it)%MinDistanceSq(npMax), source = 0._RK)
            allocate(Kmeans%Try(it)%NormedPoint(ndMax,npMax,ncMax), source = 0._RK)
        end do
    end subroutine allocateKmeansComponents

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a method of the class [KmeansOOP_type](@ref KmeansOOP_type).
    !> Perform the Kmeans clustering on the input data set represented by the array `Point(nd, np)`.
    !> Upon performing all iterations of clustering, the index `Kmeans%itopt` will point to the final optimal
    !> iteration of Kmeans clustering which resulted in the smallest potential. If only one clustering try has been requested,
    !> then obviously, `Kmeans%itopt = 1`.
    !>
    !> \param[inout]    Kmeans                  :   An object of class [KmeansOOP_type](@ref KmeansOOP_type) which will contain the clustering results.
    !> \param[in]       Point                   :   The input array of size `(nd,np)` of the points to be clustered via the Kmeans algorithm.
    !> \param[in]       nd                      :   The current number of data attributes to be used in clustering. The condition `nd<=ndMax` must hold at all times.
    !> \param[in]       np                      :   The current number of input data points to be clustered. The condition `np<=npMax` must hold at all times.
    !> \param[in]       nc                      :   The current number of Kmeans clusters to be requested. The condition `nc<=ncMax` must hold at all times.
    !> \param[in]       nt                      :   The current number of Kmeans clustering trials. The condition `nt<=ntMax` must hold at all times.
    !> \param[in]       niterMax                :   The maximum number of iterations beyond which it a lack of converge is assumed (**optional**).
    !> \param[in]       nszciMax                :   The maximum number of zero-sized cluster iterations beyond which it a lack of converge is assumed (**optional**).
    !> \param[in]       relTol                  :   The current requested relative tolerance below which should the potential fall, convergence is assumed (**optional**).
    !> \param[in]       InitCenter              :   An input array of size `(nd,nc,nt)` representing the initial centers of the Kmeans clusters (**optional**).
    !>                                              If not provided, the cluster centers will be initialized via the Kmeans++ algorithm.
    !>
    !> \return
    !> `Kmeans` :   An object of type [KmeansCluster_type](@ref kmeanscluster_type) containing the properties of the identified clusters.
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
    subroutine runKmeans(Kmeans, Point, nd, np, nc, nt, niterMax, nzsciMax, relTol, InitCenter)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: runKmeans
#endif
        use Constants_mod, only: IK, RK, NEGINF_IK

        implicit none

        real(RK)    , intent(in)            :: Point(nd,np)
        integer(IK) , intent(in)            :: nd, np, nc, nt
        integer(IK) , intent(in), optional  :: niterMax, nzsciMax
        real(RK)    , intent(in), optional  :: relTol, InitCenter(nd,nc,nt)
        class(KmeansOOP_type)               :: Kmeans

        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME//"@runKmeans()"

        real(RK)                :: potential                ! The sum of all minimum distances squared.
        real(RK)                :: distanceSq               ! distance of each point ip from a given cluster center.
        real(RK)                :: minPotential             ! The minimum of sum of all minimum distances squared across all Kmeans tries.
        real(RK)                :: SumPoint(nd,nc)          ! The sum of the points in each cluster.
        real(RK)                :: MinDistanceSq(np)        ! minimum distance of each point ip from all centers of the clusters.
        integer(IK)             :: Membership(np)           ! current (old) Memberships of the clusters.

        logical                 :: isPresentInitCenter
        logical                 :: convergenceOccurred
        integer                 :: i, j, ip, ic, it

        Kmeans%itopt = 0_IK
        Kmeans%failed = .true.

        if (present(relTol)) Kmeans%relTol = relTol
        if (present(niterMax)) Kmeans%niterMax = niterMax
        if (present(nzsciMax)) Kmeans%nzsciMax = nzsciMax

        if (nd > Kmeans%ndMax .or. np > Kmeans%npMax .or. nc > Kmeans%ncMax .or. nt > Kmeans%ntMax) then
            call Kmeans%allocate(nd, np, nc, nt)
        else
            Kmeans%nd = nd
            Kmeans%np = np
            Kmeans%nc = nc
            Kmeans%nt = nt
        end if

        minPotential = huge(minPotential)
        isPresentInitCenter = present(InitCenter)

        loopMultipleKmeansTry: do it = 1, Kmeans%nt

            Kmeans%Try(it)%Err%occurred = .false.
            Kmeans%Try(it)%Err%stat = NEGINF_IK

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Find the distances of all points to all cluster centers. If not given, initialize the centers via the Kmeans++ algorithm.
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if(isPresentInitCenter) then
                !%%%%%%%%%%%%%%%%%%%
                ic = 1
                do ip = 1, np
                    MinDistanceSq(ip) = sum( (InitCenter(1:nd,ic,it) - Point(1:nd,ip))**2 )
                    Membership(ip) = ic
                end do
                !%%%%%%%%%%%%%%%%%%%
                do ic = 2, nc - 1
                    do ip = 1, np
                        distanceSq = sum( (InitCenter(1:nd,ic,it) - Point(1:nd,ip))**2 )
                        if (distanceSq < MinDistanceSq(ip)) then
                            MinDistanceSq(ip) = distanceSq
                            Membership(ip) = ic
                        end if
                    end do
                end do
                !%%%%%%%%%%%%%%%%%%%
                ic = nc
                SumPoint = 0._RK
                potential = 0._RK
                do ip = 1, np
                    distanceSq = sum( (InitCenter(1:nd,ic,it) - Point(1:nd,ip))**2 )
                    if (distanceSq < MinDistanceSq(ip)) then
                        MinDistanceSq(ip) = distanceSq
                        Membership(ip) = ic
                    end if
                    SumPoint(1:nd,Membership(ip)) = SumPoint(1:nd,Membership(ip)) + Point(1:nd,ip)
                    Kmeans%Try(it)%Size(Membership(ip)) = Kmeans%Try(it)%Size(Membership(ip)) + 1_IK
                    potential = potential + MinDistanceSq(ip)
                end do
                !%%%%%%%%%%%%%%%%%%%
            else
                call runKPP(nc, nd, np, Point, SumPoint, Membership, Kmeans%Try(it)%Size, potential)
            end if

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! Refine all centers until convergence.
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            Kmeans%Try(it)%nzsci = 0_IK
            Kmeans%Try(it)%niter = 1_IK

            loopCenterRefinement: do

                ! compute the new cluster centers and create the new memberships.

                !%%%%%%%%%%%%%%%%%%%
                ic = 1
                Kmeans%Try(it)%Center(1:nd,ic) = SumPoint(1:nd,ic) / real(Kmeans%Try(it)%Size(ic), kind = RK)
                do ip = 1, np
                    Kmeans%Try(it)%NormedPoint(1:nd,ip,ic) = Kmeans%Try(it)%Center(1:nd,ic) - Point(1:nd,ip)
                    Kmeans%Try(it)%MinDistanceSq(ip) = sum( Kmeans%Try(it)%NormedPoint(1:nd,ip,ic)**2 )
                    Kmeans%Try(it)%Membership(ip) = ic
                end do
                !%%%%%%%%%%%%%%%%%%%
                do ic = 2, nc - 1
                    Kmeans%Try(it)%Center(1:nd,ic) = SumPoint(1:nd,ic) / real(Kmeans%Try(it)%Size(ic), kind = RK)
                    do ip = 1, np
                        Kmeans%Try(it)%NormedPoint(1:nd,ip,ic) = Kmeans%Try(it)%Center(1:nd,ic) - Point(1:nd,ip)
                        distanceSq = sum( Kmeans%Try(it)%NormedPoint(1:nd,ip,ic)**2 )
                        if (distanceSq < Kmeans%Try(it)%MinDistanceSq(ip)) then
                            Kmeans%Try(it)%MinDistanceSq(ip) = distanceSq
                            Kmeans%Try(it)%Membership(ip) = ic
                        end if
                    end do
                end do
                !%%%%%%%%%%%%%%%%%%%
                ic = nc
                convergenceOccurred = .true.
                Kmeans%Try(it)%potential = 0._RK
                Kmeans%Try(it)%Center(1:nd,ic) = SumPoint(1:nd,ic) / real(Kmeans%Try(it)%Size(ic), kind = RK)
                do ip = 1, np
                    Kmeans%Try(it)%NormedPoint(1:nd,ip,ic) = Kmeans%Try(it)%Center(1:nd,ic) - Point(1:nd,ip)
                    distanceSq = sum( Kmeans%Try(it)%NormedPoint(1:nd,ip,ic)**2 )
                    if (distanceSq < Kmeans%Try(it)%MinDistanceSq(ip)) then
                        Kmeans%Try(it)%MinDistanceSq(ip) = distanceSq
                        Kmeans%Try(it)%Membership(ip) = ic
                    end if
                    if (Kmeans%Try(it)%Membership(ip)/=Membership(ip)) then
                        convergenceOccurred = .false.
                        Kmeans%Try(it)%Size(Membership(ip)) = Kmeans%Try(it)%Size(Membership(ip)) - 1_IK
                        Kmeans%Try(it)%Size(Kmeans%Try(it)%Membership(ip)) = Kmeans%Try(it)%Size(Kmeans%Try(it)%Membership(ip)) + 1_IK
                    end if
                    Kmeans%Try(it)%potential = Kmeans%Try(it)%potential + Kmeans%Try(it)%MinDistanceSq(ip)
                end do
                !%%%%%%%%%%%%%%%%%%%

                if (any(Kmeans%Try(it)%Size==0_IK)) then

                    Kmeans%Try(it)%niter = 1_IK
                    Kmeans%Try(it)%nzsci = Kmeans%Try(it)%nzsci + 1_IK

                    if (Kmeans%Try(it)%nzsci > Kmeans%nzsciMax) then
                        !Kmeans%Try(it)%Err%msg = "Kmeans%Try(it)%nzsci > Kmeans%nzsciMax"
                        Kmeans%Try(it)%Err%occurred = .true.
                        Kmeans%Try(it)%Err%stat = 2_IK
                        exit loopCenterRefinement
                    end if

                    ! find new random Centers and start over.

                    call runKPP(nc, nd, np, Point, SumPoint, Membership, Kmeans%Try(it)%Size, potential)

                    cycle loopCenterRefinement

#if defined DEBUG_ENABLED || defined TESTING_ENABLED
                elseif (any(Kmeans%Try(it)%Size<0_IK)) then
                    block
                        use Err_mod, only: abort
                        Kmeans%Try(it)%Err%msg = PROCEDURE_NAME//": Internal error occurred - Kmeans%Try(it)%Size < 0. Please report this issue along with circumstances to the developers."
                        call abort(Kmeans%Try(it)%Err)
                        return
                    end block
#endif

                end if

                ! if `niterMax` has reached, quit.

                if (Kmeans%Try(it)%niter < Kmeans%niterMax) then

                    ! If convergence has occurred, exit loopCenterRefinement.

                    if (convergenceOccurred .or. abs(Kmeans%Try(it)%potential-potential)/Kmeans%Try(it)%potential < Kmeans%relTol) exit loopCenterRefinement

                    Kmeans%Try(it)%niter = Kmeans%Try(it)%niter + 1_IK

                    do ip = 1, np
                        SumPoint(1:nd,Membership(ip)) = SumPoint(1:nd,Membership(ip)) - Point(1:nd,ip)
                        SumPoint(1:nd,Kmeans%Try(it)%Membership(ip)) = SumPoint(1:nd,Kmeans%Try(it)%Membership(ip)) + Point(1:nd,ip)
                        Membership(ip) = Kmeans%Try(it)%Membership(ip)
                    end do
                    potential = Kmeans%Try(it)%potential

                else

                    !Kmeans%Try(it)%Err%msg = "Kmeans%Try(it)%niter > Kmeans%niterMax"
                    Kmeans%Try(it)%Err%occurred = .true.
                    Kmeans%Try(it)%Err%stat = 1_IK
                    exit loopCenterRefinement

                end if

            end do loopCenterRefinement

            if (.not. Kmeans%Try(it)%Err%occurred) then
                if (Kmeans%Try(it)%potential < minPotential) Kmeans%itopt = it
                Kmeans%failed = .false.
            end if

        end do loopMultipleKmeansTry

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
        integer(IK)                 :: randomInteger
        integer(IK)                 :: ic, ip, i

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

end module KmeansOOP_mod ! LCOV_EXCL_LINE