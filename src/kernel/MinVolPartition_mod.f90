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
!> This module contains procedures and routines for the minimum bounding ellipsoid partitioning of a given set of points.
!> \author
!> Amir Shahmoradi

module MinVolPartition_mod

    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type

    implicit none

    character(*), parameter :: MODULE_NAME = "@MinVolPartition_mod"

    type :: Ellipsoid_type
        real(RK), allocatable :: Center(:,:)                    !< An array of size `(nd,ndiv)` containing the cluster centers.
        real(RK), allocatable :: MahalSq(:,:)                   !< An array of size `(np,ndiv)` containing the maximum Mahalanobis distances of all points with respect to each bounding ellipsoids.
        real(RK), allocatable :: NormData(:,:)                  !< An array of size `(nd,np,ndiv)` containing the cluster members normalized with respect to their corresponding cluster centers.
        real(RK), allocatable :: InvCovMat(:,:,:)               !< An array of size `(nd,nd,ndiv)` containing the full symmetric inverse covariance matrices of the bounding ellipsoids.
        real(RK), allocatable :: CholDiagLower(:,:,:)           !< An array of size `(nd,0:nd,ndiv)` the Cholesky lower triangle, diagonal, and covariance matrices of the bounding ellipsoids.
        real(RK), allocatable :: LogVolNormed(:)                !< An array of size `(ndiv)` containing the log volumes of the bounding ellipsoids.
        real(RK), allocatable :: ScaleFacSq(:)                  !< An array of size `(ndiv)` containing the square of `ScaleFacSq`.
        real(RK), allocatable :: ScaleFac(:)                    !< An array of size `(ndiv)` containing the factors by which the determinants of the the bounding
                                                                !< ellipsoids must be multiplied in order for the ellipsoid to become a bounding ellipsoid.
    end type Ellipsoid_type

    !> The `MinVolPartition_type` class. Partitions an input array `Point(nd,np)` with `nd` attributes and `np` observations (points).
    type :: MinVolPartition_type
        integer(IK)                 :: ndiv                     !< The number of sample division (clusters) at each level of partitioning. A default value would be `2`.
        integer(IK)                 :: neMax                    !< The maximum possible number of clusters.
        integer(IK)                 :: neOptimal                !< The predicted optimal number of clusters.
        integer(IK)                 :: convergenceFailureCount  !< The number of times the algorithm has failed to converge.
        type(Ellipsoid_type)        :: Backup
        integer(IK) , allocatable   :: Membership(:)            !< An array of size `(np)` representing the bounding-ellipsoid membership IDs of the corresponding data points.
        integer(IK) , allocatable   :: PointIndex(:)            !< An array of size `(np)` representing the indices of the original data points before partitioning.
        integer(IK) , allocatable   :: Size(:)                  !< An array of size `(neMax)` representing the sizes of the corresponding bounding ellipsoids.
        integer(IK) , allocatable   :: EllIndex(:)              !< An array of size `(neMax)` containing the indices of the bounding ellipsoids that are terminal and cannot be further divided.
        real(RK)    , allocatable   :: Center(:,:)              !< An array of size `(nd,neMax)` representing the centers of the bounding ellipsoids.
        real(RK)    , allocatable   :: LogVolNormed(:)          !< An array of size `(neMax)` representing the log-volumes of the corresponding bounding ellipsoids.
        real(RK)    , allocatable   :: InvCovMat(:,:,:)         !< An array of size `(nd,nd,neMax)` representing the full symmetric inverse covariance matrices of the bounding ellipsoids.
        real(RK)    , allocatable   :: CholDiagLower(:,:,:)     !< An array of size `(nd,0:nd,neMax)` representing the Cholesky lower triangle, diagonal, and covariance matrices of the bounding ellipsoids.
                                                                !< The zeroth index of the second dimension contains the Cholesky diagonal elements.
        type(Err_type)              :: Err                      !< An object of class [Err_type](@ref err_mod::err_type) containing information about error occurrence.
    contains
        procedure, pass         :: get => getMinVolPartition
    end type MinVolPartition_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Construct and return an object of class [MinVolPartition_type](@ref minvolpartition_type) to be later used for (repeated)
    !> partitioning of a data array of size `Point(nd,np)` with `nd` attributes and `np` observations. This constructor solely
    !> computes the maximum number of possible bounding ellipsoids (`neMax`) and allocates the components of the object.
    !>
    !> \param[in]       nd                      :   The number of dimensions (features) of the input data array `Point`.
    !> \param[in]       npMax                   :   The maximum conceivable number of observations in the input data array `Point` (`np <= npMax`).
    !> \param[in]       stabilizationRequested  :   A logical flag indicating whether stabilization of the volumes of ellipsoids must be performed.
    !> \param[in]       maxAllowedLoopRecursion :   A logical flag representing the maximum allowed number of recursions before convergence occurs.
    !> \param[in]       maxAllowedKmeansFailure :   A logical flag representing the maximum number of times the kmeans clustering returns a cluster with `size <= nd`.
    !> \param[in]       mahalSqWeightExponent   :   The positive real value to be used to assign weights the Mahalanobis distances of the points.
    !> \param[in]       inclusionFraction       :   The fraction of non-cluster-member points that ARE inside the cluster, to be considered in the estimation of the true volume of the cluster.
    !> \param[in]       tightness               :   The factor determining how large the bounding ellipsoid can be compared to the true volume of the ellipsoid.
    !>
    !> \return
    !> `self` : An object of class [MinVolPartition_type](@ref minvolpartition_type) to be used for partitioning of the points.
    !>
    !> \warning
    !> This algorithm does not check for the consistency of the input values for `nd` and `npMax`.
    !> It is currently the user's responsibility to ensure that these input values are logical and sensible.
    !>
    !> \author
    ! Amir Shahmoradi, April 03, 2017, 2:16 PM, ICES, UTEXAS
    function constructMinVolPartition   ( nd, npMax, ndiv           &
                                        , stabilizationRequested    &
                                        , maxAllowedLoopRecursion   &
                                        , maxAllowedKmeansFailure   &
                                        , maxAllowedMinVolFailure   &
                                        , mahalSqWeightExponent     &
                                        , inclusionFraction         &
                                        , tightness                 &
                                        ) result(self)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructMinVolPartition
#endif
        use Constants_mod, only: IK, RK

        implicit none

        integer(IK) , intent(in)    :: nd, npMax, ndiv
        logical     , intent(in)    :: stabilizationRequested
        integer     , intent(in)    :: maxAllowedLoopRecursion
        integer     , intent(in)    :: maxAllowedKmeansFailure
        integer     , intent(in)    :: maxAllowedMinVolFailure
        real(RK)    , intent(in)    :: mahalSqWeightExponent
        real(RK)    , intent(in)    :: inclusionFraction
        real(RK)    , intent(in)    :: tightness
        type(MinVolPartition_type)  :: self

        self%ndiv = ndiv
        self%neMax = npMax / (nd + 1)
        self%Err%occurred = .false.

        allocate(self%Membership(nd,npMax))
        allocate(self%PointIndex(nd,npMax))
        allocate(self%Size(self%neMax))
        allocate(self%Center(nd,self%neMax))
        allocate(self%EllIndex(self%neMax))
        allocate(self%LogVolNormed(self%neMax))
        allocate(self%InvCovMat(nd,nd,self%neMax))
        allocate(self%CholDiagLower(nd,0:nd,self%neMax))

        ! allocate the stabilization backup variables

        allocate(self%Backup%Center(nd,ndiv))
        allocate(self%Backup%MahalSq(npMax,ndiv))
        allocate(self%Backup%InvCovMat(nd,nd,ndiv))
        allocate(self%Backup%NormData(nd,npMax,ndiv))
        allocate(self%Backup%CholDiagLower(nd,0:nd,ndiv))
        allocate(self%Backup%LogVolNormed(ndiv))
        allocate(self%Backup%ScaleFacSq(ndiv))
        allocate(self%Backup%ScaleFac(ndiv))

    end function constructMinVolPartition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This function is a method of the class [MinVolPartition_type](@ref minvolpartition_type).
    !> Perform the minimum-bounding-volume partitioning of the input set of `Point(nd,np)`.
    !> The partitioning result will be stored in the components of `self` upon exit.
    !>
    !> \param[inout]    self        :   An object of class [MinVolPartition_type](@ref minvolpartition_type) that will contain the partitioning information.
    !> \param[in]       nd          :   The number of dimensions (features) of the input data array `Point`.
    !> \param[in]       np          :   The number of observations in the input data array `Point`.
    !> \param[in]       parLogVol   :   The input log-estimate of the total volume of the points (the PARent ellipsoid's LOG-VOLume).
    !> \param[inout]    Point       :   The input array of size `(nd,np)` to be partitioned. On output, it is reordered.
    !>
    !> \warning
    !> The input value `np` must be always smaller than the value passed for `npMax` when constructing the parent object of this method.
    !>
    !> \author
    ! Amir Shahmoradi, April 03, 2017, 2:16 PM, ICES, UTEXAS
    subroutine getMinVolPartition(self, nd, np, parLogVol, Point)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinVolPartition
#endif
        use Statistics_mod , only: getSamCholFac, getMahalSq
        use Matrix_mod, only: getInvMatFromCholFac
        use Kmeans_mod, only: Kmeans_type
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str
        use Err_mod, only: abort

        implicit none

        type(MinVolPartition_type)  , intent(inout) :: self
        integer(IK)                 , intent(in)    :: nd, np
        real(RK)                    , intent(in)    :: parLogVol
        real(RK)                    , intent(inout) :: Point(nd,np)

        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME//"@getMinVolPartition()"

        real(RK)    :: parLogVolNormed      ! Normalized log(volume) of the parent ellipsoid
        real(RK)    :: logNormFac, scaleFac

        integer(IK) :: ipend
        integer(IK) :: ipstart
        integer(IK) :: npCurrent
        integer(IK) :: ndPlusOne
        integer(IK) :: ndEffective
        integer(IK) :: recursionCounter
        integer(IK) :: loopPartitionCounter
        integer(IK) :: kmeansFailureCounter
        integer(IK) :: minVolFailureCounter
        integer(IK) :: i, j, ic, ip

        self%Center(1:nd,1) = 0._RK
        do j = 1, np
            self%PointIndex(j) = j
            self%Center(1:nd,1) = self%Center(1:nd,1) + Point(1:nd,j)
        end do
        self%Center(1:nd,1) = self%Center(1:nd,1) / real(np,kind=RK)

        ! Get the Cholesky Factor of the partitioning sample.
        ! @todo: There is room to improve the efficiency of computing `CholDiagLower` via a recursive update.

        call getSamCholFac  ( nd = nd                                           & ! LCOV_EXCL_LINE
                            , np = np                                           & ! LCOV_EXCL_LINE
                            , Mean = self%Center(1:nd,1)                        & ! LCOV_EXCL_LINE
                            , Point = Point                                     & ! LCOV_EXCL_LINE
                            , CholeskyLower = self%CholDiagLower(1:nd,1:nd,1)   & ! LCOV_EXCL_LINE
                            , CholeskyDiago = self%CholDiagLower(1:nd,0,1)      & ! LCOV_EXCL_LINE
                            )
        if (self%CholDiagLower(1,0,1)<0._RK) then
            self%Err%msg = PROCEDURE_NAME//": Cholesky factorization failed. nd, np, neMax = "//num2str(nd)//", "//num2str(np)//", "//num2str(neMax)
            self%Err%occurred = .true.
            call abort(self%Err)
            return
        end if

        ! Get the LogVolume of the first cluster for normalizing the volumes of children clusters.

        self%InvCovMat(1:nd,1:nd,1) = getInvMatFromCholFac(nd, self%CholDiagLower(1:nd,1:nd,1), self%CholDiagLower(1:nd,0,1))
        scaleFacSq = maxval( getMahalSq(nd, np, self%Center(1:nd,1), self%InvCovMat(1:nd,1:nd,1), Point) )
        scaleFac = sqrt(scaleFacSq)

        self%LogVolNormed(1) = sum(log( self%CholDiagLower(1:nd,0,1) )) + nd * log(scaleFac)

        if (self%LogVolNormed(1) > parLogVol) then
            logNormFac = self%LogVolNormed(1)
            self%LogVolNormed(1) = 1._RK
            parLogVolNormed = exp( parLogVol - logNormFac )
        else ! enlarge the bounding ellipsoid
            !write(*,*) 'no no no no!'
            logNormFac = parLogVol
            scaleFac = scaleFac * exp( ( logNormFac - self%LogVolNormed(1) ) / nd )
            scaleFacSq = scaleFac**2
            self%LogVolNormed(1) = 1._RK
            parLogVolNormed = 1._RK
        end if

        blockPartition: if (self%neMax > 1 .and. self%LogVolNormed(1) > self%tightness * parLogVolNormed) then

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ! begin loopPartition
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ipend = np
            ipstart = 1
            ndPlusOne = nd + 1
            ndEffective = np / neMax + 1

            loopPartition: do

                npCurrent = ipend - ipstart + 1
                loopPartitionCounter = loopPartitionCounter + 1

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! Partition the current sample into ndiv clusters.
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                kmeansFailureCounter = -1
                loopKmeans: do
                    kmeansFailureCounter = kmeansFailureCounter + 1
                    Kmeans = Kmeans_type( nc = self%ndiv
                                        , nd = nd
                                        , np = npCurrent
                                        , Point = Point(ipstart:ipend)
                                        !, InitCenter
                                        !, niterMax
                                        !, nzsciMax
                                        !, reltolSq
                                        )
                    if ( Kmeans%Err%occurred .or. any(Kmeans%Size < ndPlusOne) ) then
                        if (kmeansFailureCounter < self%maxAllowedKmeansFailure) cycle loopKmeans
                        call setSingleTerminalEllipsoid()
                        cycle loopPartition
                    end if
                    exit loopKmeans
                end do loopKmeans

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! Find ndiv approximately minimum-volume sub-cllusters recursively.
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                recursionCounter = 0
                minVolFailureCounter = 0
                self%Backup%LogVolNormed = huge(0._RK)

                loopRecursiveRunMinVolKernel: do

                    recursionCounter = recursionCounter + 1
                    if(recursionCounter > maxAllowedLoopRecursion) then
                        call setSingleTerminalEllipsoid()
                        cycle loopPartition
                    end if

                    ! Compute the clusters' properties including the scaled Mahalanobis distances

                    do ip = ipstart, ipend
                        self%Backup%NormData(1:nd,ip,ic) = Point(1:nd,ip) - Kmeans%Center(ic)
                    end do

                    do ic = 1, self%ndiv

                        ! compute the Cholesky Factor of cluster samples

                        ! Only upper half of CovMat is needed
                        npMinusOneInverse = 1._RK / real(Kmeans%Size(ic)-1,kind=RK)
                        do j = 1,nd
                            do i = 1,j
                                ! Get the covariance matrix elements: only the upper half of CovMat is needed
                                Kmeans%Size(ic) = 
                                CholeskyLower(i,j) = dot_product( NormData(1:Kmeans%Size(ic),i) , NormData(1:np,j) ) * npMinusOneInverse
                            end do
                        end do

                        ! @todo: The efficiency can be improved by merging it with the above loops
                        call getCholeskyFactor(nd, CholeskyLower, CholeskyDiago)

                        call getSamCholFac  ( nd                                  &
                                            , CKSize(ic)                          &
                                            , CKCenter(1:nd,ic)                   &
                                            , Point(1:nd,IpStart(ic):IpEnd(ic))   &
                                            , CKBECholDiagLower(1:nd,1:nd,ic)     &
                                            , CKBECholDiag(1:nd,ic)               &
                                            )
                        if (CKBECholDiag(1,ic)<0._RK) then
                            write(*,*) "    FATAL: runMinVolKernel() @ getSamCholFac() @ getCholeskyFactor() failed."
                            write(*,*) "           Cluster id, dimension, size:"
                            write(*,*) ic, nd, CKSize(ic)
                            stop
                        end if

                        CKBEInvCovMat(1:nd,1:nd,ic) = getInvMatFromCholFac( nd, CKBECholDiagLower(1:nd,1:nd,ic), CKBECholDiag(1:nd,ic) )
                        CKMahalSq(1:np,ic) = getMahalSq( nd, np, CKCenter(1:nd,ic), CKBEInvCovMat(1:nd,1:nd,ic), Point )
                        ScaleFacSq(ic) = maxval(CKMahalSq(IpStart(ic):IpEnd(ic),ic))
                        ScaleFac(ic) = sqrt( ScaleFacSq(ic) )
                        CKMahalSq(1:np,ic) = CKMahalSq(1:np,ic) / ScaleFacSq(ic)
                        CKBELogVolNormed(ic) = exp( nd*log(ScaleFac(ic)) + sum(log(CKBECholDiag(1:nd,ic))) - logNormFac ) ! Get the LogVolume of the cluster
                        NumPointInCK(ic) = CKSize(ic) + inclusionFraction*count(CKMahalSq(IpStart(icOther):IpEnd(icOther),ic)<1._RK)
                        CKParLogVolNormed(ic) = real(NumPointInCK(ic),kind=RK) * parLogVolNormed / real(np,kind=RK)

                        blockStabilityCheck: if (stabilizationRequested) then

                            if ( CKBELogVolNormed(ic) > self%Backup%LogVolNormed(ic) ) then    ! This is to ensure the convergence of the algorithm

                                CKCenter(1:nd,ic) = CKCenterOld(1:nd,ic)
                                CKBECholDiagLower(1:nd,1:nd,ic) = CKBECholDiagLowerOld(1:nd,1:nd,ic)
                                CKBECholDiag(1:nd,ic) = CKBECholDiagOld(1:nd,ic)
                                CKBEInvCovMat(1:nd,1:nd,ic) = CKBEInvCovMatOld(1:nd,1:nd,ic)
                                ScaleFacSq(ic) = ScaleFacSqOld(ic)
                                ScaleFac(ic) = ScaleFacOld(ic)
                                CKBELogVolNormed(ic) = self%Backup%LogVolNormed(ic)

                                CMemberCounter(1) = 0
                                CMemberCounter(2) = CKSize(1)
                                do ip = 1,np
                                    if (CMembership(ip)==1) then
                                        CMemberCounter(1) = CMemberCounter(1) + 1
                                        CKMahalSq(CMemberCounter(1),ic) = CKMahalSqOld(ip,ic)
                                    else
                                        CMemberCounter(2) = CMemberCounter(2) + 1
                                        CKMahalSq(CMemberCounter(2),ic) = CKMahalSqOld(ip,ic)
                                    end if
                                end do
                                CKMahalSqOld(1:np,ic) = CKMahalSq(1:np,ic)

                            else  ! no stability correction is required. Everything just fine...

                                CKCenterOld(1:nd,ic) = CKCenter(1:nd,ic)
                                CKBECholDiagLowerOld(1:nd,1:nd,ic) = CKBECholDiagLower(1:nd,1:nd,ic)
                                CKBECholDiagOld(1:nd,ic) = CKBECholDiag(1:nd,ic)
                                CKBEInvCovMatOld(1:nd,1:nd,ic) = CKBEInvCovMat(1:nd,1:nd,ic)
                                CKMahalSqOld(1:np,ic) = CKMahalSq(1:np,ic)
                                ScaleFacSqOld(ic) = ScaleFacSq(ic)
                                ScaleFacOld(ic) = ScaleFac(ic)
                                self%Backup%LogVolNormed(ic) = CKBELogVolNormed(ic)

                            end if

                        end if blockStabilityCheck

                        !if (CKBELogVolNormed(ic)<CKParLogVolNormed(ic)) then  ! enlarge the bounding ellipsoid
                        !  !write(*,*) 'enlargement happening...'
                        !  scaleFacSq(ic) = scaleFacSq(ic) * (CKParLogVolNormed(ic)/CKBELogVolNormed(ic))**(2._RK/nd)  ! ATTN: (1._RK/nd) this should likely be (2._RK/nd)
                        !  scaleFac(ic) = sqrt(scaleFacSq(ic))
                        !  CKBELogVolNormed(ic) = CKParLogVolNormed(ic)
                        !end if

                        debugVolNorm: block
                            real(RK) :: dummy
                            dummy = real(CKSize(ic),kind=RK)*parLogVolNormed/real(np,kind=RK)
                            if ( CKBELogVolNormed(ic) < dummy ) then  ! enlarge the bounding ellipsoid
                                !write(*,*) 'enlargement happening...'
                                scaleFacSq(ic) = scaleFacSq(ic) * (dummy/CKBELogVolNormed(ic))**(2._RK/nd)  ! ATTN: (1._RK/nd) this should likely be (2._RK/nd)
                                scaleFac(ic) = sqrt(scaleFacSq(ic))
                                CKBELogVolNormed(ic) = dummy
                            end if
                        end block debugVolNorm

                        !end if blockStabilityCheck

                    end do

                    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    ! Reassign points to clusters if needed
                    MahalSqWeight = ( CKBELogVolNormed / CKParLogVolNormed )**mahalSqWeightExponent
                    reclusteringNeeded = .false.
                    CKCenter(1:nd,1) = CKCenter(1:nd,1) * CKSize(1)  ! Now it is sum instead of mean
                    CKCenter(1:nd,2) = CKCenter(1:nd,2) * CKSize(2)  ! Now it is sum instead of mean
                    ic = 1
                    icOther = 2
                    do ip = IpStart(ic),IpEnd(ic)
                        CMembership(ip) = ic    ! This has not been assigned correctly in the past, here is the first attempt
                        if ( MahalSqWeight(ic)*CKMahalSq(ip,ic) >= MahalSqWeight(icOther)*CKMahalSq(ip,icOther) ) then
                            CMembership(ip) = icOther
                            CKSize(ic) = CKSize(ic) - 1
                            CKSize(icOther) = CKSize(icOther) + 1
                            CKCenter(1:nd,ic) = CKCenter(1:nd,ic) - Point(1:nd,ip)
                            CKCenter(1:nd,icOther) = CKCenter(1:nd,icOther) + Point(1:nd,ip)
                            reclusteringNeeded = .true.
                        end if
                    end do

                    ic = 2
                    icOther = 1
                    do ip = IpStart(ic),IpEnd(ic)
                        CMembership(ip) = ic    ! This has not been assigned correctly in the past, here is the first attempt
                        if ( MahalSqWeight(ic)*CKMahalSq(ip,ic) > MahalSqWeight(icOther)*CKMahalSq(ip,icOther) ) then
                            CMembership(ip) = icOther
                            CKSize(ic) = CKSize(ic) - 1
                            CKSize(icOther) = CKSize(icOther) + 1
                            CKCenter(1:nd,ic) = CKCenter(1:nd,ic) - Point(1:nd,ip)
                            CKCenter(1:nd,icOther) = CKCenter(1:nd,icOther) + Point(1:nd,ip)
                            reclusteringNeeded = .true.
                        end if
                    end do

                    blockReclusteringNeeded: if (reclusteringNeeded) then  ! perform reassignment

                        blockIllegalClusterSize: if ( any(CKSize<ndPlusOne) ) then

                            minVolFailureCounter = minVolFailureCounter + 1
                            if (minVolFailureCounter>maxAllowedMinVolFailure) then
                                neOptimal = 1
                                CMembership = 1
                                CSize(1) = np
                                convergenceFailureCount = convergenceFailureCount + 1
                                !debugMinVolFailureCounter: block
                                !  write(*,*) 'minVolFailureCounter: ', minVolFailureCounter
                                !end block debugMinVolFailureCounter
                                return
                            end if

                            kmeansFailureCounter = 0
                            loopRunKmeans: do
                                call runKmeans(.true.,2,nd,np,Point,CKCenter,CKSize,CMembership,iteration,zcsIteration)
                                if ( any(CKSize<ndPlusOne) ) then
                                    kmeansFailureCounter = kmeansFailureCounter + 1
                                    if (kmeansFailureCounter>maxAllowedKmeansFailure) then
                                        neOptimal = 1
                                        CMembership = 1
                                        CSize(1) = np
                                        convergenceFailureCount = convergenceFailureCount + 1
                                        !debugKmeansFailureCounter: block
                                        !  write(*,*) 'kmeansFailureCounter: ', kmeansFailureCounter
                                        !end block debugKmeansFailureCounter
                                        return
                                    end if
                                    cycle loopRunKmeans
                                end if
                                exit loopRunKmeans
                            end do loopRunKmeans

                            self%Backup%LogVolNormed = huge(0._RK)
                            cycle loopRecursiveRunMinVolKernel

                        end if blockIllegalClusterSize

                        CKCenter(1:nd,1) = CKCenter(1:nd,1) / CKSize(1)  ! Now it is sum instead of mean
                        CKCenter(1:nd,2) = CKCenter(1:nd,2) / CKSize(2)  ! Now it is sum instead of mean

                        cycle loopRecursiveRunMinVolKernel

                    end if blockReclusteringNeeded

                    CKCenter(1:nd,1) = CKCenter(1:nd,1) / CKSize(1)  ! Now it is sum instead of mean
                    CKCenter(1:nd,2) = CKCenter(1:nd,2) / CKSize(2)  ! Now it is sum instead of mean

                    exit loopRecursiveRunMinVolKernel

                end do loopRecursiveRunMinVolKernel

                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                ! check whether the identified bounding ellipsoids are terminal
                !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                furtherClusteringCheck: if (sum(CKBELogVolNormed)<CBELogVolNormed(1) .or. CBELogVolNormed(1)>tightness*parLogVolNormed) then

                    ! at least two clusters is better than one. now search for more child clusters
                    CKneMax(1) = max( CKSize(1)/(ndEffective+1) , 1 )
                    CKneMax(2) = max( CKSize(2)/(ndEffective+1) , 1 )

                    IpStart(1) = 1             ; IpEnd(1) = CKSize(1)
                    IpStart(2) = CKSize(1) + 1 ; IpEnd(2) = np

                    do ic = 1,2   ! Search for more grand-child clusters

                        if (ic==1) then
                            icStart = 1
                            icEnd = CKneMax(1)
                        else
                            icStart = CKneOptimal(1) + 1
                            icEnd = CKneOptimal(1) + CKneMax(2)
                        end if

                        CBELogVolNormed(icStart) = CKBELogVolNormed(ic)
                        if (CKneMax(ic)>1) then ! .and. CKBELogVolNormed(ic)>1.1_RK*CKParLogVolNormed(ic)) then
                            call runMinVolKernel( CKneMax(ic), nd, CKSize(ic)               &
                                                , stabilizationRequested                    &
                                                , maxAllowedLoopRecursion                   &
                                                , maxAllowedKmeansFailure                   &
                                                , maxAllowedMinVolFailure                   &
                                                , mahalSqWeightExponent                     &
                                                , tightness                                 &
                                                , inclusionFraction                         &
                                                , CKParLogVolNormed(ic)                        &
                                                , logNormFac                                &
                                                , Point(1:nd,IpStart(ic):IpEnd(ic))         &
                                                , CKneOptimal(ic)                           &
                                                , CCenter(1:nd,icStart:icEnd)               &
                                                , CSize(icStart:icEnd)                      &
                                                , CBEInvCovMat(1:nd,1:nd,icStart:icEnd)     &
                                                , CBECholDiagLower(1:nd,1:nd,icStart:icEnd)       &
                                                , CBECholDiag(1:nd,icStart:icEnd)           &
                                                , CBELogVolNormed(icStart:icEnd)               &
                                                , CMembership(IpStart(ic):IpEnd(ic))        &
                                                , PointIndex(IpStart(ic):IpEnd(ic))         &
                                                , convergenceFailureCount                 &
                                                )
                            if (ic==2) CMembership(IpStart(ic):IpEnd(ic)) = CMembership(IpStart(ic):IpEnd(ic)) + icStart - 1
                        else
                            CKneOptimal(ic)=1
                        end if

                        if ( CKneOptimal(ic)==1 .or. CKneMax(ic)==1 ) then ! There is only one cluster here
                            CKneOptimal(ic) = 1
                            CCenter(1:nd,icStart) = CKCenter(1:nd,ic)
                            CSize(icStart) = CKSize(ic)
                            CBEInvCovMat(1:nd,1:nd,icStart) = CKBEInvCovMat(1:nd,1:nd,ic) / ScaleFacSq(ic)
                            CBECholDiagLower(1:nd,1:nd,icStart) = CKBECholDiagLower(1:nd,1:nd,ic)
                            CBECholDiag(1:nd,icStart) = CKBECholDiag(1:nd,ic)
                            ! Now rescale Cholesky factor and its diagonals, such that they describe the bounding ellipsoid of the cluster
                            do j = 1,nd-1
                                CBECholDiag(j,icStart) = CBECholDiag(j,icStart) * ScaleFac(ic)
                                do i = j+1,nd
                                    CBECholDiagLower(i,j,icStart) = CBECholDiagLower(i,j,icStart) * ScaleFac(ic)
                                end do
                            end do
                            CBECholDiag(nd,icStart) = CBECholDiag(nd,icStart) * ScaleFac(ic)
                        end if

                    end do

                    neOptimal = sum(CKneOptimal)

                else furtherClusteringCheck  ! one cluster is better

                    neOptimal = 1
                    CMembership = 1
                    CSize(1) = np

                end if furtherClusteringCheck


            end do loopPartition

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            !%%%% end loopPartition
            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        else blockPartition ! There is only one cluster

            neOptimal = 1

        end if blockPartition

        if (neOptimal==1) then
            do j = 1,nd-1 ! Now rescale Cholesky factor and its diagonals, such that they describe the bounding ellipsoid of the cluster
                self%CholDiagLower(j,0,1) = self%CholDiagLower(j,0,1) * scaleFac
                do i = j+1,nd
                    self%CholDiagLower(i,j,1) = self%CholDiagLower(i,j,1) * scaleFac
                end do
            end do
            self%CholDiagLower(nd,0,1) = self%CholDiagLower(nd,0,1) * scaleFac
            self%InvCovMat(1:nd,1:nd,1) = self%InvCovMat(1:nd,1:nd,1) / scaleFac**2
        end if

    contains

        function setSingleTerminalEllipsoid()
            implicit none
            self%Size(1) = npCurrent
            self%Membership = self%neOptimal
            self%neOptimal = self%neOptimal + 1
            self%convergenceFailureCount = self%convergenceFailureCount + 1
        end function setSingleTerminalEllipsoid

    end subroutine getMinVolPartition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module MinVolPartition_mod ! LCOV_EXCL_LINE