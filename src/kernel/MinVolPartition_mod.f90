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

    !> The `MinVolPartition_type` class. Partitions an input array `Point(nd,np)` with `nd` attributes and `np` observations (points).
    type :: MinVolPartition_type
        !type(Ellipsoid_type)        :: Backup
        integer(IK)                 :: neOptimal                !< The predicted optimal number of clusters.
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
    end type MinVolPartition_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine runMinVolPartition   ( ncMax, nd, np & ! LCOV_EXCL_LINE
                                    , maxKvolumeLoopRecursion & ! LCOV_EXCL_LINE
                                    , maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                    , maxAllowedMinVolFailure & ! LCOV_EXCL_LINE
                                    , stabilizationRequested & ! LCOV_EXCL_LINE
                                    , mahalSqWeightExponent & ! LCOV_EXCL_LINE
                                    , tightness & ! LCOV_EXCL_LINE
                                    , inclusionFraction & ! LCOV_EXCL_LINE
                                    , parLogVol & ! LCOV_EXCL_LINE
                                    ! inout
                                    , Point & ! LCOV_EXCL_LINE
                                    ! output
                                    , ncOptimal & ! LCOV_EXCL_LINE
                                    , CCenter & ! LCOV_EXCL_LINE
                                    , CSize & ! LCOV_EXCL_LINE
                                    , CBEInvCovMat & ! LCOV_EXCL_LINE
                                    , CBECholCov & ! LCOV_EXCL_LINE
                                    , CBECholDiag & ! LCOV_EXCL_LINE
                                    , CBEVolNormed & ! LCOV_EXCL_LINE
                                    , CMembership & ! LCOV_EXCL_LINE
                                    , PointIndex & ! LCOV_EXCL_LINE
                                    , convergenceFailureCounter & ! LCOV_EXCL_LINE
                                    )
        use Statistics_mod, only: getSamCholFac, getMahalSq
        use Matrix_mod, only: getInvMatFromCholFac
        implicit none
        integer(IK) , intent(in)    :: ncMax,nd,np                 ! maximum number of clusters, number of dimensions, number of points
        integer(IK) , intent(in)    :: maxKvolumeLoopRecursion
        integer(IK) , intent(in)    :: maxAllowedKmeansFailure
        integer(IK) , intent(in)    :: maxAllowedMinVolFailure
        logical     , intent(in)    :: stabilizationRequested
        real(RK)    , intent(in)    :: mahalSqWeightExponent
        real(RK)    , intent(in)    :: tightness                    ! the factor deremining how large the bounding ellipsoid can be, compared to the true volume of the ellipsoid.
        real(RK)    , intent(in)    :: inclusionFraction            ! the fraction of non-cluster-member points that ARE inside the cluster, to be considered in the estimation of the true volume of the cluster.
        real(RK)    , intent(in)    :: parLogVol                    ! input log-estimate of the total volume of the points.
        real(RK)    , intent(inout) :: Point(nd,np)                 ! the input data (points) to be clustered. on output, it is ordered.
        real(RK)    , intent(out)   :: CCenter(nd,ncMax)            ! Centers of the clusters (nd,ncMax)
        real(RK)    , intent(out)   :: CBEInvCovMat(nd,nd,ncMax)    ! CBE: Cluster Bounding Ellipsoid - On output: Full symmetric Inverse Covariance Matrices of the Bounding Ellipsoids of all clusters.
        real(RK)    , intent(out)   :: CBECholCov(nd,nd,ncMax)      ! CBE: Cluster Bounding Ellipsoid - On output: The lower triangle contains the Cholesky factor of the covariance matrix of the input points, the upper triangle maybe the CovMat itself (not needed). On output contains in the lower triangle, the Cholesky factor of all clusters found.
        real(RK)    , intent(out)   :: CBECholDiag(nd,ncMax)        ! CBE: Cluster Bounding Ellipsoid - On output: The diagonal elements of the Cholesky factor of the covariance matrix of the input points, on output, contains the cholesky factor of all found clusters.
        real(RK)    , intent(out)   :: CBEVolNormed(ncMax)          ! CBE: Cluster Bounding Ellipsoid - On output: Each element of the array contains the volume of the corresponding Cluster Bounding Ellipsoid.
        integer     , intent(out)   :: ncOptimal                    ! this variable is saved throughout recursive iterations of subroutine
        integer     , intent(out)   :: CSize(ncMax)                 ! cluster sizes (nCMax)
        integer     , intent(out)   :: CMembership(np)              ! the CMembership id of each point (representing the cluster number to which the point belongs)
        integer     , intent(out)   :: PointIndex(np)               ! Vector containing the scrambled indices of Point elements. It represent the mapping from the input Point elements to output scrambled Point vector.
        integer     , intent(inout) :: convergenceFailureCounter    ! number of times the algorithm failed to converge.
        integer                     :: i,j
        real(RK)                    :: parVolNormed                 ! LogVolume of parent ellipsoid
        real(RK)                    :: logNormFac, scaleFac
        !ncMax = np/(nd+1)

        CCenter(1:nd,1) = 0._RK
        do j = 1,np
          CCenter(1:nd,1) = CCenter(1:nd,1) + Point(1:nd,j)
          PointIndex(j) = j
        end do
        CCenter(1:nd,1) = CCenter(1:nd,1) / real(np,kind=RK)

        call getSamCholFac( nd, np, CCenter(1:nd,1), Point, CBECholCov(1:nd,1:nd,1), CBECholDiag(1:nd,1) )   ! get Cholesky Factor of cluster samples
        if (CBECholDiag(1,1)<0._RK) then
            write(*,*) 'getCholeskyFactor() failed in getSamCholFac() in runMinVolPartition()'
            write(*,*) 'Max cluster number, dimension, total points:'
            write(*,*) ncMax, nd, np
            stop
        end if
        CBEInvCovMat(1:nd,1:nd,1) = getInvMatFromCholFac(nd,CBECholCov(1:nd,1:nd,1),CBECholDiag(1:nd,1))
        scaleFac = sqrt( maxval(getMahalSq(nd, np, CCenter(1:nd,1), CBEInvCovMat(1:nd,1:nd,1), Point)) )
        CBEVolNormed(1) = sum(log(CBECholDiag(1:nd,1))) + nd*log(scaleFac)   ! Get the LogVolume of the first cluster. This will be used to normalize the volumes of all child clusters.
        if (CBEVolNormed(1)>parLogVol) then
            logNormFac = CBEVolNormed(1)
            CBEVolNormed(1) = 1._RK
            parVolNormed = exp( parLogVol - logNormFac )
        else    ! enlarge the bounding ellipsoid
            !write(*,*) 'no no no no!'
            logNormFac = parLogVol
            !scaleFac = scaleFac * exp( 0.5_RK*(logNormFac-CBEVolNormed(1))/nd )
            scaleFac = scaleFac * exp( (logNormFac-CBEVolNormed(1))/nd )
            CBEVolNormed(1) = 1._RK
            parVolNormed = 1._RK
        end if

        if (ncMax > 1 .and. CBEVolNormed(1)>tightness*parVolNormed) then
            !convergenceFailureCounter = 0
            call runMinVolKernel( ncMax = ncMax & ! LCOV_EXCL_LINE
                                , nd = nd & ! LCOV_EXCL_LINE
                                , np = np & ! LCOV_EXCL_LINE
                                , stabilizationRequested = stabilizationRequested & ! LCOV_EXCL_LINE
                                , maxKvolumeLoopRecursion = maxKvolumeLoopRecursion & ! LCOV_EXCL_LINE
                                , maxAllowedKmeansFailure = maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                , maxAllowedMinVolFailure = maxAllowedMinVolFailure & ! LCOV_EXCL_LINE
                                , mahalSqWeightExponent = mahalSqWeightExponent & ! LCOV_EXCL_LINE
                                , tightness = tightness & ! LCOV_EXCL_LINE
                                , inclusionFraction = inclusionFraction & ! LCOV_EXCL_LINE
                                , parVolNormed = parVolNormed & ! LCOV_EXCL_LINE
                                , logNormFac = logNormFac & ! LCOV_EXCL_LINE
                                , Point = Point & ! LCOV_EXCL_LINE
                                , ncOptimal = ncOptimal & ! LCOV_EXCL_LINE
                                , CCenter = CCenter & ! LCOV_EXCL_LINE
                                , CSize = CSize & ! LCOV_EXCL_LINE
                                , CBEInvCovMat = CBEInvCovMat & ! LCOV_EXCL_LINE
                                , CBECholCov = CBECholCov & ! LCOV_EXCL_LINE
                                , CBECholDiag = CBECholDiag & ! LCOV_EXCL_LINE
                                , CBEVolNormed = CBEVolNormed & ! LCOV_EXCL_LINE
                                , CMembership = CMembership & ! LCOV_EXCL_LINE
                                , PointIndex = PointIndex & ! LCOV_EXCL_LINE
                                , convergenceFailureCounter = convergenceFailureCounter & ! LCOV_EXCL_LINE
                                )
        else ! There is only one cluster
            ncOptimal = 1
        end if

        if (ncOptimal==1) then
            do j = 1,nd-1 ! Now rescale Cholesky factor and its diagonals, such that they describe the bounding ellipsoid of the cluster
                CBECholDiag(j,1) = CBECholDiag(j,1) * scaleFac
                do i = j+1,nd
                    CBECholCov(i,j,1) = CBECholCov(i,j,1) * scaleFac
                end do
            end do
            CBECholDiag(nd,1) = CBECholDiag(nd,1) * scaleFac
            CBEInvCovMat(1:nd,1:nd,1) = CBEInvCovMat(1:nd,1:nd,1) / scaleFac**2
        end if

    end subroutine runMinVolPartition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This code assumes that ncMax > 1 on input.
    ! important input requirements:
    recursive subroutine runMinVolKernel( ncMax & ! LCOV_EXCL_LINE
                                        , nd & ! LCOV_EXCL_LINE
                                        , np & ! LCOV_EXCL_LINE
                                        , stabilizationRequested & ! LCOV_EXCL_LINE
                                        , maxKvolumeLoopRecursion & ! LCOV_EXCL_LINE
                                        , maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                        , maxAllowedMinVolFailure & ! LCOV_EXCL_LINE
                                        , mahalSqWeightExponent & ! LCOV_EXCL_LINE
                                        , tightness & ! LCOV_EXCL_LINE
                                        , inclusionFraction & ! LCOV_EXCL_LINE
                                        , parVolNormed & ! LCOV_EXCL_LINE
                                        , logNormFac & ! LCOV_EXCL_LINE
                                        , Point & ! LCOV_EXCL_LINE
                                        , ncOptimal & ! LCOV_EXCL_LINE
                                        , CCenter & ! LCOV_EXCL_LINE
                                        , CSize & ! LCOV_EXCL_LINE
                                        , CBEInvCovMat & ! LCOV_EXCL_LINE
                                        , CBECholCov & ! LCOV_EXCL_LINE
                                        , CBECholDiag & ! LCOV_EXCL_LINE
                                        , CBEVolNormed & ! LCOV_EXCL_LINE
                                        , CMembership & ! LCOV_EXCL_LINE
                                        , PointIndex & ! LCOV_EXCL_LINE
                                        , convergenceFailureCounter & ! LCOV_EXCL_LINE
                                        )

        use Statistics_mod, only: getSamCholFac, getMahalSq
        use Matrix_mod, only: getInvMatFromCholFac
        use Kmeans_mod, only: runKmeans
        use Err_mod, only:  Err_type
        implicit none
        integer , save          :: ncall = 0
        integer , intent(in)    :: ncMax                             ! maximum number of clusters determined from the nd and np
        integer , intent(in)    :: nd,np                             ! number of dimensions, number of points
        logical , intent(in)    :: stabilizationRequested
        integer , intent(in)    :: maxKvolumeLoopRecursion
        integer , intent(in)    :: maxAllowedKmeansFailure
        integer , intent(in)    :: maxAllowedMinVolFailure
        real(RK), intent(in)    :: mahalSqWeightExponent
        real(RK), intent(in)    :: tightness                         ! the factor deremining how large the bounding ellipsoid can be, compared to the true volume of the ellipsoid.
        real(RK), intent(in)    :: inclusionFraction                 ! the fraction of non-cluster-member points that ARE inside the cluster, to be considered in the estimation of the true volume of the cluster.
        real(RK), intent(in)    :: parVolNormed                      ! parent's sqrt(det(invCovMat)), where invCovMat is the inverse covariance matrix of the bounding ellipsoid of the input points.
        real(RK), intent(in)    :: logNormFac                        ! parent's sqrt(det(invCovMat)), where invCovMat is the inverse covariance matrix of the bounding ellipsoid of the input points.
        integer , intent(out)   :: ncOptimal                         ! this variable is saved throughout recursive iterations of subroutine
        integer , intent(out)   :: CSize(ncMax)                      ! cluster sizes (nCMax)
        integer , intent(inout) :: convergenceFailureCounter         ! returns the number of times the algorithm failed to converge.
        integer , intent(inout) :: CMembership(np)                   ! the CMembership id of each point (representing the cluster number to which the point belongs)
        integer , intent(inout) :: PointIndex(np)                    ! Vector containing the scrambled indices of Point elements. It represent the mapping from the input Point elements to output scrambled Point vector.
        real(RK), intent(inout) :: Point(nd,np)                      ! the input data (points) to be clustered.
        real(RK), intent(inout) :: CCenter(nd,ncMax)                 ! Centers of the clusters (nd,ncMax)
        real(RK), intent(inout) :: CBECholCov(nd,nd,ncMax)           ! CBE: Cluster Bounding Ellipsoid - On input: The lower triangle contains the Cholesky factor of the covariance matrix of the input points, the upper triangle maybe the CovMat itself (not needed). On output contains in the lower triangle, the Cholesky factor of all clusters found.
        real(RK), intent(inout) :: CBEInvCovMat(nd,nd,ncMax)         ! CBE: Cluster Bounding Ellipsoid - On input: The lower triangle contains the Cholesky factor of the covariance matrix of the input points, the upper triangle maybe the CovMat itself (not needed). On output contains in the lower triangle, the Cholesky factor of all clusters found.
        real(RK), intent(inout) :: CBECholDiag(nd,ncMax)             ! CBE: Cluster Bounding Ellipsoid - On input: The diagonal elements of the Cholesky factor of the covariance matrix of the input points, on output, contains the cholesky factor of all found clusters.
        real(RK), intent(inout) :: CBEVolNormed(ncMax)               ! CBE: Cluster Bounding Ellipsoid - On input: The volume of the Bounding Ellipsoid of the parent cluster. On output: The Bounding volume of all child clusters.
        real(RK)                :: CKCenter(nd,2)                    ! Cluster Centers for nc=2 Kmeans algorithm
        real(RK)                :: CKBECholCov(nd,nd,2)              ! Covariance Matrices of the bounding ellipsoid of the Kmeans' 2 clusters, reported in the upper triangle. The lower triangle contains the Cholesky factors
        real(RK)                :: CKBEInvCovMat(nd,nd,2)            ! Inverse of the Covariance Matrices of the bounding ellipsoid of the Kmeans' 2 clusters.
        real(RK)                :: CKBECholDiag(nd,2)                ! determinants of the inverse covariance matrices of the bounding ellipsoids of the Kmeans' 2 clusters.
        real(RK)                :: CKMahalSq(np,2)                   ! Vector of Maximum Mahalanobis distances for each cluster
        real(RK)                :: CKBEVolNormed(2)                  ! New Cluster Bounding Ellisoid Volumes
        real(RK)                :: CKParVolNormed(2)                 ! New Cluster Volume estimates
        real(RK)                :: ScaleFacSq(2),ScaleFac(2)         ! The factor by which determinant must be multiplied in order for the ellipsoid to become a bounding ellipsoid.
        real(RK)                :: CKPoint(nd,np)                    ! Cluster Points in 2-means partition result
        real(RK)                :: MahalSqWeight(2)                  ! MahalSqWeight definition of Lu et al 2007
        integer                 :: CKPointIndex(np)                  ! Vector containing the scrambled indices of Point elements. It represent the mapping from the input Point elements to output scrambled Point vector.
        integer                 :: CKSize(2)                         ! cluster sizes for nc=2 Kmeans algorithm
        integer                 :: NumPointInCK(2)                   ! cluster sizes for nc=2 Kmeans algorithm
        integer                 :: CKSizeDummy(2)                    ! cluster sizes for nc=2 Kmeans algorithm
        integer                 :: CKncMax(2)                        ! maximum allowed number of clusters in each of the two K-means clusters.
        integer                 :: CKncOptimal(2)                    ! optimal number of clusters in each of the child clusters
        integer                 :: CMemberCounter(2)                 ! cluster member counter
        integer                 :: zcsIteration                      ! Counts the number of times zero sized clusters are found, zcs stands for zero cluster size
        integer                 :: IpStart(2),IpEnd(2)
        integer                 :: i,j,ip,ic,icOther,iteration,ndEffective
        integer                 :: icStart,icEnd,recursionCounter,kmeansFailureCounter,minVolFailureCounter
        logical                 :: reclusteringNeeded
        integer                 :: minClusterSize

        ! stabilization variables:

        real(RK)                :: CKCenterOld(nd,2)            ! Cluster Centers for nc=2 Kmeans algorithm
        real(RK)                :: CKBECholCovOld(nd,nd,2)      ! Covariance Matrices of the bounding ellipsoid of the Kmeans' 2 clusters, reported in the upper triangle. The lower triangle contains the Cholesky factors
        real(RK)                :: CKBEInvCovMatOld(nd,nd,2)    ! Inverse of the Covariance Matrices of the bounding ellipsoid of the Kmeans' 2 clusters.
        real(RK)                :: CKBECholDiagOld(nd,2)        ! determinants of the inverse covariance matrices of the bounding ellipsoids of the Kmeans' 2 clusters.
        real(RK)                :: CKNormedPoint(nd,np,2)       ! Vector of Maximum Mahalanobis distances for each cluster
        real(RK)                :: CKMinDistanceSq(np)          ! Vector of Maximum Mahalanobis distances for each cluster
        real(RK)                :: CKMahalSqOld(np,2)           ! Vector of Maximum Mahalanobis distances for each cluster
        real(RK)                :: CKBEVolNormedOld(2)          ! New Cluster Bounding Ellisoid Volumes
        real(RK)                :: ScaleFacSqOld(2)             ! The factor by which determinant must be multiplied in order for the ellipsoid to become a bounding ellipsoid.
        real(RK)                :: ScaleFacOld(2)               ! The factor by which determinant must be multiplied in order for the ellipsoid to become a bounding ellipsoid.
        real(RK)                :: potential                    ! The factor by which determinant must be multiplied in order for the ellipsoid to become a bounding ellipsoid.
        type(Err_type)          :: Err

        minClusterSize = nd + 1
        ndEffective = np / ncMax + 1
        ncall = ncall + 1

        kmeansFailureCounter = 0
        loopRunKmeans: do
            call runKmeans  ( nd = nd & ! LCOV_EXCL_LINE
                            , np = np & ! LCOV_EXCL_LINE
                            , nc = 2_IK & ! LCOV_EXCL_LINE
                            , Point = Point & ! LCOV_EXCL_LINE
                            , Size = CKSize & ! LCOV_EXCL_LINE
                            , Center = CKCenter & ! LCOV_EXCL_LINE
                            , Membership = CMembership & ! LCOV_EXCL_LINE
                            , NormedPoint = CKNormedPoint & ! LCOV_EXCL_LINE
                            , MinDistanceSq = CKMinDistanceSq & ! LCOV_EXCL_LINE
                            , potential = potential & ! LCOV_EXCL_LINE
                            , niter = iteration & ! LCOV_EXCL_LINE
                            , nzsci = zcsIteration & ! LCOV_EXCL_LINE
                            , Err = Err & ! LCOV_EXCL_LINE
                            !, niterMax & ! LCOV_EXCL_LINE
                            !, nzsciMax & ! LCOV_EXCL_LINE
                            !, relTol & ! LCOV_EXCL_LINE
                            )
            if ( any(CKSize<minClusterSize) ) then
                kmeansFailureCounter = kmeansFailureCounter + 1
                if (kmeansFailureCounter>maxAllowedKmeansFailure) then
                    ncOptimal = 1
                    CMembership = 1
                    CSize(1) = np
                    convergenceFailureCounter = convergenceFailureCounter + 1
                    return
                end if
                cycle loopRunKmeans
            end if
            exit loopRunKmeans
        end do loopRunKmeans

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        CKBEVolNormedOld = huge(0._RK)
        recursionCounter = 0
        minVolFailureCounter = 0

        loopRecursiveRunMinVolKernel: do    ! find two approximately minimum-volume sub-cllusters recursively

            recursionCounter = recursionCounter + 1
            if(mod(recursionCounter,maxKvolumeLoopRecursion)==0) then
                write(*,*) 'recursionCounter: ', recursionCounter
                ncOptimal = 1
                CMembership = 1
                CSize(1) = np
                convergenceFailureCounter = convergenceFailureCounter + 1
                return
            end if

            CMemberCounter(1) = 0
            CMemberCounter(2) = CKSize(1)
            do ip = 1,np
                if (CMembership(ip)==1) then
                    CMemberCounter(1) = CMemberCounter(1) + 1
                    CKPoint(1:nd,CMemberCounter(1)) = Point(1:nd,ip)
                    CKPointIndex(CMemberCounter(1)) = PointIndex(ip)
                else
                    CMemberCounter(2) = CMemberCounter(2) + 1
                    CKPoint(1:nd,CMemberCounter(2)) = Point(1:nd,ip)
                    CKPointIndex(CMemberCounter(2)) = PointIndex(ip)
                end if
            end do
            Point = CKPoint           ! Point is now ordered
            PointIndex = CKPointIndex ! PointIndex is now ordered

            IpStart(1) = 1             ; IpEnd(1) = CKSize(1)
            IpStart(2) = CKSize(1) + 1 ; IpEnd(2) = np

            do ic = 1,2   ! Get scaled Mahalanobis distances

                icOther = 3 - ic

                ! get Cholesky Factor of cluster samples

                call getSamCholFac  ( nd &
                                    , CKSize(ic) &
                                    , CKCenter(1:nd,ic) &
                                    , Point(1:nd,IpStart(ic):IpEnd(ic)) &
                                    , CKBECholCov(1:nd,1:nd,ic) &
                                    , CKBECholDiag(1:nd,ic) &
                                    )
                if (CKBECholDiag(1,ic)<0._RK) then
                    write(*,*) "    FATAL: runMinVolKernel() @ getSamCholFac() @ getCholeskyFactor() failed."
                    write(*,*) "           Cluster id, dimension, size:"
                    write(*,*) ic, nd, CKSize(ic)
                    stop
                end if

                CKBEInvCovMat(1:nd,1:nd,ic) = getInvMatFromCholFac( nd, CKBECholCov(1:nd,1:nd,ic), CKBECholDiag(1:nd,ic) )
                CKMahalSq(1:np,ic) = getMahalSq( nd, np, CKCenter(1:nd,ic), CKBEInvCovMat(1:nd,1:nd,ic), Point )
                ScaleFacSq(ic) = maxval(CKMahalSq(IpStart(ic):IpEnd(ic),ic))
                ScaleFac(ic) = sqrt( ScaleFacSq(ic) )
                CKMahalSq(1:np,ic) = CKMahalSq(1:np,ic) / ScaleFacSq(ic)
                CKBEVolNormed(ic) = exp( nd*log(ScaleFac(ic)) + sum(log(CKBECholDiag(1:nd,ic))) - logNormFac ) ! Get the LogVolume of the cluster
                NumPointInCK(ic) = CKSize(ic) + inclusionFraction*count(CKMahalSq(IpStart(icOther):IpEnd(icOther),ic)<1._RK)
                CKParVolNormed(ic) = real(NumPointInCK(ic),kind=RK) * parVolNormed / real(np,kind=RK)

                blockStabilityCheck: if (stabilizationRequested) then

                    if ( CKBEVolNormed(ic) > CKBEVolNormedOld(ic) ) then    ! This is to ensure the convergence of the algorithm

                        CKCenter(1:nd,ic) = CKCenterOld(1:nd,ic)
                        CKBECholCov(1:nd,1:nd,ic) = CKBECholCovOld(1:nd,1:nd,ic)
                        CKBECholDiag(1:nd,ic) = CKBECholDiagOld(1:nd,ic)
                        CKBEInvCovMat(1:nd,1:nd,ic) = CKBEInvCovMatOld(1:nd,1:nd,ic)
                        ScaleFacSq(ic) = ScaleFacSqOld(ic)
                        ScaleFac(ic) = ScaleFacOld(ic)
                        CKBEVolNormed(ic) = CKBEVolNormedOld(ic)

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
                        CKBECholCovOld(1:nd,1:nd,ic) = CKBECholCov(1:nd,1:nd,ic)
                        CKBECholDiagOld(1:nd,ic) = CKBECholDiag(1:nd,ic)
                        CKBEInvCovMatOld(1:nd,1:nd,ic) = CKBEInvCovMat(1:nd,1:nd,ic)
                        CKMahalSqOld(1:np,ic) = CKMahalSq(1:np,ic)
                        ScaleFacSqOld(ic) = ScaleFacSq(ic)
                        ScaleFacOld(ic) = ScaleFac(ic)
                        CKBEVolNormedOld(ic) = CKBEVolNormed(ic)

                    end if

                end if blockStabilityCheck

                !if (CKBEVolNormed(ic)<CKParVolNormed(ic)) then  ! enlarge the bounding ellipsoid
                !  !write(*,*) 'enlargement happening...'
                !  scaleFacSq(ic) = scaleFacSq(ic) * (CKParVolNormed(ic)/CKBEVolNormed(ic))**(2._RK/nd)  ! ATTN: (1._RK/nd) this should likely be (2._RK/nd)
                !  scaleFac(ic) = sqrt(scaleFacSq(ic))
                !  CKBEVolNormed(ic) = CKParVolNormed(ic)
                !end if

#if defined DEBUG_ENABLED
                debugVolNorm: block
                    real(RK) :: dummy
                    dummy = real(CKSize(ic),kind=RK)*parVolNormed/real(np,kind=RK)
                    if ( CKBEVolNormed(ic) < dummy ) then  ! enlarge the bounding ellipsoid
                        !write(*,*) 'enlargement happening...'
                        scaleFacSq(ic) = scaleFacSq(ic) * (dummy/CKBEVolNormed(ic))**(2._RK/nd)  ! ATTN: (1._RK/nd) this should likely be (2._RK/nd)
                        scaleFac(ic) = sqrt(scaleFacSq(ic))
                        CKBEVolNormed(ic) = dummy
                    end if
                end block debugVolNorm
#endif
              !end if blockStabilityCheck

            end do

            !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            ! Reassign points to clusters if needed
            MahalSqWeight = ( CKBEVolNormed / CKParVolNormed )**mahalSqWeightExponent
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

                blockIllegalClusterSize: if ( any(CKSize<minClusterSize) ) then

                    minVolFailureCounter = minVolFailureCounter + 1
                    if (minVolFailureCounter>maxAllowedMinVolFailure) then
                        ncOptimal = 1
                        CMembership = 1
                        CSize(1) = np
                        convergenceFailureCounter = convergenceFailureCounter + 1
                        !debugMinVolFailureCounter: block
                        !  write(*,*) 'minVolFailureCounter: ', minVolFailureCounter
                        !end block debugMinVolFailureCounter
                        return
                    end if

                    kmeansFailureCounter = 0
                    loopRerunKmeans: do
                        call runKmeans  ( nd = nd & ! LCOV_EXCL_LINE
                                        , np = np & ! LCOV_EXCL_LINE
                                        , nc = 2_IK & ! LCOV_EXCL_LINE
                                        , Point = Point & ! LCOV_EXCL_LINE
                                        , Size = CKSize & ! LCOV_EXCL_LINE
                                        , Center = CKCenter & ! LCOV_EXCL_LINE
                                        , Membership = CMembership & ! LCOV_EXCL_LINE
                                        , NormedPoint = CKNormedPoint & ! LCOV_EXCL_LINE
                                        , MinDistanceSq = CKMinDistanceSq & ! LCOV_EXCL_LINE
                                        , potential = potential & ! LCOV_EXCL_LINE
                                        , niter = iteration & ! LCOV_EXCL_LINE
                                        , nzsci = zcsIteration & ! LCOV_EXCL_LINE
                                        , Err = Err & ! LCOV_EXCL_LINE
                                        !, niterMax & ! LCOV_EXCL_LINE
                                        !, nzsciMax & ! LCOV_EXCL_LINE
                                        !, relTol & ! LCOV_EXCL_LINE
                                        )
                        if ( any(CKSize<minClusterSize) ) then
                        kmeansFailureCounter = kmeansFailureCounter + 1
                        if (kmeansFailureCounter>maxAllowedKmeansFailure) then
                            ncOptimal = 1
                            CMembership = 1
                            CSize(1) = np
                            convergenceFailureCounter = convergenceFailureCounter + 1
                            !debugKmeansFailureCounter: block
                            !  write(*,*) 'kmeansFailureCounter: ', kmeansFailureCounter
                            !end block debugKmeansFailureCounter
                            return
                        end if
                        cycle loopRerunKmeans
                        end if
                        exit loopRerunKmeans
                    end do loopRerunKmeans

                    CKBEVolNormedOld = huge(0._RK)
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

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        furtherClusteringCheck: if (sum(CKBEVolNormed)<CBEVolNormed(1) .or. CBEVolNormed(1)>tightness*parVolNormed) then ! at least two clusters is better than one. now search for more child clusters

            CKncMax(1) = max( CKSize(1)/(ndEffective+1) , 1 )
            CKncMax(2) = max( CKSize(2)/(ndEffective+1) , 1 )

            IpStart(1) = 1             ; IpEnd(1) = CKSize(1)
            IpStart(2) = CKSize(1) + 1 ; IpEnd(2) = np

            loopSearchForChildren: do ic = 1,2   ! Search for more grand-child clusters

                if (ic==1) then
                    icStart = 1
                    icEnd = CKncMax(1)
                else
                    icStart = CKncOptimal(1) + 1
                    icEnd = CKncOptimal(1) + CKncMax(2)
                end if

                CBEVolNormed(icStart) = CKBEVolNormed(ic)
                if (CKncMax(ic)>1) then ! .and. CKBEVolNormed(ic)>1.1_RK*CKParVolNormed(ic)) then
                    call runMinVolKernel( ncMax = CKncMax(ic) & ! LCOV_EXCL_LINE
                                        , nd = nd & ! LCOV_EXCL_LINE
                                        , np = CKSize(ic) & ! LCOV_EXCL_LINE
                                        , stabilizationRequested = stabilizationRequested & ! LCOV_EXCL_LINE
                                        , maxKvolumeLoopRecursion = maxKvolumeLoopRecursion & ! LCOV_EXCL_LINE
                                        , maxAllowedKmeansFailure = maxAllowedKmeansFailure & ! LCOV_EXCL_LINE
                                        , maxAllowedMinVolFailure = maxAllowedMinVolFailure & ! LCOV_EXCL_LINE
                                        , mahalSqWeightExponent = mahalSqWeightExponent & ! LCOV_EXCL_LINE
                                        , tightness = tightness & ! LCOV_EXCL_LINE
                                        , inclusionFraction = inclusionFraction & ! LCOV_EXCL_LINE
                                        , parVolNormed = CKParVolNormed(ic) & ! LCOV_EXCL_LINE
                                        , logNormFac = logNormFac & ! LCOV_EXCL_LINE
                                        , Point = Point(1:nd,IpStart(ic):IpEnd(ic)) & ! LCOV_EXCL_LINE
                                        , ncOptimal = CKncOptimal(ic) & ! LCOV_EXCL_LINE
                                        , CCenter = CCenter(1:nd,icStart:icEnd) & ! LCOV_EXCL_LINE
                                        , CSize = CSize(icStart:icEnd) & ! LCOV_EXCL_LINE
                                        , CBEInvCovMat = CBEInvCovMat(1:nd,1:nd,icStart:icEnd) & ! LCOV_EXCL_LINE
                                        , CBECholCov = CBECholCov(1:nd,1:nd,icStart:icEnd) & ! LCOV_EXCL_LINE
                                        , CBECholDiag = CBECholDiag(1:nd,icStart:icEnd) & ! LCOV_EXCL_LINE
                                        , CBEVolNormed = CBEVolNormed(icStart:icEnd) & ! LCOV_EXCL_LINE
                                        , CMembership = CMembership(IpStart(ic):IpEnd(ic)) & ! LCOV_EXCL_LINE
                                        , PointIndex = PointIndex(IpStart(ic):IpEnd(ic)) & ! LCOV_EXCL_LINE
                                        , convergenceFailureCounter = convergenceFailureCounter & ! LCOV_EXCL_LINE
                                        )
                    if (ic==2) CMembership(IpStart(ic):IpEnd(ic)) = CMembership(IpStart(ic):IpEnd(ic)) + icStart - 1
                else
                    CKncOptimal(ic)=1
                end if

                if ( CKncOptimal(ic)==1 .or. CKncMax(ic)==1 ) then ! There is only one cluster here
                    CKncOptimal(ic) = 1
                    CCenter(1:nd,icStart) = CKCenter(1:nd,ic)
                    CSize(icStart) = CKSize(ic)
                    CBEInvCovMat(1:nd,1:nd,icStart) = CKBEInvCovMat(1:nd,1:nd,ic) / ScaleFacSq(ic)
                    CBECholCov(1:nd,1:nd,icStart) = CKBECholCov(1:nd,1:nd,ic)
                    CBECholDiag(1:nd,icStart) = CKBECholDiag(1:nd,ic)
                    do j = 1,nd-1 ! Now rescale Cholesky factor and its diagonals, such that they describe the bounding ellipsoid of the cluster
                        CBECholDiag(j,icStart) = CBECholDiag(j,icStart) * ScaleFac(ic)
                        do i = j+1,nd
                            CBECholCov(i,j,icStart) = CBECholCov(i,j,icStart) * ScaleFac(ic)
                        end do
                    end do
                    CBECholDiag(nd,icStart) = CBECholDiag(nd,icStart) * ScaleFac(ic)
                end if

            end do loopSearchForChildren

            ncOptimal = sum(CKncOptimal)

        else furtherClusteringCheck  ! one cluster is better

            ncOptimal = 1
            CMembership = 1
            CSize(1) = np

        end if furtherClusteringCheck

        !write(*,*) 'ncMax: ', ncMax
        !write(*,*) 'ncOptimal: ', ncOptimal
        !write(*,*) 'CSize: ', CSize(1:ncOptimal)
        !write(*,*) 'CCenter: ', CCenter(:,1:ncOptimal)
        !write(*,*) 'Membership: ', CMembership(:)
        !read(*,*)

    end subroutine runMinVolKernel

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  end module MinVolPartition_mod