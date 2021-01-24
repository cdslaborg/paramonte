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

!> \brief This module contains the classes and procedures for various statistical computations.
!> \author Amir Shahmoradi

module Statistics_mod

    use Constants_mod, only: RK, IK
    use Err_mod, only: Err_type

    implicit none

    !logical, save :: paradramPrintEnabled = .false.
    !logical, save :: paradisePrintEnabled = .false.

    character(len=*), parameter :: MODULE_NAME = "@Statistics_mod"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface getSNormCDF
        module procedure :: getSNormCDF_SPR, getSNormCDF_DPR
    end interface getSNormCDF

    interface getBetaCDF
        module procedure :: getBetaCDF_SPR, getBetaCDF_DPR
    end interface getBetaCDF

    interface getBetaContinuedFraction
        module procedure :: getBetaContinuedFraction_SPR, getBetaContinuedFraction_DPR
    end interface getBetaContinuedFraction

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface getMean
        module procedure :: getMean_2D
    end interface getMean

    interface flatten
        module procedure :: flatten_2D
    end interface flatten

    interface getNormData
        module procedure :: getNormData_2D, normalizeWeightedData_2d
    end interface getNormData

    interface getVariance
        module procedure :: getVariance_1D, getVariance_2D
    end interface getVariance

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface getLogProbLogn
        module procedure :: getLogProbLognSP, getLogProbLognMP
    end interface getLogProbLogn

    interface getLogProbLognorm
        module procedure :: getLogProbLognSP, getLogProbLognMP
    end interface getLogProbLognorm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface getLogProbNormSP
        module procedure :: getLogProbNormSP_RK, getLogProbNormSP_CK
    end interface getLogProbNormSP

    interface getLogProbNormMP
        module procedure :: getLogProbNormMP_RK, getLogProbNormMP_CK
    end interface getLogProbNormMP

    interface getLogProbMVNSP
        module procedure :: getLogProbMVNSP_RK, getLogProbMVNSP_CK
    end interface getLogProbMVNSP

    interface getLogProbMVNMP
        module procedure :: getLogProbMVNMP_RK, getLogProbMVNMP_CK
    end interface getLogProbMVNMP

    interface getLogProbNorm
        module procedure :: getLogProbNormSP_RK, getLogProbNormMP_RK
        module procedure :: getLogProbNormSP_CK, getLogProbNormMP_CK
    end interface getLogProbNorm

    interface getLogProbMVN
        module procedure :: getLogProbMVNSP_RK, getLogProbMVNMP_RK
        module procedure :: getLogProbMVNSP_CK, getLogProbMVNMP_CK
    end interface getLogProbMVN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface getLogProbMixNorm
        module procedure :: getLogProbMixNormSP_RK, getLogProbMixNormSP_CK
        module procedure :: getLogProbMixNormMP_RK, getLogProbMixNormMP_CK
    end interface getLogProbMixNorm

    interface getLogProbMixMVN
        module procedure :: getLogProbMixMVNSP_RK, getLogProbMixMVNSP_CK
        module procedure :: getLogProbMixMVNMP_RK, getLogProbMixMVNMP_CK
    end interface getLogProbMixMVN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface getMahalSqSP
        module procedure :: getMahalSqSP_RK, getMahalSqSP_CK
    end interface getMahalSqSP

    interface getMahalSqMP
        module procedure :: getMahalSqMP_RK, getMahalSqMP_CK
    end interface getMahalSqMP

    interface getMahalSq
        module procedure :: getMahalSqSP_RK, getMahalSqMP_RK
        module procedure :: getMahalSqSP_CK, getMahalSqMP_CK
    end interface getMahalSq

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type :: ClusteredPoint_type
        integer(IK)                 :: nd
        integer(IK)                 :: np
        integer(IK)                 :: nc
        integer(IK)                 :: ndmin
        integer(IK)                 :: ndmax
        integer(IK)                 :: ncmin
        integer(IK)                 :: ncmax
        integer(IK)                 :: sizeMin
        integer(IK)                 :: sizeMax
        real(RK)                    :: etaMin
        real(RK)                    :: etaMax
        real(RK)                    :: stdMin
        real(RK)                    :: stdMax
        real(RK)                    :: centerMin
        real(RK)                    :: centerMax
        real(RK)    , allocatable   :: Eta(:)
        real(RK)    , allocatable   :: Std(:,:)
        real(RK)    , allocatable   :: Point(:,:)
        real(RK)    , allocatable   :: Center(:,:)
        real(RK)    , allocatable   :: ChoDia(:,:)
        real(RK)    , allocatable   :: LogVolume(:)
        real(RK)    , allocatable   :: ChoLowCovUpp(:,:,:)
        integer(IK) , allocatable   :: Membership(:)
        integer(IK) , allocatable   :: Size(:)
        character(:), allocatable   :: dist
        type(Err_type)              :: Err
    contains
        procedure, pass :: get => getClusteredPoint
        procedure, pass :: write => writeClusteredPoint
    end type ClusteredPoint_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    !> \brief
    !> Return the square of Mahalanobis distance for a single point. The output is a scalar variable.
    !>
    !> \param[in]   nd          :   The number of dimensions of the input `Point`.
    !> \param[in]   MeanVec     :   The mean vector of the sample.
    !> \param[in]   InvCovMat   :   The inverse covariance matrix of the sample.
    !> \param[in]   Point       :   The `Point` whose distance from the sample is to computed.
    !>
    !> \return
    !> `mahalSq` : The Mahalanobis distance squared of the point from
    !> the sample characterized by the input `MeanVec` and `InvCovMat`.
    pure function getMahalSqSP_RK(nd,MeanVec,InvCovMat,Point) result(mahalSq)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMahalSqSP_RK
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: MeanVec(nd)
        real(RK)   , intent(in) :: InvCovMat(nd,nd)        ! Inverse of the covariance matrix
        real(RK)   , intent(in) :: Point(nd)               ! input data points
        real(RK)                :: mahalSq
        real(RK)                :: NormedPoint(nd)
        NormedPoint = Point - MeanVec
        mahalSq = dot_product( NormedPoint , matmul(InvCovMat,NormedPoint) )
    end function getMahalSqSP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the square of Mahalanobis distances for an row-wise array of points.
    !>
    !> \param[in]   nd          :   The number of dimensions of the input `Point` array.
    !> \param[in]   np          :   The number of points in the input input `Point` array.
    !> \param[in]   MeanVec     :   The mean vector of length `nd` of the sample.
    !> \param[in]   InvCovMat   :   The inverse covariance matrix `(nd,nd)` of the sample.
    !> \param[in]   Point       :   The `(nd,np)` array of points whose distances squared from the sample are to computed.
    !>
    !> \return
    !> `MahalSq` : A vector of length `np` containing the squares of the Mahalanobis distances
    !> of the input points from the sample characterized by the input `MeanVec` and `InvCovMat`.
    !>
    !> \warning
    !> For the sake of preserving the purity and computational efficiency of the function,
    !> if the computation fails at any stage, the first element of output will be returned negative.
    pure function getMahalSqMP_RK(nd,np,MeanVec,InvCovMat,Point) result(MahalSq)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMahalSqMP_RK
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK), intent(in) :: nd,np
        real(RK), intent(in)    :: MeanVec(nd)
        real(RK), intent(in)    :: InvCovMat(nd,nd) ! Inverse of the covariance matrix
        real(RK), intent(in)    :: Point(nd,np)     ! input data points
        real(RK)                :: MahalSq(np)      ! function return
        real(RK)                :: NormedPoint(nd)
        integer(IK)             :: ip
        do ip = 1,np
            NormedPoint = Point(1:nd,ip) - MeanVec
            MahalSq(ip) = dot_product( NormedPoint , matmul(InvCovMat,NormedPoint) )
            if (MahalSq(ip)<0._RK) then
            ! LCOV_EXCL_START
                MahalSq(1) = -1._RK
                return
            end if
            ! LCOV_EXCL_STOP
        end do
    end function getMahalSqMP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the square of Mahalanobis distance for a single complex point. The output is a scalar variable.
    !>
    !> \param[in]   nd          :   The number of dimensions of the input `Point`.
    !> \param[in]   MeanVec     :   The mean vector of the sample.
    !> \param[in]   InvCovMat   :   The inverse covariance matrix of the sample.
    !> \param[in]   Point       :   The `Point` whose distance from the sample is to computed.
    !>
    !> \return
    !> `mahalSq` : The Mahalanobis distance squared of the point from
    !> the sample characterized by the input `MeanVec` and `InvCovMat`.
    pure function getMahalSqSP_CK(nd,MeanVec,InvCovMat,Point) result(mahalSq)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMahalSqSP_CK
#endif
        use Constants_mod, only: IK, RK, CK
        implicit none
        integer(IK), intent(in)  :: nd
        complex(CK), intent(in)  :: MeanVec(nd)
        complex(CK), intent(in)  :: InvCovMat(nd,nd) ! Inverse of the covariance matrix
        complex(CK), intent(in)  :: Point(nd)        ! input data points
        complex(CK)              :: mahalSq             ! function return
        mahalSq = sum( (Point-MeanVec) * matmul(InvCovMat,Point-MeanVec) )
    end function getMahalSqSP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the square of Mahalanobis distances for an row-wise array of complex-valued points.
    !>
    !> \param[in]   nd          :   The number of dimensions of the input `Point` array.
    !> \param[in]   np          :   The number of points in the input input `Point` array.
    !> \param[in]   MeanVec     :   The mean vector of length `nd` of the sample.
    !> \param[in]   InvCovMat   :   The inverse covariance matrix `(nd,nd)` of the sample.
    !> \param[in]   Point       :   The `(nd,np)` array of points whose distances squared from the sample are to computed.
    !>
    !> \return
    !> `MahalSq` : A vector of length `np` containing the squares of the Mahalanobis distances
    !> of the input points from the sample characterized by the input `MeanVec` and `InvCovMat`.
    !>
    !> \warning
    !> For the sake of preserving the purity and computational efficiency of the function,
    !> if the computation fails at any stage, the first element of output will be returned negative.
    pure function getMahalSqMP_CK(nd,np,MeanVec,InvCovMat,Point) result(MahalSq)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMahalSqMP_CK
#endif
        use Constants_mod, only: IK, RK, CK
        implicit none
        integer(IK), intent(in)  :: nd,np
        complex(CK), intent(in)  :: MeanVec(nd)
        complex(CK), intent(in)  :: InvCovMat(nd,nd) ! Inverse of the covariance matrix
        complex(CK), intent(in)  :: Point(nd,np)     ! input data points
        complex(CK)              :: MahalSq(np)         ! function return
        integer(IK)              :: ip
        do ip = 1,np
            MahalSq(ip) = sum( (Point(1:nd,ip)-MeanVec) * &
            matmul(InvCovMat,Point(1:nd,ip)-MeanVec) )
            if (real(MahalSq(ip))<0._RK) then
            ! LCOV_EXCL_START
                MahalSq(1) = (-1._RK, -1._RK)
                return
            end if
            ! LCOV_EXCL_STOP
        end do
    end function getMahalSqMP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogProbNormSP_RK(mean,inverseVariance,logSqrtInverseVariance,point) result(logProbNorm)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbNormSP_RK
#endif
        use Constants_mod, only: RK, LOGINVSQRT2PI
        implicit none
        real(RK), intent(in) :: mean,inverseVariance,logSqrtInverseVariance,point
        real(RK)             :: logProbNorm
        logProbNorm = LOGINVSQRT2PI + logSqrtInverseVariance - 0.5_RK * inverseVariance * (point-mean)**2
    end function getLogProbNormSP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogProbNormMP_RK(np,mean,inverseVariance,logSqrtInverseVariance,Point) result(logProbNorm)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbNormMP_RK
#endif
        use Constants_mod, only: LOGINVSQRT2PI
        implicit none
        integer(IK), intent(in) :: np
        real(RK)   , intent(in) :: mean,inverseVariance,logSqrtInverseVariance,Point(np)
        real(RK)                :: logProbNorm(np)
        logProbNorm = LOGINVSQRT2PI + logSqrtInverseVariance - 0.5_RK * inverseVariance * (Point-mean)**2
    end function getLogProbNormMP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! NOTE: if MahalSq computation fails, output probability will be returned as NullVal%RK from module Constant_mod.
    pure function getLogProbMVNSP_RK(nd,MeanVec,InvCovMat,logSqrtDetInvCovMat,Point) result(logProbNorm)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMVNSP_RK
#endif
        use Constants_mod, only: LOGINVSQRT2PI, NullVal
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: MeanVec(nd)
        real(RK)   , intent(in) :: InvCovMat(nd,nd)
        real(RK)   , intent(in) :: logSqrtDetInvCovMat
        real(RK)   , intent(in) :: Point(nd)
        real(RK)                :: logProbNorm, dummy
        dummy = getMahalSqSP_RK(nd,MeanVec,InvCovMat,Point)
        if (dummy<0._RK) then
            logProbNorm = NullVal%RK
        else
            logProbNorm = nd*LOGINVSQRT2PI + logSqrtDetInvCovMat - 0.5_RK * dummy
        end if
    end function getLogProbMVNSP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! NOTE: if MahalSq computation fails, output probability will be returned as NullVal%RK from module Constant_mod.
    pure function getLogProbMVNMP_RK(nd,np,MeanVec,InvCovMat,logSqrtDetInvCovMat,Point) result(logProbNorm)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMVNMP_RK
#endif
        use Constants_mod, only: LOGINVSQRT2PI, NullVal
        implicit none
        integer(IK), intent(in) :: nd,np
        real(RK)   , intent(in) :: MeanVec(nd)
        real(RK)   , intent(in) :: InvCovMat(nd,nd)
        real(RK)   , intent(in) :: logSqrtDetInvCovMat
        real(RK)   , intent(in) :: Point(nd,np)
        real(RK)                :: logProbNorm(np), Dummy(np)
        Dummy = getMahalSqMP_RK(nd,np,MeanVec,InvCovMat,Point)
        if (Dummy(1)<0._RK) then
            logProbNorm = NullVal%RK
        else
            logProbNorm = nd*LOGINVSQRT2PI + logSqrtDetInvCovMat - 0.5_RK * Dummy
        end if
    end function getLogProbMVNMP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getLogProbNormSP_CK(mean,inverseVariance,logSqrtInverseVariance,point) result(logProbNorm)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbNormSP_CK
#endif
        use Constants_mod, only: RK, CK, LOGINVSQRT2PI
        implicit none
        complex(CK), intent(in) :: mean,inverseVariance,logSqrtInverseVariance,point
        complex(CK)             :: logProbNorm
        logProbNorm = LOGINVSQRT2PI + logSqrtInverseVariance - 0.5_RK * inverseVariance * (point-mean)**2
    end function getLogProbNormSP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getLogProbNormMP_CK(np,mean,inverseVariance,logSqrtInverseVariance,Point) result(logProbNorm)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbNormMP_CK
#endif
        use Constants_mod, only: IK, RK, CK, LOGINVSQRT2PI
        implicit none
        integer(IK), intent(in) :: np
        complex(CK)   , intent(in) :: mean,inverseVariance,logSqrtInverseVariance,Point(np)
        complex(CK)                :: logProbNorm(np)
        logProbNorm = LOGINVSQRT2PI + logSqrtInverseVariance - 0.5_RK * inverseVariance * (Point-mean)**2
    end function getLogProbNormMP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getLogProbMVNSP_CK(nd,MeanVec,InvCovMat,logSqrtDetInvCovMat,Point) result(logProbMVN)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMVNSP_CK
#endif
        use Constants_mod, only: IK, RK, CK, LOGINVSQRT2PI, NullVal
        implicit none
        integer(IK), intent(in) :: nd
        complex(CK), intent(in) :: MeanVec(nd)
        complex(CK), intent(in) :: InvCovMat(nd,nd)
        complex(CK), intent(in) :: logSqrtDetInvCovMat
        complex(CK), intent(in) :: Point(nd)
        complex(CK)             :: logProbMVN, dummy
        dummy = getMahalSqSP(nd,MeanVec,InvCovMat,Point)
        if (real(dummy)<0._RK) then
            logProbMVN = NullVal%RK
        else
            logProbMVN  = nd*LOGINVSQRT2PI + logSqrtDetInvCovMat - 0.5_RK * dummy
        end if
    end function getLogProbMVNSP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getLogProbMVNMP_CK(nd,np,MeanVec,InvCovMat,logSqrtDetInvCovMat,Point) result(logProbMVN)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMVNMP_CK
#endif
        use Constants_mod, only: IK, RK, CK, LOGINVSQRT2PI, NullVal
        implicit none
        integer(IK), intent(in) :: nd,np
        complex(CK), intent(in) :: MeanVec(nd)
        complex(CK), intent(in) :: InvCovMat(nd,nd)
        complex(CK), intent(in) :: logSqrtDetInvCovMat
        complex(CK), intent(in) :: Point(nd,np)
        complex(CK)             :: logProbMVN(np), Dummy(np)
        Dummy = getMahalSqMP(nd,np,MeanVec,InvCovMat,Point)
        if (real(Dummy(1))<0._RK) then
            logProbMVN = NullVal%RK
        else
            logProbMVN  = nd*LOGINVSQRT2PI + logSqrtDetInvCovMat - 0.5_RK * Dummy
        end if
    end function getLogProbMVNMP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! SDSP stands for Single-Dimensional Gaussian mixture with Single Point input
    ! For a proper probability normalization, the sum of the amplitudes must equal one.
    function getLogProbMixNormSP_RK(nmode,LogAmplitude,MeanVec,InvCovMat,LogSqrtDetInvCovMat,point) result(logProbMixNorm)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMixNormSP_RK
#endif
        use Constants_mod, only: IK, RK, LOGTINY_RK
        implicit none
        integer(IK), intent(in) :: nmode
        real(RK)   , intent(in) :: LogAmplitude(nmode),MeanVec(nmode)
        real(RK)   , intent(in) :: InvCovMat(nmode),LogSqrtDetInvCovMat(nmode)
        real(RK)   , intent(in) :: point
        real(RK)                :: logProbMixNorm
        real(RK)                :: normFac,LogProb(nmode)
        integer(IK)             :: imode
        do imode = 1, nmode
            LogProb(imode)  = LogAmplitude(imode) + getLogProbNormSP_RK(MeanVec(imode),InvCovMat(imode),LogSqrtDetInvCovMat(imode),point)
        end do
        normFac = maxval(LogProb)
        LogProb = LogProb - normFac
        where(LogProb<LOGTINY_RK)
            LogProb = 0._RK
        elsewhere
            LogProb = exp(LogProb)
        end where
        logProbMixNorm = normFac + log(sum(LogProb))
    end function getLogProbMixNormSP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function getLogProbMixNormMP_RK(nmode,np,LogAmplitude,MeanVec,InvCovMat,LogSqrtDetInvCovMat,Point) result(logProbMixNorm)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMixNormMP_RK
#endif
        use Constants_mod, only: IK, RK, LOGTINY_RK
        implicit none
        integer(IK), intent(in) :: nmode,np
        real(RK)   , intent(in) :: LogAmplitude(nmode),MeanVec(nmode)
        real(RK)   , intent(in) :: InvCovMat(nmode),LogSqrtDetInvCovMat(nmode)
        real(RK)   , intent(in) :: Point(np)
        real(RK)                :: logProbMixNorm(np)
        real(RK)                :: NormFac(np),LogProb(nmode,np)
        integer(IK)             :: imode, ip
        do imode = 1, nmode
            LogProb(imode,1:np) = LogAmplitude(imode) + getLogProbNormMP_RK(np,MeanVec(imode),InvCovMat(imode),LogSqrtDetInvCovMat(imode),Point)
        end do
        NormFac = maxval(LogProb,dim=1)
        do ip = 1,np
            LogProb(1:nmode,ip) = LogProb(1:nmode,ip) - NormFac(ip)
            do imode = 1,nmode
                if ( LogProb(imode,ip) < LOGTINY_RK ) then
                    LogProb(imode,ip) = 0._RK
                else
                    LogProb(imode,ip) = exp( LogProb(imode,ip) )
                end if
            end do
            logProbMixNorm(ip) = NormFac(ip) + log(sum(LogProb(1:nmode,ip)))
        end do
    end function getLogProbMixNormMP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getLogProbMixMVNSP_RK(nmode,nd,LogAmplitude,MeanVec,InvCovMat,LogSqrtDetInvCovMat,Point) result(logProbMixMVN)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMixMVNSP_RK
#endif
        use Constants_mod, only: IK, RK, LOGTINY_RK
        implicit none
        integer(IK), intent(in) :: nmode,nd
        real(RK)   , intent(in) :: LogAmplitude(nmode), MeanVec(nd,nmode)
        real(RK)   , intent(in) :: InvCovMat(nd,nd,nmode), LogSqrtDetInvCovMat(nmode)
        real(RK)   , intent(in) :: Point(nd)
        real(RK)                :: logProbMixMVN
        real(RK)                :: normFac,LogProb(nmode)
        integer(IK)             :: imode
        do imode = 1, nmode
            LogProb(imode) = LogAmplitude(imode) + getLogProbMVNSP_RK(nd,MeanVec(1:nd,imode),InvCovMat(1:nd,1:nd,imode),LogSqrtDetInvCovMat(imode),Point)
        end do
        normFac = maxval(LogProb)
        LogProb = LogProb - normFac
        where(LogProb<LOGTINY_RK)
            LogProb = 0._RK
        elsewhere
            LogProb = exp(LogProb)
        end where
        logProbMixMVN = normFac + log(sum(LogProb))
    end function getLogProbMixMVNSP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function getLogProbMixMVNMP_RK(nmode,nd,np,LogAmplitude,MeanVec,InvCovMat,LogSqrtDetInvCovMat,Point) result(logProbMixMVN)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMixMVNMP_RK
#endif
        use Constants_mod, only: IK, RK, LOGTINY_RK
        implicit none
        integer(IK), intent(in) :: nmode,nd,np
        real(RK)   , intent(in) :: LogAmplitude(nmode),MeanVec(nd,nmode)
        real(RK)   , intent(in) :: InvCovMat(nd,nd,nmode), LogSqrtDetInvCovMat(nmode)
        real(RK)   , intent(in) :: Point(nd,np)
        real(RK)                :: logProbMixMVN(np)
        real(RK)                :: NormFac(np),LogProb(nmode,np)
        integer(IK)             :: imode, ip
        do imode = 1, nmode
            LogProb(imode,1:np) = LogAmplitude(imode) + &
            getLogProbMVNMP_RK(nd,np,MeanVec(1:nd,imode),InvCovMat(1:nd,1:nd,imode),LogSqrtDetInvCovMat(imode),Point)
        end do
        NormFac = maxval(LogProb,dim=1)
        do ip = 1,np
            LogProb(1:nmode,ip) = LogProb(1:nmode,ip) - NormFac(ip)
            do imode = 1,nmode
                if ( LogProb(imode,ip)<LOGTINY_RK ) then
                    LogProb(imode,ip) = 0._RK
                else
                    LogProb(imode,ip) = exp( LogProb(imode,ip) )
                end if
            end do
        logProbMixMVN(ip) = NormFac(ip) + log(sum(LogProb(1:nmode,ip)))
        end do
    end function getLogProbMixMVNMP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! SDSP stands for 1-dimensional Gaussian mixture with scalar input point
    function getLogProbMixNormSP_CK(nmode,LogAmplitude,MeanVec,InvCovMat,LogSqrtDetInvCovMat,point) result(logProbMixNorm)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMixNormSP_CK
#endif
        use Constants_mod, only: IK, RK, CK, LOGTINY_RK
        implicit none
        integer(IK), intent(in) :: nmode
        complex(CK), intent(in) :: LogAmplitude(nmode),MeanVec(nmode)
        complex(CK), intent(in) :: InvCovMat(nmode),LogSqrtDetInvCovMat(nmode)
        complex(CK), intent(in) :: point
        complex(CK)             :: logProbMixNorm
        complex(CK)             :: normFac,LogProb(nmode)
        integer(IK)             :: imode
        do imode = 1, nmode
            LogProb(imode) = LogAmplitude(imode) + getLogProbNorm(MeanVec(imode),InvCovMat(imode),LogSqrtDetInvCovMat(imode),point)
        end do
        normFac = maxval(real(LogProb))
        LogProb = LogProb - normFac
        where(real(LogProb)<LOGTINY_RK)
            LogProb = 0._RK
        elsewhere
            LogProb = exp(LogProb)
        end where
        logProbMixNorm = normFac + log(sum(LogProb))
    end function getLogProbMixNormSP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getLogProbMixNormMP_CK(nmode,np,LogAmplitude,MeanVec,InvCovMat,LogSqrtDetInvCovMat,Point) result(logProbMixNorm)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMixNormMP_CK
#endif
        use Constants_mod, only: IK, RK, CK, LOGTINY_RK
        implicit none
        integer(IK), intent(in) :: nmode,np
        complex(CK), intent(in) :: LogAmplitude(nmode),MeanVec(nmode)
        complex(CK), intent(in) :: InvCovMat(nmode),LogSqrtDetInvCovMat(nmode)
        complex(CK), intent(in) :: Point(np)
        complex(CK)             :: logProbMixNorm(np)
        complex(CK)             :: normFac(np),LogProb(nmode,np)
        integer(IK)             :: imode, ip
        do imode = 1, nmode
            LogProb(imode,1:np) = LogAmplitude(imode) + getLogProbNorm(np,MeanVec(imode),InvCovMat(imode),LogSqrtDetInvCovMat(imode),Point)
        end do
        normFac = maxval(real(LogProb),dim=1)
        do ip = 1,np
            LogProb(1:nmode,ip) = LogProb(1:nmode,ip) - normFac(ip)
            do imode = 1,nmode
                if ( real(LogProb(imode,ip)) < LOGTINY_RK ) then
                    LogProb(imode,ip) = 0._RK
                else
                    LogProb(imode,ip) = exp( LogProb(imode,ip) )
                end if
            end do
            logProbMixNorm(ip) = normFac(ip) + log(sum(LogProb(1:nmode,ip)))
        end do
    end function getLogProbMixNormMP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function getLogProbMixMVNSP_CK(nmode,nd,LogAmplitude,MeanVec,InvCovMat,LogSqrtDetInvCovMat,Point) result(logProbMixMVN)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMixMVNSP_CK
#endif
        use Constants_mod, only: IK, RK, CK, LOGTINY_RK
        implicit none
        integer(IK), intent(in) :: nmode,nd
        complex(CK), intent(in) :: LogAmplitude(nmode), MeanVec(nd,nmode)
        complex(CK), intent(in) :: InvCovMat(nd,nd,nmode), LogSqrtDetInvCovMat(nmode)
        complex(CK), intent(in) :: Point(nd)
        complex(CK)             :: logProbMixMVN
        complex(CK)             :: normFac,LogProb(nmode)
        integer(IK)             :: imode
        do imode = 1, nmode
            LogProb(imode) = LogAmplitude(imode) + getLogProbMVN(nd,MeanVec(1:nd,imode),InvCovMat(1:nd,1:nd,imode),LogSqrtDetInvCovMat(imode),Point)
        end do
        normFac = maxval(real(LogProb))
        LogProb = LogProb - normFac
        where(real(LogProb)<LOGTINY_RK)
            LogProb = 0._RK
        elsewhere
            LogProb = exp(LogProb)
        end where
        logProbMixMVN = normFac + log(sum(LogProb))
    end function getLogProbMixMVNSP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getLogProbMixMVNMP_CK(nmode,nd,np,LogAmplitude,MeanVec,InvCovMat,LogSqrtDetInvCovMat,Point) result(logProbMixMVN)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMixMVNMP_CK
#endif
        use Constants_mod, only: IK, RK, CK, LOGTINY_RK
        implicit none
        integer(IK), intent(in) :: nmode,nd,np
        complex(CK), intent(in) :: LogAmplitude(nmode),MeanVec(nd,nmode)
        complex(CK), intent(in) :: InvCovMat(nd,nd,nmode), LogSqrtDetInvCovMat(nmode)
        complex(CK), intent(in) :: Point(nd,np)
        complex(CK)             :: logProbMixMVN(np)
        complex(CK)             :: normFac(np),LogProb(nmode,np)
        integer(IK)             :: imode, ip
        do imode = 1, nmode
            LogProb(imode,1:np) = LogAmplitude(imode) + getLogProbMVN(nd,np,MeanVec(1:nd,imode),InvCovMat(1:nd,1:nd,imode),LogSqrtDetInvCovMat(imode),Point)
        end do
        normFac = maxval(real(LogProb),dim=1)
        do ip = 1,np
            LogProb(1:nmode,ip) = LogProb(1:nmode,ip) - normFac(ip)
            do imode = 1,nmode
                if ( real(LogProb(imode,ip))<LOGTINY_RK ) then
                    LogProb(imode,ip) = 0._RK
                else
                    LogProb(imode,ip) = exp( LogProb(imode,ip) )
                end if
            end do
            logProbMixMVN(ip) = normFac(ip) + log(sum(LogProb(1:nmode,ip)))
        end do
    end function getLogProbMixMVNMP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !include "Statistics_mod@MahalSq_RK.inc.f90"
    !include "Statistics_mod@MahalSq_CK.inc.f90"
    !include "Statistics_mod@LogProbGaus_RK.inc.f90"
    !include "Statistics_mod@LogProbGaus_CK.inc.f90"
    !include "Statistics_mod@LogProbGausMix_RK.inc.f90"
    !include "Statistics_mod@LogProbGausMix_CK.inc.f90"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Flatten the input `Point` array such that each element of the output
    !> `FlattenedPoint` array has the same unity weight as elements in the array.
    !>
    !> \param[in]       nd      :   The number of dimensions of the input sample.
    !> \param[in]       np      :   The number of points in the sample.
    !> \param[in]       Point   :   The array of shape `(nd,np)` containing the sample.
    !> \param[in]       Weight  :   The vector of length `np` containing the weights of points in the sample.
    !>                              The values of elements of Weight are allowed to be negative, in which case,
    !>                              the corresponding elements will be excluded from the output `FlattenedPoint`.
    !>
    !> \warning
    !> Note the shape of the input argument `Point(nd,np)`.
    !>
    !> \return
    !> `FlattenedPoint` : The flattened array whose elements all have the same weight.
    pure function flatten_2D(nd,np,Point,Weight) result(FlattenedPoint)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: flatten_2D
#endif
        ! the implementation for one-dimension is very concise and nice: mean = sum(Weight*Point) / sum(Weight)
        implicit none
        integer(IK), intent(in) :: np,nd            ! np: number of observations, nd: number of variables for each observation
        real(RK)   , intent(in) :: Point(nd,np)     ! Point is the data matrix
        integer(IK), intent(in) :: Weight(np)       ! sample weight
        integer(IK)             :: ip, iweight, sumWeight, counter
        real(RK), allocatable   :: FlattenedPoint(:,:)
        sumWeight = 0_IK
        do ip = 1, np
            if (Weight(ip)>0_IK) sumWeight = sumWeight + Weight(ip)
        end do
        allocate(FlattenedPoint(nd,sumWeight))
        counter = 0_IK
        do ip = 1, np
            do iweight = 1, Weight(ip)
                counter = counter + 1_IK
                FlattenedPoint(1:nd,counter) = Point(1:nd,ip)
            end do
        end do
    end function flatten_2D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the mean of a sample of multidimensional points.
    !>
    !> \param[in]       nd      :   The number of dimensions of the input sample.
    !> \param[in]       np      :   The number of points in the sample.
    !> \param[in]       Point   :   The array of shape `(nd,np)` containing the sample.
    !> \param[in]       Weight  :   The vector of length `np` containing the weights of points in the sample (**optional**, default = vector of ones).
    !>
    !> \warning
    !> Note the shape of the input argument `Point(nd,np)`.
    !>
    !> \return
    !> `Mean` : The output mean vector of length `nd`.
    pure function getMean_2D(nd,np,Point,Weight) result(Mean)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMean_2D
#endif
        ! the implementation for one-dimension is very concise and nice: mean = sum(Weight*Point) / sum(Weight)
        implicit none
        integer(IK), intent(in)             :: np,nd            ! np: number of observations, nd: number of variables for each observation
        real(RK)   , intent(in)             :: Point(nd,np)     ! Point is the data matrix
        integer(IK), intent(in), optional   :: Weight(np)       ! sample weight
        real(RK)                            :: Mean(nd)         ! output mean vector
        integer(IK)                         :: ip, sumWeight
        Mean = 0._RK
        if (present(Weight)) then
            sumWeight = 0_IK
            do ip = 1,np
                sumWeight = sumWeight + Weight(ip)
                Mean = Mean + Weight(ip) * Point(1:nd,ip)
            end do
            Mean = Mean / sumWeight
        else
            do ip = 1,np
                Mean = Mean + Point(1:nd,ip)
            end do
            Mean = Mean / real(np,kind=RK)
        end if
    end function getMean_2D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the normalized data with respect to the input mean vector of length `nd`.
    !>
    !> \param[in]       nd      :   The number of dimensions of the input sample.
    !> \param[in]       np      :   The number of points in the sample.
    !> \param[in]       Mean    :   The mean vector of length `nd`.
    !> \param[in]       Point   :   The array of shape `(nd,np)` containing the sample.
    !>
    !> \return
    !> `NormData` : The output normalized points array of shape `(np,nd)`.
    !>
    !> \remark
    !> Note the difference in the shape of the input `Point` vs. the output `NormData`.
    pure function getNormData_2D(nd,np,Mean,Point) result(NormData)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormData_2D
#endif
        implicit none
        integer(IK), intent(in)  :: np,nd           ! np is the number of observations, nd is the number of parameters for each observation
        real(RK)   , intent(in)  :: Mean(nd)        ! Mean vector
        real(RK)   , intent(in)  :: Point(nd,np)    ! Point is the matrix of the data, CovMat contains the elements of the sample covariance matrix
        real(RK)                 :: NormData(np,nd)
        integer(IK)              :: i
        do i = 1,np
            NormData(i,1:nd) = Point(1:nd,i) - Mean
        end do
    end function getNormData_2D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the normalized 2D data of size `(nd,np)` with respect to the mean of the data along the second dimension of length `np`.
    !>
    !> \param[in]   nd          :   The length of the input `Data` matrix along the first dimension.
    !> \param[in]   np          :   The length of the input `Data` matrix along the second dimension.
    !> \param[in]   Data        :   The input data series data vector.
    !> \param[in]   Weight      :   The vector of weights of the input data points (**optional**, default = array of ones).
    !> \param[in]   tenabled    :   A logical value that, if `.true.` will cause the output `NormData` to have transposed shape
    !>                              of the input `Point(nd,np)` matrix, that is `(np,nd)` (**optional**, default = `.false.`).
    !>
    !> \return
    !> `NormData` : The integrated autocorrelation (IAC) via the BatchMeans method.
    !>
    !> \remark
    !> Note that np must be large enough to get a meaningful answer.
    pure function normalizeWeightedData_2D(nd, np, Data, Weight, tenabled) result(NormData)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: normalizeWeightedData_2D
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK) , intent(in)            :: nd, np
        real(RK)    , intent(in)            :: Data(nd,np)
        integer(IK) , intent(in), optional  :: Weight(np)
        logical     , intent(in), optional  :: tenabled
        real(RK)    , allocatable           :: NormData(:,:)
        logical                             :: tenabledDefault
        real(RK)                            :: Mean(nd)
        real(RK)                            :: sumWeight
        integer(IK)                         :: id, ip

        tenabledDefault = .false.
        if (present(tenabled)) tenabledDefault = tenabled

        Mean = 0._RK
        if (present(Weight)) then
            sumWeight = 0._RK
            do ip = 1, np
                sumWeight = sumWeight + Weight(ip)
                Mean = Mean + Weight(ip) * Data(1:nd,ip)
            end do
        else
            sumWeight = np
            do ip = 1, np
                Mean = Mean + Data(1:nd,ip)
            end do
        end if
        Mean = Mean / sumWeight

        if (tenabledDefault) then
            allocate(NormData(np,nd))
            do concurrent(id = 1:nd)
                NormData(1:np,id) = Data(id,1:np) - Mean(id)
            end do
        else
            allocate(NormData(nd,np))
            do concurrent(id = 1:nd)
                NormData(id,1:np) = Data(id,1:np) - Mean(id)
            end do
        end if

    end function normalizeWeightedData_2D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the variance of the input vector of values of length `np`.
    !>
    !> \param[in]       np          :   The number of points in the sample.
    !> \param[in]       mean        :   The mean of the vector.
    !> \param[in]       Point       :   The array of shape `np` containing the sample.
    !> \param[in]       Weight      :   The vector of weights of the points in the sample (**optional**).
    !> \param[in]       sumWeight   :   The sum of the vector of weights (**optional**, if `Weight` is also missing).
    !>
    !> \return
    !> `variance` : The variance of the input sample.
    !>
    !> \warning
    !> If `Weight` is present, then `sumWeight` must be also present.
    !> Why? if mean is already given, it implies that sumWeight is also computed a priori.
    !>
    !> \remark
    !> One also use the concise Fortran syntax to achieve the same goal as this function:
    !> ```
    !> mean = sum(Weight*Point) / sum(Weight)
    !> variance = sum( (Weight*(Point-mean))**2 ) / (sum(Weight)-1)
    !> ```
    !> But the above concise version will be slightly slower as it involves three loops instead of two.
    pure function getVariance_1D(np,mean,Point,Weight,sumWeight) result(variance)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVariance_1D
#endif
        implicit none
        integer(IK), intent(in)             :: np                       ! np is the number of observations (points) whose variance is to be computed
        real(RK)   , intent(in)             :: mean                     ! the mean value of the vector Point
        real(RK)   , intent(in)             :: Point(np)                ! Point is the vector of data
        integer(IK), intent(in), optional   :: Weight(np), sumWeight    ! sample weight
        real(RK)                            :: variance                 ! output mean vector
        integer(IK)                         :: ip
        variance = 0._RK
        if (present(Weight)) then
            do ip = 1,np
                variance = variance + Weight(ip) * ( Point(ip) - mean )**2
            end do
            variance = variance / real(sumWeight-1_IK,kind=RK)
        else
            do ip = 1,np
                variance = variance + ( Point(ip) - mean )**2
            end do
            variance = variance / real(np-1_IK,kind=RK)
        end if
    end function getVariance_1D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the vector of variance of a set of `np` points of `nd` dimensions.
    !>
    !> \param[in]       nd          :   The number of dimensions of the input sample.
    !> \param[in]       np          :   The number of points in the sample.
    !> \param[in]       Mean        :   The mean vector of the sample.
    !> \param[in]       Point       :   The array of shape `(nd,np)` containing the sample.
    !> \param[in]       Weight      :   The vector of weights of the points in the sample of shape `(nd,np)` (**optional**).
    !>
    !> \return
    !> `Variance` : The vector of length `nd` of the variances of the input sample.
    pure function getVariance_2D(nd,np,Mean,Point,Weight) result(Variance)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getVariance_2D
#endif
        ! returns the variance of each row in Point weighted by the corresponding Weight if provided
        ! pass the Mean vector by calling getMean() or getMean_2D()
        implicit none
        integer(IK), intent(in)             :: nd, np           ! np is the number of observations (points) whose variance is to be computed
        real(RK)   , intent(in)             :: Mean(nd)         ! the Mean value of the vector Point
        real(RK)   , intent(in)             :: Point(nd,np)     ! Point is the vector of data
        integer(IK), intent(in), optional   :: Weight(np)       ! sample weight
        real(RK)                            :: Variance(nd)     ! output Mean vector
        integer(IK)                         :: ip, sumWeight
        Variance = 0._RK
        if (present(Weight)) then
            sumWeight = 0_IK
            do ip = 1,np
                sumWeight = sumWeight + Weight(ip)
                Variance = Variance + Weight(ip) * ( Point(1:nd,ip) - Mean )**2
            end do
            Variance = Variance / real(sumWeight-1_IK,kind=RK)
        else
            do ip = 1,np
                Variance = Variance + ( Point(1:nd,ip) - Mean )**2
            end do
            Variance = Variance / real(np-1_IK,kind=RK)
        end if
    end function getVariance_2D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the lower triangle Cholesky Factor of the covariance matrix of a set of points in the lower part of `CholeskyLower`.
    !> The upper part of `CholeskyLower`, including the diagonal elements of it, will contain the covariance matrix of the sample.
    !> The output argument `CholeskyDiago`, contains the diagonal terms of Cholesky Factor.
    !>
    !> \param[in]       nd              :   The number of dimensions of the input sample.
    !> \param[in]       np              :   The number of points in the sample.
    !> \param[in]       Mean            :   The mean vector of the sample.
    !> \param[in]       Point           :   The array of shape `(nd,np)` containing the sample.
    !> \param[out]      CholeskyLower   :   The output matrix of shape `(nd,nd)` whose lower triangle contains elements of the Cholesky factor.
    !>                                      The upper triangle of the matrix contains the covariance matrix of the sample.
    !> \param[out]      CholeskyDiago   :   The diagonal elements of the Cholesky factor.
    !>
    !> \todo
    !> The efficiency of this code can be further improved.
    subroutine getSamCholFac(nd,np,Mean,Point,CholeskyLower,CholeskyDiago)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSamCholFac
#endif
        use Matrix_mod, only: getCholeskyFactor
        implicit none
        integer(IK), intent(in)  :: nd,np                  ! np is the number of observations, nd is the number of parameters for each observation
        real(RK)   , intent(in)  :: Mean(nd)               ! Mean vector
        real(RK)   , intent(in)  :: Point(nd,np)           ! Point is the matrix of the data, CovMat contains the elements of the sample covariance matrix
        real(RK)   , intent(out) :: CholeskyLower(nd,nd)   ! Lower Cholesky Factor of the covariance matrix
        real(RK)   , intent(out) :: CholeskyDiago(nd)      ! Diagonal elements of the Cholesky factorization
        real(RK)                 :: NormData(np,nd), npMinusOneInverse
        integer(IK)              :: i,j

        do i = 1,np
            NormData(i,1:nd) = Point(1:nd,i) - Mean
        end do

        ! Only upper half of CovMat is needed
        npMinusOneInverse = 1._RK / real(np-1,kind=RK)
        do j = 1,nd
            do i = 1,j
                ! Get the covariance matrix elements: only the upper half of CovMat is needed
                CholeskyLower(i,j) = dot_product( NormData(1:np,i) , NormData(1:np,j) ) * npMinusOneInverse
            end do
        end do

        ! @todo: The efficiency can be improved by merging it with the above loops
        call getCholeskyFactor(nd, CholeskyLower, CholeskyDiago)

  end subroutine getSamCholFac

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the lower triangle Cholesky Factor of the covariance matrix of a set of points in the lower part of `CholeskyLower`.
    !> The upper part of `CholeskyLower`, including the diagonal elements of it, will contain the covariance matrix of the sample.
    !> The output argument `CholeskyDiago`, contains the diagonal terms of Cholesky Factor.
    !> This routine delivers the same functionality of [getSamCholFac()](@ref getsamcholfac).
    !> However, it is considerably faster, by a few factors, for `nd >> 10`.
    !> Furthermore, it requires significantly less runtime memory, by about half.
    !> The exact amount of speedup depends heavily on the architecture and memory layout of the runtime system.
    !>
    !> \param[in]       nd              :   The number of dimensions of the input sample.
    !> \param[in]       np              :   The number of points in the sample.
    !> \param[in]       Mean            :   The mean vector of the sample.
    !> \param[in]       Point           :   The array of shape `(nd,np)` containing the sample.
    !> \param[out]      CholeskyLower   :   The output matrix of shape `(nd,nd)` whose lower triangle contains elements of the Cholesky factor.
    !>                                      The upper triangle of the matrix contains the covariance matrix of the sample.
    !> \param[out]      CholeskyDiago   :   The diagonal elements of the Cholesky factor.
    !>
    !> \todo
    !> The efficiency of this code can be further improved.
    subroutine getSamCholFacHighDim(nd,np,Mean,Point,CholeskyLower,CholeskyDiago)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSamCholFacHighDim
#endif
        use Matrix_mod, only: getCholeskyFactor
        implicit none
        integer(IK), intent(in)  :: nd,np                  ! np is the number of observations, nd is the number of parameters for each observation
        real(RK)   , intent(in)  :: Mean(nd)               ! Mean vector
        real(RK)   , intent(in)  :: Point(nd,np)           ! Point is the matrix of the data, CovMat contains the elements of the sample covariance matrix
        real(RK)   , intent(out) :: CholeskyLower(nd,nd)   ! Lower Cholesky Factor of the covariance matrix
        real(RK)   , intent(out) :: CholeskyDiago(nd)      ! Diagonal elements of the Cholesky factorization
        real(RK)                 :: npMinusOneInverse
        real(RK)                 :: NormData(nd)
        integer(IK)              :: i, j, ip

        npMinusOneInverse = 1._RK / real(np-1,kind=RK)

        ! Compute the covariance matrix upper.

        CholeskyLower = 0._RK
        do ip = 1, np
            do j = 1, nd
                NormData(j) = Point(j,ip) - Mean(j)
                CholeskyLower(1:j,j) = CholeskyLower(1:j,j) + NormData(1:j) * NormData(j)
            end do
        end do
        CholeskyLower = CholeskyLower * npMinusOneInverse

        ! Compute the Cholesky factor lower.

        call getCholeskyFactor(nd, CholeskyLower, CholeskyDiago)

  end subroutine getSamCholFacHighDim

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the sample mean, covariance matrix, the Mahalanobis distances squared of the points with respect to the sample,
    !> and the square root of the determinant of the inverse covariance matrix of the sample.
    !>
    !> \param[in]       nd                  :   The number of dimensions of the input sample.
    !> \param[in]       np                  :   The number of points in the sample.
    !> \param[in]       Point               :   The array of shape `(np,nd)` containing the sample.
    !> \param[out]      CovMat              :   The output matrix of shape `(nd,nd)` representing the covariance matrix of the input data.
    !> \param[out]      Mean                :   The output mean vector of the sample.
    !> \param[out]      MahalSq             :   The output Mahalanobis distances squared of the points with respect to the sample (**optional**).
    !> \param[out]      InvCovMat           :   The output inverse covariance matrix of the input data (**optional**).
    !> \param[out]      sqrtDetInvCovMat    :   The output square root of the determinant of the inverse covariance matrix of the sample (**optional**).
    !>
    !> \warning
    !> If `sqrtDetInvCovMat` is present, then `MahalSq` and `InvCovMat` must be also present.
    !>
    !> \remark
    !> Note the shape of the input `Point(np,nd)`.
    !>
    !> \remark
    !> See also, [getSamCovMeanTrans](@ref getsamcovmeantrans).
    !>
    !> \remark
    !> For more information, see Geisser & Cornfield (1963) "Posterior distributions for multivariate normal parameters".
    !> Also, Box and Tiao (1973), "Bayesian Inference in Statistical Analysis" Page 421.
    !>
    !> \author
    !> Amir Shahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    subroutine getSamCovMean(np,nd,Point,CovMat,Mean,MahalSq,InvCovMat,sqrtDetInvCovMat)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSamCovMean
#endif
        use Matrix_mod, only: getInvPosDefMatSqrtDet
        implicit none
        integer(IK), intent(in)             :: np,nd                  ! np is the number of observations, nd is the number of parameters for each observation
        real(RK)   , intent(in)             :: Point(np,nd)           ! Point is the matrix of the data, CovMat contains the elements of the sample covariance matrix
        real(RK)   , intent(out)            :: CovMat(nd,nd)          ! Covariance matrix of the input data
        real(RK)   , intent(out)            :: Mean(nd)               ! Mean vector
        real(RK)   , intent(out), optional  :: MahalSq(np)            ! Vector of Mahalanobis Distances Squared, with respect to the mean position of the sample
        real(RK)   , intent(out), optional  :: InvCovMat(nd,nd)       ! Inverse Covariance matrix of the input data
        real(RK)   , intent(out), optional  :: sqrtDetInvCovMat       ! sqrt determinant of the inverse covariance matrix
        real(RK)                            :: NormData(np,nd)
        real(RK)                            :: DummyVec(nd)
        integer(IK)                         :: i,j

        do j = 1,nd
            Mean(j) = sum(Point(1:np,j)) / real(np,kind=RK)
            NormData(1:np,j) = Point(1:np,j) - Mean(j)
        end do
        do i = 1,nd
            do j = 1,nd
                CovMat(i,j) = dot_product(NormData(1:np,i),NormData(1:np,j))/real(np-1,kind=RK)
            end do
        end do

        if (present(sqrtDetInvCovMat)) then   ! Calculate InCovMat, MahalSq, and sqrt Determinant of InCovMat
            do j = 1,nd
                do i = 1,j
                InvCovMat(i,j) = CovMat(i,j)    ! Only the upper half of CovMat is needed
                end do
            end do
            call getInvPosDefMatSqrtDet(nd,InvCovMat,sqrtDetInvCovMat)
            do i = 1,np
                do j = 1,nd
                DummyVec(j) = dot_product(InvCovMat(1:nd,j),NormData(i,1:nd))
                end do
                MahalSq(i) = dot_product(NormData(i,1:nd),DummyVec)
            end do
        end if

    end subroutine getSamCovMean

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the sample mean, covariance matrix, the Mahalanobis distances squared of the points with respect to the sample,
    !> and the square root of the determinant of the inverse covariance matrix of the sample.
    !>
    !> \param[in]       nd                  :   The number of dimensions of the input sample.
    !> \param[in]       np                  :   The number of points in the sample.
    !> \param[in]       Point               :   The array of shape `(nd,np)` containing the sample.
    !> \param[out]      CovMat              :   The output matrix of shape `(nd,nd)` representing the covariance matrix of the input data.
    !> \param[out]      Mean                :   The output mean vector of the sample.
    !> \param[out]      MahalSq             :   The output Mahalanobis distances squared of the points with respect to the sample (**optional**).
    !> \param[out]      InvCovMat           :   The output inverse covariance matrix of the input data (**optional**).
    !> \param[out]      sqrtDetInvCovMat    :   The output square root of the determinant of the inverse covariance matrix of the sample (**optional**).
    !>
    !> \warning
    !> If `sqrtDetInvCovMat` is present, then `MahalSq` and `InvCovMat` must be also present.
    !>
    !> \remark
    !> Note the shape of the input `Point(nd,np)`.
    !>
    !> \remark
    !> This subroutine has the same functionality as [getSamCovMean](@ref getsamcovmean), with the only difference that input data is transposed here on input.
    !> Based on the preliminary benchmarks with Intel 2017 ifort, `getSamCovMean()` is slightly faster than `getSamCovMeanTrans()`.
    !>
    !> \remark
    !> For more information, see Geisser & Cornfield (1963) "Posterior distributions for multivariate normal parameters".
    !> Also, Box and Tiao (1973), "Bayesian Inference in Statistical Analysis" Page 421.
    !>
    !> \author
    !> Amir Shahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    subroutine getSamCovMeanTrans(np,nd,Point,CovMat,Mean,MahalSq,InvCovMat,sqrtDetInvCovMat)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSamCovMeanTrans
#endif

        use Matrix_mod, only: getInvPosDefMatSqrtDet
        implicit none
        integer(IK), intent(in)            :: np,nd                  ! np is the number of observations, nd is the number of parameters for each observation
        real(RK)   , intent(in)            :: Point(nd,np)           ! Point is the matrix of the data, CovMat contains the elements of the sample covariance matrix
        real(RK)   , intent(out)           :: CovMat(nd,nd)          ! Covariance matrix of the input data
        real(RK)   , intent(out)           :: Mean(nd)               ! Mean vector
        real(RK)   , intent(out), optional :: MahalSq(np)            ! Vector of Mahalanobis Distances Squared, with respect to the mean position of the sample
        real(RK)   , intent(out), optional :: InvCovMat(nd,nd)       ! Inverse Covariance matrix of the input data
        real(RK)   , intent(out), optional :: sqrtDetInvCovMat       ! sqrt determinant of the inverse covariance matrix
        real(RK)   , dimension(nd)         :: DummyVec
        real(RK)   , dimension(nd,np)      :: NormData
        integer(IK)                        :: i,j

        Mean = 0._RK
        do i = 1,np
            do j = 1,nd
                Mean(j) = Mean(j) + Point(j,i)
            end do
        end do
        Mean = Mean / real(np,kind=RK)

        do i = 1,np
            NormData(1:nd,i) = Point(1:nd,i) - Mean
        end do

        do i = 1,nd
            do j = 1,nd
                CovMat(i,j) = dot_product(NormData(i,1:np),NormData(j,1:np)) / real(np-1,kind=RK)
            end do
        end do

        if (present(sqrtDetInvCovMat)) then   ! Calculate InCovMat, MahalSq, and sqrt Determinant of InCovMat
            do j = 1,nd
                do i = 1,j
                InvCovMat(i,j) = CovMat(i,j)    ! Only the upper half of CovMat is needed
                end do
            end do
            call getInvPosDefMatSqrtDet(nd,InvCovMat,sqrtDetInvCovMat)
            do i = 1,np
                do j = 1,nd
                DummyVec(j) = dot_product(InvCovMat(1:nd,j),NormData(1:nd,i))
                end do
                MahalSq(i) = dot_product(NormData(1:nd,i),DummyVec)
                !MahalSq = dot_product(NormData(1:nd,i),DummyVec)
                !if (maxMahal<MahalSq) maxMahal = MahalSq
            end do
            !maxMahal = maxval(MahalSq)
        end if

    end subroutine getSamCovMeanTrans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the sample mean and the upper triangle of the covariance matrix of the input sample.
    !>
    !> \param[in]       nd                  :   The number of dimensions of the input sample.
    !> \param[in]       np                  :   The number of points in the sample.
    !> \param[in]       Point               :   The array of shape `(nd,np)` containing the sample.
    !> \param[out]      CovMatUpper         :   The output matrix of shape `(nd,nd)` whose upper triangle represents the covariance matrix of the input data.
    !> \param[out]      Mean                :   The output mean vector of the sample.
    !>
    !> \remark
    !> Note the shape of the input `Point(nd,np)`.
    !>
    !> \remark
    !> This subroutine has the same functionality as [getSamCovMeanTrans](@ref getsamcovmeantrans), with the only difference
    !> only the upper triangle of the covariance matrix is returned. Also, optional arguments are not available.
    !> This subroutine is specifically optimized for use in the ParaMonte samplers.
    !>
    !> \remark
    !> For more information, see Geisser & Cornfield (1963) "Posterior distributions for multivariate normal parameters".
    !> Also, Box and Tiao (1973), "Bayesian Inference in Statistical Analysis" Page 421.
    !>
    !> \author
    !> Amir Shahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    subroutine getSamCovUpperMeanTrans(np,nd,Point,CovMatUpper,Mean)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSamCovUpperMeanTrans
#endif
        use Matrix_mod, only: getInvPosDefMatSqrtDet
        implicit none
        integer(IK), intent(in)            :: np,nd                  ! np is the number of observations, nd is the number of parameters for each observation
        real(RK)   , intent(in)            :: Point(nd,np)           ! Point is the matrix of the data, CovMatUpper contains the elements of the sample covariance matrix
        real(RK)   , intent(out)           :: CovMatUpper(nd,nd)     ! Covariance matrix of the input data
        real(RK)   , intent(out)           :: Mean(nd)               ! Mean vector
        real(RK)                           :: npMinusOneInvReal, NormData(nd,np)
        integer(IK)                        :: i,j

        Mean = 0._RK
        do i = 1,np
            do j = 1,nd
                Mean(j) = Mean(j) + Point(j,i)
            end do
        end do
        Mean = Mean / real(np,kind=RK)

        do i = 1,np
            NormData(1:nd,i) = Point(1:nd,i) - Mean
        end do

        npMinusOneInvReal = 1._RK / real(np-1,kind=RK)
        do j = 1,nd
            do i = 1,j
                CovMatUpper(i,j) = dot_product(NormData(i,1:np),NormData(j,1:np)) * npMinusOneInvReal
            end do
        end do

  end subroutine getSamCovUpperMeanTrans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the mean and the upper triangle of the covariance matrix of the input *weighted* sample.
    !>
    !> \param[in]       nd                  :   The number of dimensions of the input sample.
    !> \param[in]       sumWeight           :   The sum of all sample weights.
    !> \param[in]       np                  :   The number of points in the sample.
    !> \param[in]       Point               :   The array of shape `(nd,np)` containing the sample.
    !> \param[in]       Weight              :   The integer array of length `np`, representing the weights of individual points in the sample.
    !> \param[out]      CovMatUpper         :   The output matrix of shape `(nd,nd)` whose upper triangle represents the covariance matrix of the input data.
    !> \param[out]      Mean                :   The output mean vector of the sample.
    !>
    !> \warning
    !> Note the shape of the input argument `Point(nd,np)`.
    !>
    !> \remark
    !> This subroutine has the same functionality as [getSamCovUpperMeanTrans](@ref getsamcovuppermeantrans), with the only difference
    !> only the upper triangle of the covariance matrix is returned. Also, optional arguments are not available.
    !> This subroutine is specifically optimized for use in the ParaMonte samplers.
    !>
    !> \remark
    !> For more information, see Geisser & Cornfield (1963) "Posterior distributions for multivariate normal parameters".
    !> Also, Box and Tiao (1973), "Bayesian Inference in Statistical Analysis" Page 421.
    !>
    !> \author
    !> Amir Shahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    subroutine getWeiSamCovUppMeanTrans(np,sumWeight,nd,Point,Weight,CovMatUpper,Mean)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getWeiSamCovUppMeanTrans
#endif

        use Matrix_mod, only: getInvPosDefMatSqrtDet
        implicit none
        integer(IK), intent(in)             :: np,nd,sumWeight        ! np is the number of observations, nd is the number of parameters for each observation
        integer(IK), intent(in)             :: Weight(np)             ! weights of the points
        real(RK)   , intent(in)             :: Point(nd,np)           ! Point is the matrix of the data, CovMatUpper contains the elements of the sample covariance matrix
        real(RK)   , intent(out)            :: CovMatUpper(nd,nd)     ! Covariance matrix of the input data
        real(RK)   , intent(out)            :: Mean(nd)               ! Mean vector
        real(RK)                            :: sumWeightMinusOneInvReal
        real(RK)                            :: NormData(nd,np)
        integer(IK)                         :: i,j,ip

        Mean = 0._RK
        do i = 1,np
            do j = 1,nd
                Mean(j) = Mean(j) + Weight(i)*Point(j,i)
            end do
        end do
        Mean = Mean / real(sumWeight,kind=RK)

        do i = 1,np
            NormData(1:nd,i) = Point(1:nd,i) - Mean
        end do

        sumWeightMinusOneInvReal = 1._RK / real(sumWeight-1,kind=RK)
        do j = 1,nd
            do i = 1,j
                CovMatUpper(i,j) = 0
                do ip = 1,np
                    CovMatUpper(i,j) = CovMatUpper(i,j) + Weight(ip)*NormData(i,ip)*NormData(j,ip)
                end do
                CovMatUpper(i,j) = CovMatUpper(i,j) * sumWeightMinusOneInvReal
            end do
        end do

    end subroutine getWeiSamCovUppMeanTrans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Given two input sample means and covariance matrices, return the combination of them as a single mean and covariance matrix.
    !>
    !> \param[in]       nd          :   The number of dimensions of the input sample.
    !> \param[in]       npA         :   The number of points in sample `A`.
    !> \param[in]       MeanVecA    :   The mean vector of sample `A`.
    !> \param[in]       CovMatA     :   The covariance matrix of sample `A`.
    !> \param[in]       npB         :   The number of points in sample `B`.
    !> \param[in]       MeanVecB    :   The mean vector of sample `B`.
    !> \param[in]       CovMatB     :   The covariance matrix of sample `B`.
    !> \param[out]      MeanVecAB   :   The output mean vector of the combined sample.
    !> \param[out]      CovMatAB    :   The output covariance matrix of the combined sample.
    !>
    !> \author
    !> Amir Shahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    !>
    !> \remark
    !> An exact implementation of this algorithm which needs only the upper triangles of the input matrices and
    !> yields only the upper triangle of the covariance matrix is given in [mergeMeanCovUpper](@ref mergemeancovupper).
    !> The alternative implementation is much more efficient, by a factor of 6-7 with all compiler optimization flags on.
    !>
    ! This subroutine uses a recursion equation similar to http://stats.stackexchange.com/questions/97642/how-to-combine-sample-means-and-sample-variances
    ! See also: O’Neill, (2014), "Some Useful Moment Results in Sampling Problems".
    ! See also: https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Covariance
    ! See also: https://stats.stackexchange.com/questions/43159/how-to-calculate-pooled-variance-of-two-groups-given-known-group-variances-mean
    subroutine mergeMeanCov(nd,npA,MeanVecA,CovMatA,npB,MeanVecB,CovMatB,MeanVecAB,CovMatAB)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: mergeMeanCov
#endif
        implicit none
        integer(IK), intent(in)       :: nd
        integer(IK), intent(in)       :: npA,npB
        real(RK)   , intent(in)       :: MeanVecA(nd),CovMatA(nd,nd)
        real(RK)   , intent(in)       :: MeanVecB(nd),CovMatB(nd,nd)
        real(RK)   , intent(out)      :: MeanVecAB(nd),CovMatAB(nd,nd)
        real(RK)   , dimension(nd,1)  :: MeanMatA,MeanMatB,MeanMat
        real(RK)                      :: DistanceSq(nd,nd)
        real(RK)                      :: npABinverse
        integer(IK)                   :: npAB

        npAB = npA + npB
        npABinverse = 1._RK / real(npAB, kind=RK)
        MeanMatA(1:nd,1) = MeanVecA
        MeanMatB(1:nd,1) = MeanVecB

        ! Compute the new Mean

        MeanVecAB = ( npA * MeanVecA + npB * MeanVecB ) * npABinverse
        MeanMat(1:nd,1) = MeanVecAB

        ! Compute the new Covariance matrix

        DistanceSq = matmul( (MeanMatA-MeanMatB), transpose((MeanMatA-MeanMatB)) ) * npA * npB * npABinverse
        CovMatAB = ( (npA-1) * CovMatA + (npB-1) * CovMatB + DistanceSq ) / real(npAB-1, kind=RK)

    end subroutine mergeMeanCov

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    !> \brief
!    !> Given two input sample means and covariance matrices, return the combination of them as a single mean and covariance matrix.
!    !>
!    !> \param[in]       nd              :   The number of dimensions of the input sample.
!    !> \param[in]       npA             :   The number of points in sample `A`.
!    !> \param[in]       MeanVecA        :   The mean vector of sample `A`.
!    !> \param[in]       CovMatUpperA    :   The covariance matrix of sample `A`.
!    !> \param[in]       npB             :   The number of points in sample `B`.
!    !> \param[in]       MeanVecB        :   The mean vector of sample `B`.
!    !> \param[in]       CovMatUpperB    :   The covariance matrix of sample `B`.
!    !> \param[out]      MeanVecAB       :   The output mean vector of the combined sample.
!    !> \param[out]      CovMatUpperAB   :   The output covariance matrix of the combined sample.
!    !>
!    !> \todo
!    !> The efficiency of this algorithm might still be improved by converting the upper triangle covariance matrix to a packed vector.
!    !>
!    !> \remark
!    !> This subroutine is the same as [mergeMeanCov](@ref mergeMeanCov), with the **important difference** that only the
!    !> upper triangles and diagonals of the input covariance matrices need to be given by the user: `CovMatUpperA`, `CovMatUpperB`
!    !> This alternative implementation is 6-7 times faster, with all compiler optimization flags on.
!    !>
!    !> \author
!    !> Amir Shahmoradi, Nov 24, 2020, 4:19 AM, Dallas, TX
!    subroutine mergeMeanCovUpperSlow(nd,npA,MeanVecA,CovMatUpperA,npB,MeanVecB,CovMatUpperB,MeanVecAB,CovMatUpperAB)
!#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
!        !DEC$ ATTRIBUTES DLLEXPORT :: mergeMeanCovUpperSlow
!#endif
!        implicit none
!        integer(IK), intent(in)             :: nd
!        integer(IK), intent(in)             :: npA,npB
!        real(RK)   , intent(in)             :: MeanVecA(nd),CovMatUpperA(nd,nd)
!        real(RK)   , intent(in)             :: MeanVecB(nd),CovMatUpperB(nd,nd)
!        real(RK)   , intent(out)            :: MeanVecAB(nd),CovMatUpperAB(nd,nd)
!        real(RK)                            :: npABinverse, npAnpB2npAB
!        real(RK)                            :: npA2npAB, npB2npAB
!        real(RK)                            :: MeanVecDiffAB(nd)
!        integer(IK)                         :: npAB, i, j
!
!        npAB = npA + npB
!        npABinverse = 1._RK / real(npAB, kind=RK)
!        npAnpB2npAB = npA * npB * npABinverse
!        npA2npAB = npA * npABinverse
!        npB2npAB = npB * npABinverse
!
!        ! Compute the new Mean and Covariance matrix
!
!        do j = 1, nd
!            MeanVecDiffAB(j) = MeanVecA(j) - MeanVecB(j)
!            !MeanVecAB(j) = ( npA * MeanVecA(j) + npB * MeanVecB(j) ) * npABinverse
!            MeanVecAB(j) = npA2npAB * MeanVecA(j) + npB2npAB * MeanVecB(j)
!            do i = 1, j
!                CovMatUpperAB(i,j) = ( (npA-1) * CovMatUpperA(i,j) + (npB-1) * CovMatUpperB(i,j) + MeanVecDiffAB(i) * MeanVecDiffAB(j) * npAnpB2npAB ) / real(npAB-1, kind=RK)
!            end do
!        end do
!
!
!    end subroutine mergeMeanCovUpperSlow

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Given two input sample means and covariance matrices, return the combination of them as a single mean and covariance matrix.
    !>
    !> \param[in]       nd              :   The number of dimensions of the input sample.
    !> \param[in]       npA             :   The number of points in sample `A`.
    !> \param[in]       MeanVecA        :   The mean vector of sample `A`.
    !> \param[in]       CovMatUpperA    :   The covariance matrix of sample `A`.
    !> \param[in]       npB             :   The number of points in sample `B`.
    !> \param[in]       MeanVecB        :   The mean vector of sample `B`.
    !> \param[in]       CovMatUpperB    :   The covariance matrix of sample `B`.
    !> \param[out]      MeanVecAB       :   The output mean vector of the combined sample.
    !> \param[out]      CovMatUpperAB   :   The output covariance matrix of the combined sample.
    !>
    !> \remark
    !> This subroutine is the same as [mergeMeanCov](@ref mergemeancov), with the **important difference** that only the
    !> upper triangles and diagonals of the input covariance matrices need to be given by the user: `CovMatUpperA`, `CovMatUpperB`
    !> This alternative implementation is 6-7 times faster, with all compiler optimization flags on.
    !> In addition, all computational coefficients are predefined in this implementation,
    !> resulting in an extra 10%-15% efficiency gain.
    !>
    !> \todo
    !> The efficiency of this algorithm might still be improved by converting the upper triangle covariance matrix to a packed vector.
    !>
    !> \author
    !> Amir Shahmoradi, Nov 24, 2020, 4:19 AM, Dallas, TX
    subroutine mergeMeanCovUpper(nd,npA,MeanVecA,CovMatUpperA,npB,MeanVecB,CovMatUpperB,MeanVecAB,CovMatUpperAB)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: mergeMeanCovUpper
#endif
        implicit none
        integer(IK), intent(in)             :: nd
        integer(IK), intent(in)             :: npA,npB
        real(RK)   , intent(in)             :: MeanVecA(nd),CovMatUpperA(nd,nd)
        real(RK)   , intent(in)             :: MeanVecB(nd),CovMatUpperB(nd,nd)
        real(RK)   , intent(out)            :: MeanVecAB(nd),CovMatUpperAB(nd,nd)
        real(RK)                            :: MeanVecDiffAB(nd)
        real(RK)                            :: npAnpB2npAB2npABMinusOne
        real(RK)                            :: npAMinusOne2npABMinusOne
        real(RK)                            :: npBMinusOne2npABMinusOne
        real(RK)                            :: npABMinusOneInverse
        real(RK)                            :: npABinverse
        real(RK)                            :: npA2npAB
        real(RK)                            :: npB2npAB
        integer(IK)                         :: npAB, i, j

        npAB = npA + npB
        npABinverse = 1._RK / real(npAB, kind=RK)
        npABMinusOneInverse = 1._RK / real(npAB-1, kind=RK)
        npAnpB2npAB2npABMinusOne = npA * npB * npABinverse * npABMinusOneInverse
        npA2npAB = npA * npABinverse
        npB2npAB = npB * npABinverse
        npAMinusOne2npABMinusOne = (npA - 1_IK) * npABMinusOneInverse
        npBMinusOne2npABMinusOne = (npB - 1_IK) * npABMinusOneInverse

        ! Compute the new Mean and Covariance matrix

        do j = 1, nd
            MeanVecDiffAB(j) = MeanVecA(j) - MeanVecB(j)
            MeanVecAB(j) = npA2npAB * MeanVecA(j) + npB2npAB * MeanVecB(j)
            do i = 1, j
                CovMatUpperAB(i,j)  = npAMinusOne2npABMinusOne * CovMatUpperA(i,j) &
                                    + npBMinusOne2npABMinusOne * CovMatUpperB(i,j) &
                                    + npAnpB2npAB2npABMinusOne * MeanVecDiffAB(i) * MeanVecDiffAB(j)
            end do
        end do


    end subroutine mergeMeanCovUpper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Given two input sample means and covariance matrices, return the combination of them as a single mean and covariance matrix.
    !>
    !> \param[in]       nd              :   The number of dimensions of the input sample.
    !> \param[in]       npA             :   The number of points in sample `A`.
    !> \param[in]       MeanVecA        :   The mean vector of sample `A`.
    !> \param[in]       CovMatUpperA    :   The covariance matrix of sample `A`.
    !> \param[in]       npB             :   The number of points in sample `B`.
    !> \param[inout]    MeanVecB        :   The mean vector of sample `B`.
    !> \param[inout]    CovMatUpperB    :   The covariance matrix of sample `B`.
    !>
    !> \remark
    !> This subroutine is the same as [mergeMeanCovUpper](@ref mergemeancovupper), with the **important difference** that
    !> the resulting output mean and covariance matrices are written to the input arguments `MeanVecB`, `CovMatUpperB`.
    !> This alternative implementation results in another extra 15%-20% efficiency gain. This result is based on
    !> the benchmarks with Intel Fortran compiler 19.4 with all compiler optimization flags on.
    !>
    !> \todo
    !> The efficiency of this algorithm might still be improved by converting the upper triangle covariance matrix to a packed vector.
    !>
    !> \author
    !> Amir Shahmoradi, Nov 25, 2020, 1:00 AM, Dallas, TX
    subroutine mergeMeanCovUpperDense(nd,npA,MeanVecA,CovMatUpperA,npB,MeanVecB,CovMatUpperB)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: mergeMeanCovUpperDense
#endif
        implicit none
        integer(IK), intent(in)             :: nd
        integer(IK), intent(in)             :: npA,npB
        real(RK)   , intent(in)             :: MeanVecA(nd),CovMatUpperA(nd,nd)
        real(RK)   , intent(inout)          :: MeanVecB(nd),CovMatUpperB(nd,nd)
        real(RK)                            :: MeanVecDiffAB(nd)
        real(RK)                            :: npAnpB2npAB2npABMinusOne
        real(RK)                            :: npAMinusOne2npABMinusOne
        real(RK)                            :: npBMinusOne2npABMinusOne
        real(RK)                            :: npABMinusOneInverse
        real(RK)                            :: npABinverse
        real(RK)                            :: npA2npAB
        real(RK)                            :: npB2npAB
        integer(IK)                         :: npAB, i, j

        npAB = npA + npB
        npABinverse = 1._RK / real(npAB, kind=RK)
        npABMinusOneInverse = 1._RK / real(npAB-1, kind=RK)
        npAnpB2npAB2npABMinusOne = npA * npB * npABinverse * npABMinusOneInverse
        npA2npAB = npA * npABinverse
        npB2npAB = npB * npABinverse
        npAMinusOne2npABMinusOne = (npA - 1_IK) * npABMinusOneInverse
        npBMinusOne2npABMinusOne = (npB - 1_IK) * npABMinusOneInverse

        ! Compute the new Mean and Covariance matrix

        do j = 1, nd
            MeanVecDiffAB(j) = MeanVecA(j) - MeanVecB(j)
            MeanVecB(j) = npA2npAB * MeanVecA(j) + npB2npAB * MeanVecB(j)
            do i = 1, j
                CovMatUpperB(i,j)   = npAMinusOne2npABMinusOne * CovMatUpperA(i,j) &
                                    + npBMinusOne2npABMinusOne * CovMatUpperB(i,j) &
                                    + npAnpB2npAB2npABMinusOne * MeanVecDiffAB(i) * MeanVecDiffAB(j)
            end do
        end do

    end subroutine mergeMeanCovUpperDense

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return a random standard Gaussian deviate with zero mean and unit variance.
    !>
    !> \return
    !> `randGaus` : The random standard Gaussian deviate with zero mean and unit variance.
    !>
    !> \remark
    !> See also, Numerical Recipes in Fortran, by Press et al. (1990)
    !>
    !> \author
    !> Amir Shahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    function getRandGaus() result(randGaus)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandGaus
#endif

        implicit none
        integer(IK), save :: iset=0
        real(RK)   , save :: gset
        real(RK)          :: fac,rsq,vec(2)
        real(RK)          :: randGaus

        if (iset == 0_IK) then
            do
                !block
                !integer :: i, n
                !real(RK) :: unifrnd(30)
                !!integer, dimension(:), allocatable :: seed
                !if (paradramPrintEnabled .or. paradisePrintEnabled) then
                !    !do i = 1, 22
                !    call random_number(unifrnd)
                !    write(*,"(*(g0,:,'"//new_line("a")//"'))") unifrnd
                !    !end do
                !    !call random_seed(size = n); allocate(seed(n))
                !    !call random_seed(get = seed)
                !    !write(*,"(*(g0,:,' '))") seed
                !    !write(*,"(*(g0,:,' '))") StateOld
                !    !write(*,"(*(g0,:,' '))") StateNew
                !    !write(*,"(*(g0,:,' '))") CholeskyLower
                !    !write(*,"(*(g0,:,' '))") domainCheckCounter
                !    paradisePrintEnabled = .false.
                !    paradramPrintEnabled = .false.
                !end if
                !end block
                call random_number(vec)
                vec = 2._RK*vec - 1._RK
                rsq = vec(1)**2 + vec(2)**2
                if ( rsq > 0._RK .and. rsq < 1._RK ) exit
            end do
            fac = sqrt(-2._RK*log(rsq)/rsq)
            gset = vec(1)*fac
            randGaus = vec(2)*fac
            iset = 1_IK
        else
            randGaus = gset
            iset = 0_IK
        endif

    end function getRandGaus

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return a random Gaussian deviate with the given mean and standard deviation.
    !>
    !> \param[in]       mean    :   The mean of the Gaussian distribution.
    !> \param[in]       std     :   The standard deviation of the Gaussian distribution. It must be a positive real number.
    !>
    !> \return
    !> `randNorm` : A normally distributed deviate with the given mean and standard deviation.
    !>
    !> \author
    !> Amir Shahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    function getRandNorm(mean,std) result(randNorm)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandNorm
#endif
        implicit none
        real(RK), intent(in)    :: mean, std
        real(RK)                :: randNorm
        randNorm = mean + std * getRandGaus()
    end function getRandNorm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return a log-normally distributed deviate with the given mean and standard deviation.
    !>
    !> \param[in]       mean    :   The mean of the Lognormal distribution.
    !> \param[in]       std     :   The standard deviation of the Lognormal distribution. It must be a positive real number.
    !>
    !> \return
    !> `randLogn` : A Lognormally distributed deviate with the given mean and standard deviation.
    !>
    !> \author
    !> Amir Shahmoradi, Oct 16, 2009, 11:14 AM, Michigan
    function getRandLogn(mean,std) result(randLogn)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandLogn
#endif
        implicit none
        real(RK), intent(in)    :: mean, std
        real(RK)                :: randLogn
        randLogn = exp( mean + std*getRandGaus() )
    end function getRandLogn

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This subroutine is legacy and slow. use getRandMVN() in this same module.
    ! Given the mean vector MeanVec and the covariance matrix CovMat, this subroutine generates a random vector x (of length nd>=2)
    ! from an nd-dimensional multivariate normal distribution.
    ! First a vector of nd standard normal random deviates is generated,
    ! then this vector is multiplied by the Cholesky decomposition of the covariance matrix,
    ! then the vector MeanVec is added to the resulting vector, and is stored in the output vector x.
    ! ATTENTION: Only the upper half of the covariance matrix (plus the diagonal terms) need to be given in the input.
    ! in the ouput, the upper half and diagonal part will still be the covariance matrix, while the lower half will be
    ! the Cholesky decomposition elements (excluding its diagonal terms that are provided only in the vector Diagonal).
    ! USES choldc.f90, getRandGaus.f90
    !> Amir Shahmoradi, March 22, 2012, 2:21 PM, IFS, UTEXAS
    subroutine getMVNDev(nd,MeanVec,CovMat,X)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMVNDev
#endif

        use iso_fortran_env, only: output_unit
        use Matrix_mod, only: getCholeskyFactor

        implicit none
        integer(IK), intent(in)  :: nd
        real(RK)   , intent(in)  :: MeanVec(nd), CovMat(nd,nd)
        real(RK)   , intent(out) :: X(nd)
        real(RK)                 :: CholeskyLower(nd,nd), Diagonal(nd), DummyVec(nd)
        integer(IK)              :: i

        CholeskyLower = CovMat
        call getCholeskyFactor(nd,CholeskyLower,Diagonal)
        if (Diagonal(1)<0._RK) then
        ! LCOV_EXCL_START
            write(output_unit,"(A)") "getCholeskyFactor() failed in getMVNDev()"
            error stop
        end if
        ! LCOV_EXCL_STOP
        do i=1,nd
            DummyVec(i) = getRandGaus()
            x(i) = DummyVec(i) * Diagonal(i)
        end do
        do i=2,nd
            x(i) = x(i) + dot_product(CholeskyLower(i,1:i-1),DummyVec(1:i-1))
        end do
        x = x + MeanVec

  end subroutine getMVNDev

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This subroutine is legacy and slow. use getRandMVU() in this same module.
    ! Given the mean vector MeanVec and the covariance matrix CovMat, this subroutine generates a random vector X (of length nd>=2)
    ! from an nd-dimensional multivariate ellipsoidal uniform distribution, such that getMVUDev() is randomly distributed inside the nd-dimensional ellipsoid.
    ! ATTENTION: Only the upper half of the covariance matrix (plus the diagonal terms) need to be given in the input.
    ! in the ouput, the upper half and diagonal part will still be the covariance matrix, while the lower half will be
    ! the Cholesky decomposition elements (excluding its diagonal terms that are provided only in the vector Diagonal).
    ! USES getCholeskyFactor.f90, getRandGaus.f90
    !> Amir Shahmoradi, April 25, 2016, 2:21 PM, IFS, UTEXAS
    subroutine getMVUDev(nd,MeanVec,CovMat,X)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMVUDev
#endif

        use Matrix_mod, only: getCholeskyFactor

        implicit none
        integer(IK), intent(in)  :: nd
        real(RK)   , intent(in)  :: MeanVec(nd)
        real(RK)   , intent(in)  :: CovMat(nd,nd)
        real(RK)   , intent(out) :: X(nd)
        real(RK)                 :: Diagonal(nd), DummyVec(nd), CholeskyLower(nd,nd), dummy
        integer(IK)              :: i

        CholeskyLower = CovMat
        call getCholeskyFactor(nd,CholeskyLower,Diagonal)
        if (Diagonal(1)<0._RK) then
        ! LCOV_EXCL_START
            error stop
            !call abortProgram( output_unit , 1 , 1 , 'Statitistics@getMVUDev()@getCholeskyFactor() failed.' )
        end if
        ! LCOV_EXCL_STOP
        do i=1,nd
            DummyVec(i) = getRandGaus()
        end do
        call random_number(dummy)
        dummy = (dummy**(1._RK/real(nd,kind=RK)))/norm2(DummyVec)  ! Now DummyVec is a uniformly drawn point from inside of nd-D sphere.
        DummyVec = dummy*DummyVec

        ! Now transform this point to a point inside the ellipsoid.
        do i=1,nd
            X(i) = DummyVec(i)*Diagonal(i)
        end do

        do i=2,nd
            X(i) = X(i) + dot_product(CholeskyLower(i,1:i-1),DummyVec(1:i-1))
        end do

        X = X + MeanVec

    end subroutine getMVUDev

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return a MultiVariate Normal (MVN) random vector with the given mean and
    !> covariance matrix represented by the the input Cholesky factorization.
    !>
    !> \param[in]   nd              :   The number of dimensions of the MVN distribution.
    !> \param[in]   MeanVec         :   The mean vector of the MVN distribution.
    !> \param[in]   CholeskyLower   :   The Cholesky lower triangle of the covariance matrix of the MVN distribution.
    !> \param[in]   Diagonal        :   The Diagonal elements of the Cholesky lower triangle of the covariance matrix of the MVN distribution.
    !>
    !> \return
    !> `RandMVN` : The randomly generated MVN vector.
    !>
    !> \author
    !> Amir Shahmoradi, April 23, 2017, 12:36 AM, ICES, UTEXAS
    function getRandMVN(nd,MeanVec,CholeskyLower,Diagonal) result(RandMVN)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandMVN
#endif
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: MeanVec(nd)
        real(RK)   , intent(in) :: CholeskyLower(nd,nd), Diagonal(nd)   ! Cholesky lower triangle and its diagonal terms, calculated from the input CovMat.
        real(RK)                :: RandMVN(nd), dummy
        integer(IK)             :: j,i
        RandMVN = MeanVec
        do j = 1,nd
            dummy = getRandGaus()
            RandMVN(j) = RandMVN(j) + Diagonal(j) * dummy
            do i = j+1,nd
                RandMVN(i) = RandMVN(i) + CholeskyLower(i,j) * dummy
            end do
        end do
    end function getRandMVN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    ! Given an input Mean vector of length nd, Covariance Matrix of dimension (nd,nd), as well as a vector of integers representing
!    ! the variables (corresponding to CovMat columns) that are given
!    ! This subroutine gives out a conditional Multivariate Normal Random deviate.
!    ! random p-tivariate normal deviate, given that the first pg variables x1 are given (i.e. fixed).
!    ! For a review of Multivariate Normal distribution: Applied Multivariate Statistical Analysis, Johnson, Wichern, 1998, 4th ed.
!    !> Amir Shahmoradi, Oct 20, 2009, 9:12 PM, MTU
!    function getCondRandMVN(nd,MeanVec,CovMat,nIndIndx,IndIndx) result(CondRandMVN)
!        use Matrix_mod, only: getRegresCoef
!        implicit none
!        integer(IK), intent(in) :: nd, nIndIndx, IndIndx(nIndIndx)
!        real(RK)   , intent(in) :: MeanVec(nd), CovMat(nd,nd)
!        real(RK)                :: CondRandMVN(nd),
!        integer(IK)             :: j, i
!    end function getCondRandMVN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Given the Cholesky Lower triangle and diagonals of a given covariance matrix, this function return one point uniformly
    !> randomly drawn from inside of an `nd`-ellipsoid, whose `nd` elements `RandMVU(i), i=1:nd` are guaranteed
    !> to be in the range:
    !> ```
    !> MeanVec(i) - sqrt(CovMat(i,i)) < RandMVU(i) < MeanVec(i) + sqrt(CovMat(i,i))
    !> ```
    !>
    !> \param[in]   nd              :   The number of dimensions of the MVU distribution.
    !> \param[in]   MeanVec         :   The mean vector of the MVU distribution.
    !> \param[in]   CholeskyLower   :   The Cholesky lower triangle of the covariance matrix of the MVU distribution.
    !> \param[in]   Diagonal        :   The Diagonal elements of the Cholesky lower triangle of the covariance matrix of the MVU distribution.
    !>
    !> \return
    !> `RandMVU` : The randomly generated MVU vector.
    !>
    !> \author
    !> Amir Shahmoradi, April 23, 2017, 1:36 AM, ICES, UTEXAS
    function getRandMVU(nd,MeanVec,CholeskyLower,Diagonal) result(RandMVU)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandMVU
#endif
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: MeanVec(nd)
        real(RK)   , intent(in) :: CholeskyLower(nd,nd) ! Cholesky lower triangle, calculated from the input CovMat.
        real(RK)   , intent(in) :: Diagonal(nd)         ! Cholesky diagonal terms, calculated from the input CovMat.
        real(RK)                :: RandMVU(nd), dummy, DummyVec(nd), sumSqDummyVec
        integer(IK)             :: i,j
        sumSqDummyVec = 0._RK
        do j=1,nd
            DummyVec(j) = getRandGaus()
            sumSqDummyVec = sumSqDummyVec + DummyVec(j)**2
        end do
        call random_number(dummy)
        dummy = dummy**(1._RK/nd) / sqrt(sumSqDummyVec)
        DummyVec = DummyVec * dummy  ! a uniform random point from inside of nd-sphere
        RandMVU = MeanVec
        do j = 1,nd
            RandMVU(j) = RandMVU(j) + Diagonal(j) * DummyVec(j)
            do i = j+1,nd
                RandMVU(i) = RandMVU(i) + CholeskyLower(i,j) * DummyVec(j)
            end do
        end do
    end function getRandMVU

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return `.true.` if the input `NormedPoint` (normalized with respect to the center of the target ellipsoid) is
    !> within or on the boundary of the ellipsoid whose boundary is described by the representative matrix
    !> \f$ \Sigma \f$ (`RepMat`), such that,
    !> \f{equation}{
    !> X^T ~ \Sigma^{-1} ~ X = 1 ~,
    !> \f}
    !> for all \f$X\f$ on the boundary.
    !>
    !> \param[in]   nd          :   The number of dimensions of the ellipsoid (the size of `NormedPoint`).
    !> \param[in]   NormedPoint :   The input point, normalized with respect to the center of the ellipsoid,
    !>                                  whose location with respect to the ellipsoid boundary is to be determined.
    !> \param[in]   InvRepMat   :   The inverse of the representative matrix of the target ellipsoid.
    !>
    !> \return
    !> `isInsideEllipsoid` : The logical value indicating whether the input point is inside or on the boundary of the target ellipsoid.
    !>
    !> \remark
    !> Note that the input matrix is the inverse of `RepMat`: `InvRepMat`.
    !>
    !> \author
    !> Amir Shahmoradi, April 23, 2017, 1:36 AM, ICES, UTEXAS
    pure function isInsideEllipsoid(nd,NormedPoint,InvRepMat)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: isInsideEllipsoid
#endif
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: NormedPoint(nd)
        real(RK)   , intent(in) :: InvRepMat(nd,nd)
        logical                 :: isInsideEllipsoid
        isInsideEllipsoid = dot_product(NormedPoint,matmul(InvRepMat,NormedPoint)) <= 1._RK
    end function isInsideEllipsoid

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the natural logarithm of the probability density function value of a point uniformly distributed within an ellipsoid,
    !> whose logarithm of the square root of the determinant of its representative covariance matrix is given by `logSqrtDetCovMat`.
    !>
    !> \param[in]   nd                  :   The number of dimensions of the MVU distribution.
    !> \param[in]   logSqrtDetCovMat    :   The logarithm of the square root of the determinant of
    !>                                      the inverse of the representative covariance matrix of the ellipsoid.
    !>
    !> \return
    !> `logProbMVU` :   The natural logarithm of the probability density function value
    !>                  of a point uniformly distributed within the target ellipsoid.
    !>
    !> \author
    !> Amir Shahmoradi, April 23, 2017, 1:36 AM, ICES, UTEXAS
    pure function getLogProbMVU(nd,logSqrtDetCovMat) result(logProbMVU)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMVU
#endif
        use Math_mod, only: getLogVolEllipsoid
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: logSqrtDetCovMat
        real(RK)                :: logProbMVU
        logProbMVU = -getLogVolEllipsoid(nd,logSqrtDetCovMat)
    end function getLogProbMVU

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return a random point on the target ellipsoid by projecting a random point uniformly distributed within the ellipsoid on its boundary.
    !>
    !> \param[in]   nd              :   The number of dimensions of the ellipsoid.
    !> \param[in]   MeanVec         :   The mean vector (center) of the ellipsoid.
    !> \param[in]   CholeskyLower   :   The Cholesky lower triangle of the representative covariance matrix of the ellipsoid.
    !> \param[in]   Diagonal        :   The Diagonal elements of the Cholesky lower triangle of the representative covariance matrix of the ellipsoid.
    !>
    !> \return
    !> `RandPointOnEllipsoid` : A random point on the target ellipsoid by projecting a random
    !>                          point uniformly distributed within the ellipsoid on its boundary.
    !>
    !> \remark
    !> This is algorithm is similar to [getRandMVU](@ref getrandmvu), with the only difference that
    !> points are drawn randomly from the surface of the ellipsoid instead of inside of its interior.
    !>
    !> \remark
    !> Note that the distribution of points on the surface of the ellipsoid is **NOT** uniform.
    !> Regions of high curvature will have more points randomly sampled from them.
    !> Generating uniform random points on arbitrary-dimension ellipsoids is not a task with trivial solution!
    !>
    !> \todo
    !> The performance of this algorithm can be further improved.
    !>
    !> \author
    !> Amir Shahmoradi, April 23, 2017, 1:36 AM, ICES, UTEXAS
    function getRandPointOnEllipsoid(nd,MeanVec,CholeskyLower,Diagonal) result(RandPointOnEllipsoid)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandPointOnEllipsoid
#endif
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: MeanVec(nd)
        real(RK)   , intent(in) :: CholeskyLower(nd,nd) ! Cholesky lower triangle, calculated from the MVN CovMat.
        real(RK)   , intent(in) :: Diagonal(nd)         ! Cholesky diagonal terms, calculated from the MVN CovMat.
        real(RK)                :: RandPointOnEllipsoid(nd), DummyVec(nd), sumSqDummyVec
        integer(IK)             :: i,j
        sumSqDummyVec = 0._RK
        do j=1,nd
            DummyVec(j) = getRandGaus()
            sumSqDummyVec = sumSqDummyVec + DummyVec(j)**2
        end do
        DummyVec = DummyVec / sqrt(sumSqDummyVec)    ! DummyVec is a random point on the surface of nd-sphere.
        RandPointOnEllipsoid = 0._RK
        do j = 1,nd
            RandPointOnEllipsoid(j) = RandPointOnEllipsoid(j) + Diagonal(j) * DummyVec(j)
            do i = j+1,nd
                RandPointOnEllipsoid(i) = RandPointOnEllipsoid(i) + CholeskyLower(i,j) * DummyVec(j)
            end do
        end do
        RandPointOnEllipsoid = RandPointOnEllipsoid + MeanVec
    end function getRandPointOnEllipsoid

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the natural logarithm of the Lognormal probability density function.
    !>
    !> \param[in]   logMean                 :   The mean of the Lognormal distribution.
    !> \param[in]   inverseVariance         :   The inverse variance of the Lognormal distribution.
    !> \param[in]   logSqrtInverseVariance  :   The natural logarithm of the square root of the inverse variance of the Lognormal distribution.
    !> \param[in]   logPoint                :   The natural logarithm of the point at which the Lognormal PDF must be computed.
    !>
    !> \return
    !> `logProbLogn` : The natural logarithm of the Lognormal probability density function.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    pure function getLogProbLognSP(logMean,inverseVariance,logSqrtInverseVariance,logPoint) result(logProbLogn)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbLognSP
#endif
        use Constants_mod, only: LOGINVSQRT2PI
        implicit none
        real(RK), intent(in) :: logMean,inverseVariance,logSqrtInverseVariance,logPoint
        real(RK)             :: logProbLogn
        logProbLogn = LOGINVSQRT2PI + logSqrtInverseVariance - logPoint - 0.5_RK * inverseVariance * (logPoint-logMean)**2
    end function getLogProbLognSP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the natural logarithm of the Lognormal probability density function.
    !>
    !> \param[in]   np                      :   The size of the input vector of points represented by `LogPoint`.
    !> \param[in]   logMean                 :   The mean of the Lognormal distribution.
    !> \param[in]   inverseVariance         :   The inverse variance of the Lognormal distribution.
    !> \param[in]   logSqrtInverseVariance  :   The natural logarithm of the square root of the inverse variance of the Lognormal distribution.
    !> \param[in]   LogPoint                :   The natural logarithm of the vector of points at which the Lognormal PDF must be computed.
    !>
    !> \return
    !> `logProbLogn` : The natural logarithm of the Lognormal probability density function.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    pure function getLogProbLognMP(np,logMean,inverseVariance,logSqrtInverseVariance,LogPoint) result(logProbLogn)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbLognMP
#endif
        use Constants_mod, only: LOGINVSQRT2PI
        implicit none
        integer(IK), intent(in) :: np
        real(RK)   , intent(in) :: logMean,inverseVariance,logSqrtInverseVariance,LogPoint(np)
        real(RK)                :: logProbLogn(np)
        logProbLogn = LOGINVSQRT2PI + logSqrtInverseVariance - LogPoint - 0.5_RK * inverseVariance * (LogPoint-logMean)**2
    end function getLogProbLognMP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return a single-precision uniformly-distributed random real-valued number in the range `[0,1]`.
    !>
    !> \param[inout]    idum    :   The input integer random seed with the `save` attribute.
    !>
    !> \return
    !> `randRealLecuyer`        :   A single-precision uniformly-distributed random real-valued number in the range `[0,1]`.
    !>
    !> \warning
    !> Do not change the value of `idum` between calls.
    !>
    !> \remark
    !> This routine is guaranteed to random numbers with priodicity larger than `~2*10**18` random numbers.
    !> For more info see Press et al. (1990) Numerical Recipes.
    !>
    !> \remark
    !> This routine is solely kept for backward compatibility reasons.
    !> The Fortran intrinsic subroutine `random_number()` is recommended to be used against this function.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    function getRandRealLecuyer(idum) result(randRealLecuyer)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandRealLecuyer
#endif
        implicit none
        integer(IK), intent(inout) :: idum
        integer(IK), parameter     :: im1=2147483563, im2=2147483399, imm1=im1-1, ia1=40014, ia2=40692
        integer(IK), parameter     :: iq1=53668, iq2=52774, ir1=12211, ir2=3791, ntab=32, ndiv=1+imm1/ntab
        real(RK)   , parameter     :: am=1._RK/real(im1,kind=RK), eps=1.2e-7_RK, rnmx=1._RK-eps
        real(RK)                   :: randRealLecuyer
        integer(IK)                :: idum2,j,k,iv(ntab),iy
        save                       :: iv, iy, idum2
        data idum2/123456789/, iv/ntab*0/, iy/0/
        if (idum <= 0) then
            idum = max(-idum,1)
            idum2 = idum
            do j = ntab+8,1,-1
                k = idum/iq1
                idum = ia1*(idum-k*iq1)-k*ir1
                if (idum < 0) idum = idum+im1
                if (j <= ntab) iv(j) = idum
            end do
            iy = iv(1)
        endif
        k = idum/iq1
        idum = ia1*(idum-k*iq1)-k*ir1
        if (idum < 0) idum=idum+im1
        k = idum2/iq2
        idum2 = ia2*(idum2-k*iq2)-k*ir2
        if (idum2 < 0) idum2=idum2+im2
        j = 1+iy/ndiv
        iy = iv(j)-idum2
        iv(j) = idum
        if(iy < 1)iy = iy+imm1
        randRealLecuyer = min(am*iy,rnmx)
    end function getRandRealLecuyer

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return an integer uniformly-distributed random integer-valued number in the range `[lowerBound , upperBound]`.
    !>
    !> \param[in]       lowerBound  :   The inclusive integer lower bound of the integer uniform distribution.
    !> \param[in]       upperBound  :   The inclusive integer upper bound of the integer uniform distribution.
    !> \param[inout]    idum        :   The input integer random seed with the `save` attribute.
    !>
    !> \return
    !> `randRealLecuyer`    : A uniformly-distributed random integer-valued number in the range `[lowerBound , upperBound]`.
    !>
    !> \warning
    !> Do not change the value of `idum` between calls.
    !>
    !> \remark
    !> This routine is guaranteed to random numbers with priodicity larger than `~2*10**18` random numbers.
    !> For more info see Press et al. (1990) Numerical Recipes.
    !>
    !> \remark
    !> This routine is solely kept for backward compatibility reasons.
    !> The [getRandInt](@ref getrandint) is recommended to be used against this routine.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    function getRandIntLecuyer(lowerBound,upperBound,idum)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandIntLecuyer
#endif
        implicit none
        integer(IK), intent(in)    :: lowerBound,upperBound
        integer(IK), intent(inout) :: idum
        integer(IK)                :: getRandIntLecuyer
        getRandIntLecuyer = lowerBound + nint( getRandRealLecuyer(idum)*real(upperBound-lowerBound,kind=RK) )
    end function getRandIntLecuyer

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return an integer uniformly-distributed random integer-valued number in the range `[lowerBound , upperBound]`.
    !>
    !> \param[in]       lowerBound  :   The lower bound of the integer uniform distribution.
    !> \param[in]       upperBound  :   The upper bound of the integer uniform distribution.
    !>
    !> \return
    !> `randInt` : A uniformly-distributed random integer-valued number in the range `[lowerBound , upperBound]`.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    function getRandInt(lowerBound,upperBound) result(randInt)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandInt
#endif
        implicit none
        integer(IK), intent(in) :: lowerBound,upperBound
        real(RK)                :: dummy
        integer(IK)             :: randInt
        call random_number(dummy)
        randInt = lowerBound + nint( dummy*real(upperBound-lowerBound,kind=RK) )
    end function getRandInt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return an integer uniformly-distributed random integer-valued number in the range `[lowerBound , upperBound]` using
    !> the built-in random number generator of Fortran.
    !>
    !> \param[in]       lowerBound  :   The lower bound of the integer uniform distribution.
    !> \param[in]       upperBound  :   The upper bound of the integer uniform distribution.
    !>
    !> \return
    !> `randUniform`    : A uniformly-distributed random real-valued number in the range `[lowerBound , upperBound]`.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    function getRandUniform(lowerBound,upperBound) result(randUniform)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandUniform
#endif
        implicit none
        real(RK), intent(in)    :: lowerBound, upperBound
        real(RK)                :: randUniform
        call random_number(randUniform)
        randUniform = lowerBound + randUniform * (upperBound - lowerBound)
    end function getRandUniform

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return a Gamma-distributed random number, following the prescription in the GSL library.
    !>
    !> \param[in]       alpha   :   The shape parameter of the Gamma distribution.
    !>
    !> \return
    !> `randGamma`    : A Gamma-distributed random real-valued number in the range `[0,+Infinity]`.
    !>
    !> \warning
    !> The condition `alpha > 0` must hold, otherwise the negative value `-1` will be returned to indicate the occurrence of error.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    function getRandGamma(alpha) result(randGamma)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandGamma
#endif
        implicit none
        real(RK), intent(in) :: alpha
        real(RK)             :: randGamma
        real(RK)             :: c,u,v,z
        if (alpha<=0._RK) then  ! illegal value of alpha
            randGamma = -1._RK
            return
        else
            randGamma = alpha
            if (randGamma<1._RK) randGamma = randGamma + 1._RK
            randGamma = randGamma - 0.3333333333333333_RK
            c = 3._RK*sqrt(randGamma)
            c = 1._RK / c
            do
                do
                    z = getRandGaus()
                    v = 1._RK + c*z
                    if (v<=0._RK) cycle
                    exit
                end do
                v = v**3
                call random_number(u)
                if ( log(u) >= 0.5_RK * z**2 + randGamma * ( 1._RK - v + log(v) ) ) cycle
                randGamma = randGamma * v
                exit
            end do
            if (alpha<1._RK) then
                call random_number(u)
                randGamma = randGamma * u**(1._RK/alpha)
            end if
        end if
    end function getRandGamma

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return a Gamma-distributed random number, whose shape parameter `alpha` is an integer.
    !>
    !> \param[in]       alpha   :   The shape integer parameter of the Gamma distribution.
    !>
    !> \return
    !> `randGamma`    : A Gamma-distributed random real-valued number in the range `[0,+Infinity]`.
    !>
    !> \warning
    !> The condition `alpha > 1` must hold, otherwise the negative value `-1` will be returned to indicate the occurrence of error.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    function getRandGammaIntShape(alpha) result(randGamma)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandGammaIntShape
#endif
        implicit none
        integer(IK), intent(in) :: alpha
        real(RK)                :: randGamma
        real(RK)                :: am,e,h,s,x,y,Vector(2),Array(5)
        if (alpha < 1) then  ! illegal value of alpha
            randGamma = -1._RK
            return
        elseif (alpha < 6) then
            call random_number(Array(1:alpha))
            x = -log(product(Array(1:alpha)))
        else    ! use rejection sampling
            do
                call random_number(Vector)
                Vector(2) = 2._RK*Vector(2)-1._RK
                if (dot_product(Vector,Vector) > 1._RK) cycle
                y = Vector(2) / Vector(1)
                am = alpha - 1
                s = sqrt(2._RK*am + 1._RK)
                x = s*y + am
                if (x <= 0.0) cycle
                e = (1._RK+y**2) * exp(am*log(x/am)-s*y)
                call random_number(h)
                if (h <= e) exit
            end do
        end if
        randGamma = x
    end function getRandGammaIntShape

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return a random Beta-distributed variable.
    !>
    !> \param[in]       alpha   :   The first shape parameter of the Beta distribution.
    !> \param[in]       beta    :   The second shape parameter of the Beta distribution.
    !>
    !> \return
    !> `randBeta`   : A Beta-distributed random real-valued number in the range `[0,1]`.
    !>
    !> \warning
    !> The conditions `alpha > 0` and `beta > 0` must hold, otherwise the negative
    !> value `-1` will be returned to indicate the occurrence of error.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    function getRandBeta(alpha,beta) result(randBeta)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandBeta
#endif
        implicit none
        real(RK), intent(in) :: alpha,beta
        real(RK)             :: randBeta
        real(RK)             :: x
        if ( alpha>0._RK .and. beta>0._RK ) then
            x = getRandGamma(alpha)
            randBeta = x / ( x + getRandGamma(beta) )
        else ! illegal value of alpha or beta
            randBeta = -1._RK
        end if
    end function getRandBeta

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return a random Exponential-distributed value whose inverse mean is given as input.
    !>
    !> \param[in]   invMean :   The inverse of the mean of the Exponential distribution.
    !>
    !> \return
    !> `randExp` : An Exponential-distributed random real-valued number in the range `[0,+Infinity]` with mean `1 / invMean`.
    !>
    !> \warning
    !> It is the user's onus to ensure `invMean > 0` on input. This condition will NOT be checked by this routine.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    function getRandExpWithInvMean(invMean) result(randExp)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandExpWithInvMean
#endif
        implicit none
        real(RK), intent(in)    :: invMean
        real(RK)                :: randExp
        call random_number(randExp)
        randExp = -log(randExp) * invMean
    end function getRandExpWithInvMean

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return a random Exponential-distributed value whose mean is \f$\lambda = 1\f$.
    !>
    !> \return
    !> `randExp` : A random Exponential-distributed value whose mean \f$\lambda = 1\f$.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    function getRandExp() result(randExp)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandExp
#endif
        implicit none
        real(RK) :: randExp
        call random_number(randExp)
        randExp = -log(randExp)
    end function getRandExp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return a random correlation matrix.
    !>
    !> \param[in]   nd  :   The rank of the correlation matrix.
    !> \param[in]   eta :   The parameter that roughly represents the shape parameters of the Beta distribution.
    !>                      The larger the value of `eta` is, the more homogeneous the correlation matrix will look.
    !>                      In other words, set this parameter to some small number to generate strong random correlations
    !>                      in the output random correlation matrix.
    !>
    !> \return
    !> `RandCorMat` : A random correlation matrix.
    !>
    !> \warning
    !> The conditions `nd > 1` and `eta > 0.0` must hold, otherwise the first element of
    !> output, `getRandCorMat(1,1)`, will be set to `-1` to indicate the occurrence of an error.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    function getRandCorMat(nd,eta) result(RandCorMat)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandCorMat
#endif
        use Matrix_mod, only: getCholeskyFactor
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: eta
        real(RK)                :: RandCorMat(nd,nd), dummy
        real(RK)                :: beta,sumSqDummyVec,DummyVec(nd-1),W(nd-1),Diagonal(nd-1)
        integer(IK)             :: m,j,i

        if (nd<2_IK .or. eta<=0._RK) then  ! illegal value for eta.
            RandCorMat(1,1) = -1._RK
            return
        end if

        do m = 1,nd
            RandCorMat(m,m) = 1._RK
        end do
        beta = eta + 0.5_RK*(nd-2._RK)
        dummy = getRandBeta(beta,beta)
        if (dummy<=0._RK .or. dummy>=1._RK) then
        ! LCOV_EXCL_START
            error stop
            !call abortProgram( output_unit , 1 , 1 , 'Statitistics@getRandCorMat() failed. Random Beta variable out of bound: ' // num2str(dummy) )
        end if
        ! LCOV_EXCL_STOP
        RandCorMat(1,2) = 2._RK * dummy - 1._RK ! for the moment, only the upper half of RandCorMat is needed, the lower half will contain cholesky lower triangle.

        do m = 2,nd-1
            beta = beta - 0.5_RK
            sumSqDummyVec = 0._RK
            do j=1,m
                DummyVec(j) = getRandGaus()
                sumSqDummyVec = sumSqDummyVec + DummyVec(j)**2
            end do
            DummyVec(1:m) = DummyVec(1:m) / sqrt(sumSqDummyVec)   ! DummyVec is now a uniform random point from inside of m-sphere
            dummy = getRandBeta(0.5e0_RK*m,beta)
            W(1:m) = sqrt(dummy) * DummyVec(1:m)
            call getCholeskyFactor(m,RandCorMat(1:m,1:m),Diagonal(1:m))
            if (Diagonal(1)<0._RK) then
                error stop
                !call abortProgram( output_unit , 1 , 1 , 'Statitistics@getRandCorMat()@getCholeskyFactor() failed.' )
            end if
            DummyVec(1:m) = 0._RK
            do j = 1,m
                DummyVec(j) = DummyVec(j) + Diagonal(j) * W(j)
                do i = j+1,m
                    DummyVec(i) = DummyVec(i) + RandCorMat(i,j) * DummyVec(j)
                end do
            end do
            RandCorMat(1:m,m+1) = DummyVec(1:m)
        end do
        do i=1,nd-1
            RandCorMat(i+1:nd,i) = RandCorMat(i,i+1:nd)
        end do
  end function getRandCorMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!  function getRandCorMat(nd,eta)    ! based on the idea of LKJ (2007). But there is something wrong with this routine
!    use Matrix_mod, only: getCholeskyFactor
!    implicit none
!    !integer, intent(in) :: nd,eta
!    integer, intent(in) :: nd
!    real(RK), intent(in) :: eta
!    integer :: m,mNew,j,i
!    real(RK) :: getRandCorMat(nd,nd), dummy, failureCounter
!    real(RK) :: beta,sumSqDummyVec,DummyVec(nd-1),W(nd-1),Diagonal(nd-1)
!
!    if (nd<2 .or. eta<=0._RK) then  ! illegal value for eta. set getRandCorMat=0, return
!      getRandCorMat = -1._RK
!      return
!    end if
!
!    do m = 1,nd
!      getRandCorMat(m,m) = 1._RK
!    end do
!    beta = eta + 0.5_RK*(nd-2._RK)
!
!    do
!      dummy = getRandBeta(beta,beta)
!      if (dummy>0._RK .and. dummy<1._RK) exit
!      write(*,*) "**Warning** random Beta variable out of bound.", dummy
!      write(*,*) "Something is wrong with getRandBeta()."
!      cycle
!    end do
!    getRandCorMat(1,2) = 2._RK * dummy - 1._RK    ! for the moment, only the upper half of getRandCorMat is needed, the lower half will contain cholesky lower triangle.
!
!    m = 2
!    call getCholeskyFactor(m,getRandCorMat(1:m,1:m),Diagonal(1:m))
!
!    failureCounter = 0
!    onionLayer: do
!
!      beta = beta - 0.5_RK
!
!      sumSqDummyVec = 0._RK
!      do j=1,m
!        DummyVec(j) = getRandGaus()
!        sumSqDummyVec = sumSqDummyVec + DummyVec(j)**2
!      end do
!      DummyVec(1:m) = DummyVec(1:m) / sqrt(sumSqDummyVec)   ! DummyVec is now a uniform random point from inside of m-sphere
!
!      mNew = m + 1
!      posDefCheck: do
!
!        do
!          dummy = getRandBeta(0.5_RK*m,beta)
!          if (dummy>0._RK .and. dummy<1._RK) exit
!          write(*,*) "**Warning** random Beta variable out of bound.", dummy
!          write(*,*) "Something is wrong with getRandBeta()."
!          read(*,*)
!          cycle
!        end do
!        W(1:m) = sqrt(dummy) * DummyVec(1:m)
!
!        getRandCorMat(1:m,mNew) = 0._RK
!        do j = 1,m
!          getRandCorMat(j,mNew) = getRandCorMat(j,mNew) + Diagonal(j) * W(j)
!          do i = j+1,m
!            getRandCorMat(i,mNew) = getRandCorMat(i,mNew) + getRandCorMat(i,j) * getRandCorMat(j,mNew)
!          end do
!        end do
!
!
!        call getCholeskyFactor(mNew,getRandCorMat(1:mNew,1:mNew),Diagonal(1:mNew))  ! Now check if the new matrix is positive-definite, then proceed with the next layer
!        if (Diagonal(1)<0._RK) then
!          failureCounter = failureCounter + 1
!          cycle posDefCheck
!          !write(*,*) "Cholesky factorization failed in getRandCorMat()."
!          !write(*,*) m
!          !write(*,*) getRandCorMat(1:m,1:m)
!          !stop
!        end if
!        exit posDefCheck
!
!      end do posDefCheck
!
!      if (mNew==nd) exit onionLayer
!      m = mNew
!
!    end do onionLayer
!
!    if (failureCounter>0) write(*,*) 'failureRatio: ', dble(failureCounter)/dble(nd-2)
!    do i=1,nd-1
!      getRandCorMat(i+1:nd,i) = getRandCorMat(i,i+1:nd)
!    end do
!
!  end function getRandCorMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Returns a random correlation matrix using Monte Carlo rejection method.
    !>
    !> \param[in]   nd      :   The rank of the correlation matrix.
    !> \param[in]   minRho  :   The minimum correlation coefficient to be expected in the output random correlation matrix.
    !> \param[in]   maxRho  :   The maximum correlation coefficient to be expected in the output random correlation matrix.
    !>
    !> \return
    !> `RandCorMat` : A random correlation matrix. Only the upper half of
    !> `RandCorMat` is the correlation matrix, lower half is NOT set on output.
    !>
    !> \warning
    !> The conditions `nd >= 1` and `maxRho < minRho` must hold, otherwise, `RandCorMat(1,1) = -1._RK` will be returned.
    !>
    !> \remark
    !> This subroutine is very slow for high matrix dimensions ( `nd >~ 10` ).
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    function getRandCorMatRejection(nd,minRho,maxRho) result(RandCorMat)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandCorMatRejection
#endif

        use Matrix_mod, only: isPosDef
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: minRho,maxRho
        real(RK)                :: RandCorMat(nd,nd), RhoVec(nd*(nd-1))
        integer(IK)             :: i,j,irho
        if (maxRho<minRho .or. nd<1_IK) then
            RandCorMat(1,1) = -1._RK
            return
        end if
        if (nd==1_IK) then
          RandCorMat = 1._RK
        else
            rejection: do
                call random_number(RhoVec)
                RhoVec = minRho + RhoVec*(maxRho-minRho)
                irho = 0
                do j=1,nd
                    RandCorMat(j,j) = 1._RK
                    do i=1,j-1
                        irho = irho + 1
                        RandCorMat(i,j) = RhoVec(irho)
                    end do
                end do
                if (isPosDef(nd,RandCorMat)) exit rejection
                cycle rejection
            end do rejection
        end if
        do j=1,nd-1
            RandCorMat(j+1:nd,j) = RandCorMat(j,j+1:nd)
        end do
  end function getRandCorMatRejection

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Convert the upper-triangle covariance matrix to the upper-triangle correlation matrix.
    !>
    !> \param[in]   nd          :   The rank of the covariance matrix.
    !> \param[in]   CovMatUpper :   The upper-triangle covariance matrix. The lower-triangle will not be used.
    !>
    !> \return
    !> `CorMatUpper` : An upper-triangle correlation matrix. The lower-triangle will NOT be set.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    pure function getCorMatUpperFromCovMatUpper(nd,CovMatUpper) result(CorMatUpper)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCorMatUpperFromCovMatUpper
#endif
        implicit none
        integer(IK)  , intent(in) :: nd
        real(RK)     , intent(in) :: CovMatUpper(nd,nd)
        real(RK)                  :: CorMatUpper(nd,nd)
        real(RK)                  :: InverseStdVec(nd)
        integer(IK)               :: i,j
        do j = 1, nd
            InverseStdVec(j) = 1._RK / sqrt(CovMatUpper(j,j))
            do i = 1, j
                CorMatUpper(i,j) = CovMatUpper(i,j) * InverseStdVec(j) * InverseStdVec(i)
            end do
        end do
    end function getCorMatUpperFromCovMatUpper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Convert the upper-triangle correlation matrix to the upper-triangle covariance matrix.
    !>
    !> \param[in]   nd          :   The rank of the correlation matrix.
    !> \param[in]   StdVec      :   The input standard deviation vector of length `nd`.
    !> \param[in]   CorMatUpper :   The upper-triangle correlation matrix. The lower-triangle will not be used.
    !>
    !> \return
    !> `CovMatUpper` : An upper-triangle covariance matrix. The lower-triangle will NOT be set.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    pure function getCovMatUpperFromCorMatUpper(nd,StdVec,CorMatUpper) result(CovMatUpper)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovMatUpperFromCorMatUpper
#endif
        implicit none
        integer(IK)  , intent(in) :: nd
        real(RK)     , intent(in) :: StdVec(nd), CorMatUpper(nd,nd)
        real(RK)                  :: CovMatUpper(nd,nd)
        integer(IK)               :: i,j
        do j=1,nd
            do i=1,j
                CovMatUpper(i,j) = CorMatUpper(i,j) * StdVec(j) * StdVec(i)
            end do
        end do
    end function getCovMatUpperFromCorMatUpper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Convert the lower-triangle correlation matrix to the upper-triangle covariance matrix.
    !>
    !> \param[in]   nd          :   The rank of the correlation matrix.
    !> \param[in]   StdVec      :   The input standard deviation vector of length `nd`.
    !> \param[in]   CorMatLower :   The lower-triangle correlation matrix. The upper-triangle will not be used.
    !>
    !> \return
    !> `CovMatUpper` : An upper-triangle covariance matrix. The lower-triangle will NOT be set.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    pure function getCovMatUpperFromCorMatLower(nd,StdVec,CorMatLower) result(CovMatUpper)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovMatUpperFromCorMatLower
#endif
        implicit none
        integer(IK)  , intent(in) :: nd
        real(RK)     , intent(in) :: StdVec(nd), CorMatLower(nd,nd)
        real(RK)                  :: CovMatUpper(nd,nd)
        integer(IK)               :: i,j
        do j=1,nd
            CovMatUpper(j,j) = StdVec(j)**2
            do i=1,j-1
                CovMatUpper(i,j) = CorMatLower(j,i) * StdVec(j) * StdVec(i)
            end do
        end do
    end function getCovMatUpperFromCorMatLower

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Convert the upper-triangle correlation matrix to the lower-triangle covariance matrix.
    !>
    !> \param[in]   nd          :   The rank of the correlation matrix.
    !> \param[in]   StdVec      :   The input standard deviation vector of length `nd`.
    !> \param[in]   CorMatUpper :   The upper-triangle correlation matrix. The lower-triangle will not be used.
    !>
    !> \return
    !> `CovMatLower` : An lower-triangle covariance matrix. The upper-triangle will NOT be set.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    pure function getCovMatLowerFromCorMatUpper(nd,StdVec,CorMatUpper) result(CovMatLower)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovMatLowerFromCorMatUpper
#endif
        implicit none
        integer(IK)  , intent(in) :: nd
        real(RK)     , intent(in) :: StdVec(nd), CorMatUpper(nd,nd)
        real(RK)                  :: CovMatLower(nd,nd)
        integer(IK)               :: i,j
        do j=1,nd
            CovMatLower(j,j) = StdVec(j)**2
            do i=1,j-1
                CovMatLower(j,i) = CorMatUpper(i,j) * StdVec(j) * StdVec(i)
            end do
        end do
    end function getCovMatLowerFromCorMatUpper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Convert the lower-triangle correlation matrix to the lower-triangle covariance matrix.
    !>
    !> \param[in]   nd          :   The rank of the correlation matrix.
    !> \param[in]   StdVec      :   The input standard deviation vector of length `nd`.
    !> \param[in]   CorMatLower :   The lower-triangle correlation matrix. The upper-triangle will not be used.
    !>
    !> \return
    !> `CovMatLower` : An lower-triangle covariance matrix. The upper-triangle will NOT be set.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    pure function getCovMatLowerFromCorMatLower(nd,StdVec,CorMatLower) result(CovMatLower)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovMatLowerFromCorMatLower
#endif
        implicit none
        integer(IK)  , intent(in) :: nd
        real(RK)     , intent(in) :: StdVec(nd), CorMatLower(nd,nd)
        real(RK)                  :: CovMatLower(nd,nd)
        integer(IK)               :: i,j
        do j=1,nd
            CovMatLower(j,j) = StdVec(j)**2
            do i=1,j-1
                CovMatLower(j,i) = CorMatLower(j,i) * StdVec(j) * StdVec(i)
            end do
        end do
    end function getCovMatLowerFromCorMatLower

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Convert the input correlation matrix to the output covariance matrix.
    !>
    !> \param[in]   nd          :   The rank of the correlation matrix.
    !> \param[in]   StdVec      :   The input standard deviation vector of length `nd`.
    !> \param[in]   CorMatUpper :   The upper-triangle correlation matrix. The lower-triangle will not be used.
    !>
    !> \return
    !> `CovMatFull` : The full covariance matrix.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    pure function getCovMatFromCorMatUpper(nd,StdVec,CorMatUpper) result(CovMatFull)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovMatFromCorMatUpper
#endif
        implicit none
        integer(IK)  , intent(in) :: nd
        real(RK)     , intent(in) :: StdVec(nd), CorMatUpper(nd,nd)   ! only upper half needed
        real(RK)                  :: CovMatFull(nd,nd)
        integer(IK)               :: i,j
        do j=1,nd
            CovMatFull(j,j) = StdVec(j)**2
            do i=1,j-1
                CovMatFull(i,j) = CorMatUpper(i,j) * StdVec(j) * StdVec(i)
                CovMatFull(j,i) = CovMatFull(i,j)
            end do
        end do
    end function getCovMatFromCorMatUpper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    !> \brief
!    !> Return the Geometric distribution PDF values for a range of trials, starting at index `1`.
!    !> If the probability of success on each trial is `successProb`, then the probability that
!    !> the `k`th trial (out of `k` trials) is the first success is `GeoLogPDF(k)`.
!    !>
!    !> \param[in]   successProb     :   The probability of success.
!    !> \param[in]   logPdfPrecision :   The precision value below which the PDF is practically considered to be zero (**optional**).
!    !> \param[in]   minSeqLen       :   The minimum length of the range of `k` values for which the PDF will be computed (**optional**).
!    !> \param[in]   seqLen          :   The length of the range of `k` values for which the PDF will be computed (**optional**).
!    !>                                  If provided, it will overwrite the the output sequence length as inferred from
!    !>                                  the combination of `minSeqLen` and `logPdfPrecision`.
!    !>
!    !> \return
!    !> `GeoLogPDF`  :   An allocatable representing the geometric PDF over a range of `k` values, whose length is
!    !>                  `seqLen`, or if not provided, is determined from the values of `logPdfPrecision` and `minSeqLen`.
!    !>
!    !> \author
!    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
!    function getGeoLogPDF_old(successProb,logPdfPrecision,minSeqLen,seqLen) result(GeoLogPDF)
!#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
!        !DEC$ ATTRIBUTES DLLEXPORT :: getGeoLogPDF_old
!#endif
!        use Constants_mod, only: IK, RK
!        implicit none
!        real(RK)    , intent(in)            :: successProb
!        real(RK)    , intent(in), optional  :: logPdfPrecision
!        integer(IK) , intent(in), optional  :: minSeqLen
!        integer(IK) , intent(in), optional  :: seqLen
!        real(RK)    , allocatable           :: GeoLogPDF(:)
!        real(RK)    , parameter             :: LOG_PDF_PRECISION = log(0.001_RK)
!        real(RK)                            :: logProbFailure
!        integer(IK)                         :: lenGeoLogPDF, i
!        logProbFailure = log(1._RK - successProb)
!        if (present(seqLen)) then
!            lenGeoLogPDF = seqLen
!        else
!            if (present(logPdfPrecision)) then
!                lenGeoLogPDF = ceiling(  logPdfPrecision / logProbFailure)
!            else
!                lenGeoLogPDF = ceiling(LOG_PDF_PRECISION / logProbFailure)
!            end if
!            if (present(minSeqLen)) lenGeoLogPDF = max(minSeqLen,lenGeoLogPDF)
!        end if
!        allocate(GeoLogPDF(lenGeoLogPDF))
!        GeoLogPDF(1) = log(successProb)
!        do i = 2, lenGeoLogPDF
!            GeoLogPDF(i) = GeoLogPDF(i-1) + logProbFailure
!        end do
!    end function getGeoLogPDF_old

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the Geometric distribution PDF values for the input number of trials,
    !> the trials at which first success happens, and the success probability.
    !>
    !> \param[in]   numTrial    :   The number of trials.
    !> \param[in]   SuccessStep :   The vector of trials of length `numTrial` at which the first success happens.
    !> \param[in]   successProb :   The probability of success.
    !>
    !> \return
    !> `LogProbGeo` :   An output vector of length `numTrial` representing the geometric PDF values.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    pure function getLogProbGeo(numTrial, SuccessStep, successProb) result(LogProbGeo)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbGeo
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK) , intent(in)    :: numTrial
        integer(IK) , intent(in)    :: SuccessStep(numTrial)
        real(RK)    , intent(in)    :: successProb
        real(RK)                    :: LogProbGeo(numTrial)
        real(RK)                    :: logProbSuccess, logProbFailure
        logProbSuccess = log(successProb)
        logProbFailure = log(1._RK - successProb)
        LogProbGeo = logProbSuccess + (SuccessStep-1_IK) * logProbFailure
    end function getLogProbGeo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Compute the natural logarithm of the Geometric distribution PDF of a limited range of Bernoulli trials,
    !> starting at index `1` up to `maxNumTrial`. In other words, upon reaching the trial `maxNumTrial`,
    !> the Bernoulli trials count restart from index `1`. This Cyclic Geometric distribution is
    !> particularly useful in the parallelization studies of Monte Carlo simulation.
    !>
    !> \param[in]   successProb :   The probability of success.
    !> \param[in]   maxNumTrial :   The maximum number of trails possible.
    !>                              After `maxNumTrial` tries, the Geometric distribution restarts from index `1`.
    !> \param[in]   numTrial    :   The length of the array `SuccessStep`. Note that `numTrial < maxNumTrial` must hold.
    !> \param[in]   SuccessStep :   A vector of length `(1:numTrial)` of integers that represent
    !>                              the steps at which the Bernoulli successes occur.
    !>
    !> \return
    !> `LogProbGeoCyclic`       :   A real-valued vector of length `(1:numTrial)` whose values represent the probabilities
    !>                              of having Bernoulli successes at the corresponding SuccessStep values.
    !>
    !> \warning
    !> Any value of SuccessStep must be an integer numbers between `1` and `maxNumTrial`.
    !> The onus is on the user to ensure this condition holds.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    pure function getLogProbGeoCyclic(successProb,maxNumTrial,numTrial,SuccessStep) result(LogProbGeoCyclic)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbGeoCyclic
#endif
        use Constants_mod, only: IK, RK, NEGLOGINF_RK
        implicit none
        real(RK)    , intent(in)    :: successProb
        integer(IK) , intent(in)    :: maxNumTrial
        integer(IK) , intent(in)    :: numTrial
        integer(IK) , intent(in)    :: SuccessStep(numTrial)
        real(RK)                    :: LogProbGeoCyclic(numTrial)
        real(RK)                    :: failureProb, logProbSuccess, logProbFailure, logDenominator, exponentiation
        if (successProb>0._RK .and. successProb<1._RK) then ! tolerate log(0)
            failureProb = 1._RK - successProb
            logProbSuccess = log(successProb)
            logProbFailure = log(failureProb)
            exponentiation = maxNumTrial * logProbFailure
            if (exponentiation<NEGLOGINF_RK) then ! tolerate underflow
                logDenominator = 0._RK
            else
                exponentiation = exp(exponentiation)
                if (exponentiation<1._RK) then ! tolerate log(0)
                    logDenominator = log(1._RK-exponentiation)
                else
                    logDenominator = NEGLOGINF_RK
                end if
            end if
            LogProbGeoCyclic = logProbSuccess + (SuccessStep-1) * logProbFailure - logDenominator
        elseif (successProb==0._RK) then
            LogProbGeoCyclic = -log(real(maxNumTrial,kind=RK))
        elseif (successProb==1._RK) then
            LogProbGeoCyclic(1) = 0._RK
            LogProbGeoCyclic(2:numTrial) = NEGLOGINF_RK
        else
            LogProbGeoCyclic = NEGLOGINF_RK
        end if
    end function getLogProbGeoCyclic

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the standard normal distribution PDF value.
    !>
    !> \param[in]   z   :   The input value at which the PDF will be computed.
    !>
    !> \return
    !> `snormPDF`   :   The standard normal distribution PDF value.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    function getSNormPDF(z) result(snormPDF)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSNormPDF
#endif
        use Constants_mod, only: INVSQRT2PI
        implicit none
        real(RK), intent(in) :: z
        real(RK)             :: snormPDF
        snormPDF = INVSQRT2PI * exp( -0.5_RK*z**2 )
    end function getSNormPDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the non-standard normal distribution PDF value.
    !>
    !> \param[in]   avg :   The mean of the Normal distribution.
    !> \param[in]   std :   The standard deviation of the Normal distribution.
    !> \param[in]   var :   The variance of the Normal distribution.
    !> \param[in]   x   :   The point at which the PDF will be computed.
    !>
    !> \return
    !> `normPDF`        :   The normal distribution PDF value at the given input point.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    function getNormPDF(avg,std,var,x) result(normPDF)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormPDF
#endif
        use Constants_mod, only: INVSQRT2PI
        implicit none
        real(RK), intent(in) :: avg,std,var,x
        real(RK)             :: normPDF
        normPDF = INVSQRT2PI * exp( -(x-avg)**2/(2._RK*var) ) / std
    end function getNormPDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the non-standard normal distribution Cumulative Probability Density function (CDF) value.
    !>
    !> \param[in]   avg :   The mean of the Normal distribution.
    !> \param[in]   std :   The standard deviation of the Normal distribution.
    !> \param[in]   x   :   The point at which the PDF will be computed.
    !>
    !> \return
    !> `normCDF`        :   The normal distribution CDF value at the given input point.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    pure function getNormCDF(avg,std,x) result(normCDF)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormCDF
#endif
        use Constants_mod, only: RK, SQRT2
        implicit none
        real(RK), intent(in) :: avg,std,x
        real(RK)             :: normCDF
        normCDF = 0.5_RK * ( 1._RK + erf( real((x-avg)/(SQRT2*std),kind=RK) ) )
    end function getNormCDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the standard normal distribution Cumulative Probability Density function (CDF) value.
    !>
    !> \param[in]   x       :   The point at which the PDF will be computed.
    !>
    !> \return
    !> `normCDF`   :   The normal distribution CDF value at the given input point.
    !>
    !> \remark
    !> This procedure performs all calculations in `real32` real kind. If 64-bit accuracy matters more than performance,
    !> then use the [getSNormCDF_DPR](@ref getsnormcdf_dpr) for a more-accurate double-precision but slower results.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    pure function getSNormCDF_SPR(x) result(normCDF)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSNormCDF_SPR
#endif
        use iso_fortran_env, only: RK => real32
        implicit none
        real(RK)    , intent(in)    :: x
        real(RK)                    :: normCDF
        real(RK), parameter         :: INVSQRT2 = 1._RK / sqrt(2._RK)
        normCDF = 0.5_RK * ( 1._RK + erf( real(x*INVSQRT2,kind=RK) ) )
    end function getSNormCDF_SPR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the standard normal distribution Cumulative Probability Density function (CDF) value.
    !>
    !> \param[in]   x   :   The point at which the PDF will be computed.
    !>
    !> \return
    !> `normCDF`   :   The normal distribution CDF value at the given input point.
    !>
    !> \remark
    !> This procedure performs all calculations in `DPR` (`real64`) real kind. If performance matters more than 64-bit accuracy,
    !> then use the [getSNormCDF_SPR](@ref getsnormcdf_spr) for a faster, but less-accurate single-precision results.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    pure function getSNormCDF_DPR(x) result(normCDF)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getSNormCDF_DPR
#endif
        use iso_fortran_env, only: RK => real64
        implicit none
        real(RK), intent(in)    :: x
        real(RK)                :: normCDF
        real(RK), parameter     :: INVSQRT2 = 1._RK / sqrt(2._RK)
        normCDF = 0.5_RK * ( 1._RK + erf( real(x*INVSQRT2,kind=RK) ) )
    end function getSNormCDF_DPR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the Beta distribution Cumulative Probability Density function (CDF) value.
    !>
    !> \param[in]   alpha   :   The first shape parameter of the Beta distribution.
    !> \param[in]   beta    :   The second shape parameter of the Beta distribution.
    !> \param[in]   x       :   The point at which the CDF will be computed.
    !>
    !> \return
    !> `betaCDF`   :   The Beta distribution CDF value at the given input point.
    !>
    !> \warning
    !> If `x` is not in the range `[0,1]`, a negative value for `betaCDF` will be returned to indicate the occurrence of error.
    !>
    !> \warning
    !> The onus is on the user to ensure that the input (`alpha`, `beta`) shape parameters are positive.
    !>
    !> \remark
    !> This procedure performs all calculations in `real32` real kind. If 64-bit accuracy matters more than performance,
    !> then use the [getBetaCDF_DPR](@ref getbetacdf_dpr) for a more-accurate double-precision but slower results.
    !>
    !> \todo
    !> The efficiency of this code can be improved by making `x` a vector on input.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    function getBetaCDF_SPR(alpha,beta,x) result(betaCDF)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaCDF_SPR
#endif
        use iso_fortran_env, only: RK => real32
        implicit none
        real(RK), intent(in)    :: alpha, beta, x
        real(RK)                :: bt
        real(RK)                :: betaCDF
        if (x < 0._RK .or. x > 1._RK) then
            betaCDF = -1._RK
            return
        end if
        if (x == 0._RK .or. x == 1._RK) then
            bt = 0._RK
        else
            bt = exp( log_gamma(real(alpha+beta,kind=RK)) &
                    - log_gamma(real(alpha,kind=RK)) - log_gamma(real(beta,kind=RK)) &
                    + alpha*log(x) + beta*log(1._RK-x) )
        end if
        if ( x < (alpha+1._RK) / (alpha+beta+2._RK) ) then
            betaCDF = bt * getBetaContinuedFraction(alpha,beta,x) / alpha
        else
            betaCDF = 1._RK - bt * getBetaContinuedFraction(beta,alpha,1._RK-x) / beta
        end if
    end function getBetaCDF_SPR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the Beta distribution Cumulative Probability Density function (CDF) value.
    !>
    !> \param[in]   alpha   :   The first shape parameter of the Beta distribution.
    !> \param[in]   beta    :   The second shape parameter of the Beta distribution.
    !> \param[in]   x       :   The point at which the CDF will be computed.
    !>
    !> \return
    !> `betaCDF`   :   The Beta distribution CDF value at the given input point.
    !>
    !> \warning
    !> If `x` is not in the range `[0,1]`, a negative value for `betaCDF` will be returned to indicate the occurrence of error.
    !>
    !> \warning
    !> The onus is on the user to ensure that the input (`alpha`, `beta`) shape parameters are positive.
    !>
    !> \remark
    !> This procedure performs all calculations in `DPR` (`real64`) real kind. If performance matters more than 64-bit accuracy,
    !> then use the [getBetaCDF_SPR](@ref getbetacdf_spr) for a faster, but less-accurate single-precision results.
    !>
    !> \todo
    !> The efficiency of this code can be improved by making `x` a vector on input.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    function getBetaCDF_DPR(alpha,beta,x) result(betaCDF)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaCDF_DPR
#endif
        use iso_fortran_env, only: RK => real64
        implicit none
        real(RK), intent(in)    :: alpha, beta, x
        real(RK)                :: bt
        real(RK)                :: betaCDF
        if (x < 0._RK .or. x > 1._RK) then
            betaCDF = -1._RK
            return
        end if
        if (x == 0._RK .or. x == 1._RK) then
            bt = 0._RK
        else
            bt = exp( log_gamma(real(alpha+beta,kind=RK)) &
                    - log_gamma(real(alpha,kind=RK)) - log_gamma(real(beta,kind=RK)) &
                    + alpha*log(x) + beta*log(1._RK-x) )
        end if
        if ( x < (alpha+1._RK) / (alpha+beta+2._RK) ) then
            betaCDF = bt * getBetaContinuedFraction(alpha,beta,x) / alpha
        else
            betaCDF = 1._RK - bt * getBetaContinuedFraction(beta,alpha,1._RK-x) / beta
        end if
    end function getBetaCDF_DPR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the Beta Continued Fraction (BCF).
    !>
    !> \param[in]   alpha   :   The first shape parameter of the Beta distribution.
    !> \param[in]   beta    :   The second shape parameter of the Beta distribution.
    !> \param[in]   x       :   The point at which the BCF will be computed.
    !>
    !> \return
    !> `betaCDF`   :   The BCF value at the given input point.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    function getBetaContinuedFraction_SPR(alpha,beta,x) result(betaContinuedFraction)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaContinuedFraction_SPR
#endif
        use iso_fortran_env, only: RK => real32
        implicit none
        real(RK)   , intent(in) :: alpha,beta,x
        real(RK)   , parameter  :: eps = epsilon(x), fpmin = tiny(x)/eps
        integer(IK), parameter  :: maxit = 100_IK
        real(RK)                :: aa,c,d,del,qab,qam,qap
        real(RK)                :: betaContinuedFraction
        integer(IK)             :: m,m2
        qab = alpha+beta
        qap = alpha+1._RK
        qam = alpha-1._RK
        c = 1._RK
        d = 1._RK-qab*x/qap
        if (abs(d) < fpmin) d = fpmin
        d = 1._RK/d
        betaContinuedFraction = d
        do m = 1,maxit
            m2 = 2*m
            aa = m*(beta-m)*x/((qam+m2)*(alpha+m2))
            d = 1._RK+aa*d
            if (abs(d) < fpmin) d = fpmin
            c = 1._RK+aa/c
            if (abs(c) < fpmin) c = fpmin
            d = 1._RK/d
            betaContinuedFraction = betaContinuedFraction*d*c
            aa = -(alpha+m)*(qab+m)*x/((alpha+m2)*(qap+m2))
            d = 1._RK+aa*d
            if (abs(d) < fpmin) d = fpmin
            c = 1._RK+aa/c
            if (abs(c) < fpmin) c = fpmin
            d = 1._RK/d
            del = d*c
            betaContinuedFraction = betaContinuedFraction*del
            if (abs(del-1._RK) <=  eps) exit
        end do
        if (m > maxit) then
        ! LCOV_EXCL_START
            error stop
            !call abortProgram( output_unit , 1 , 1 , &
            !'Statitistics@getBetaContinuedFraction_SPR() failed: alpha or beta too big, or maxit too small.' )
        end if
        ! LCOV_EXCL_STOP
    end function getBetaContinuedFraction_SPR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the Beta Continued Fraction (BCF).
    !>
    !> \param[in]   alpha   :   The first shape parameter of the Beta distribution.
    !> \param[in]   beta    :   The second shape parameter of the Beta distribution.
    !> \param[in]   x       :   The point at which the BCF will be computed.
    !>
    !> \return
    !> `betaCDF`   :   The BCF value at the given input point.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    function getBetaContinuedFraction_DPR(alpha,beta,x) result(betaContinuedFraction)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaContinuedFraction_DPR
#endif
        use iso_fortran_env, only: RK => real64
        implicit none
        real(RK)   , intent(in) :: alpha,beta,x
        real(RK)   , parameter  :: eps = epsilon(x), fpmin = tiny(x)/eps
        integer(IK), parameter  :: maxit = 100_IK
        real(RK)                :: aa,c,d,del,qab,qam,qap
        real(RK)                :: betaContinuedFraction
        integer(IK)             :: m,m2
        qab = alpha+beta
        qap = alpha+1._RK
        qam = alpha-1._RK
        c = 1._RK
        d = 1._RK-qab*x/qap
        if (abs(d) < fpmin) d = fpmin
        d = 1._RK/d
        betaContinuedFraction = d
        do m = 1,maxit
            m2 = 2*m
            aa = m*(beta-m)*x/((qam+m2)*(alpha+m2))
            d = 1._RK+aa*d
            if (abs(d) < fpmin) d = fpmin
            c = 1._RK+aa/c
            if (abs(c) < fpmin) c = fpmin
            d = 1._RK/d
            betaContinuedFraction = betaContinuedFraction*d*c
            aa = -(alpha+m)*(qab+m)*x/((alpha+m2)*(qap+m2))
            d = 1._RK+aa*d
            if (abs(d) < fpmin) d = fpmin
            c = 1._RK+aa/c
            if (abs(c) < fpmin) c = fpmin
            d = 1._RK/d
            del = d*c
            betaContinuedFraction = betaContinuedFraction*del
            if (abs(del-1._RK) <=  eps) exit
        end do
        if (m > maxit) then
        ! LCOV_EXCL_START
            error stop
            !call abortProgram( output_unit , 1 , 1 , &
            !'Statitistics@getBetaContinuedFraction_DPR() failed: alpha or beta too big, or maxit too small.' )
        end if
        ! LCOV_EXCL_STOP
    end function getBetaContinuedFraction_DPR

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the one-sample Kolmogorov–Smirnov (KS) test results for the null hypothesis that the data
    !> in vector `Point` comes from a distribution whose CDF is specified by the input `getCDF()` function.
    !>
    !> \param[in]       np      :   The number of points in the input vector `Point`.
    !> \param[inout]    Point   :   The sample. On return, this array will be sorted in Ascending order.
    !> \param[in]       getCDF  :   The function returning the Cumulative Distribution Function (CDF) of the sample.
    !> \param[out]      statKS  :   The KS test statistic.
    !> \param[out]      probKS  :   The `p`-value of the test, returned as a scalar value in the range `[0,1]`.
    !>                              The output `probKS` is the probability of observing a test statistic as extreme as,
    !>                              or more extreme than, the observed value under the null hypothesis.
    !>                              Small values of `probKS` cast doubt on the validity of the null hypothesis.
    !> \param[out]      Err     :   An object of class [Err_type](@ref err_mod::err_type) containing information
    !>                              about the occurrence of any error during the KS test computation.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    subroutine doKS1(np,Point,getCDF,statKS,probKS,Err)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: doKS1
#endif
        use Sort_mod, only : sortAscending
        use Err_mod, only: Err_type

        implicit none

        integer(IK)     , intent(in)    :: np
        real(RK)        , intent(out)   :: statKS,probKS
        real(RK)        , intent(inout) :: Point(np)
        type(Err_type)  , intent(out)   :: Err

        character(*)    , parameter     :: PROCEDURE_NAME = MODULE_NAME//"@doKS1"
        real(RK)                        :: npSqrt
        real(RK)                        :: cdf,cdfObserved,dt,frac
        integer(IK)                     :: j

        interface
            function getCDF(x)
                use Constants_mod, only: RK
                real(RK), intent(in) :: x
                real(RK)             :: getCDF
            end function getCDF
        end interface

        call sortAscending(np,Point,Err)
        if (Err%occurred) then
        ! LCOV_EXCL_START
            Err%msg = PROCEDURE_NAME//Err%msg
            return
        end if
        ! LCOV_EXCL_STOP

        statKS = 0._RK
        cdfObserved = 0._RK
        npSqrt = np
        do j = 1,np
            frac = j/npSqrt
            cdf = getCDF(Point(j))
            dt = max( abs(cdfObserved-cdf) , abs(frac-cdf) )
            if( dt > statKS ) statKS = dt
            cdfObserved = frac
        end do
        npSqrt = sqrt(npSqrt)
        probKS = getProbKS( (npSqrt+0.12_RK+0.11_RK/npSqrt)*statKS )

    end subroutine doKS1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the one-sample Kolmogorov–Smirnov (KS) test results for the assumption that the
    !> points originate from a uniform distribution in [0,1]. So, all Point must be in [0,1] on input.
    !>
    !> \param[in]       np      :   The number of points in the input vector `Point`.
    !> \param[inout]    Point   :   The sample. On return, this array will be sorted in Ascending order.
    !> \param[out]      statKS  :   The KS test statistic.
    !> \param[out]      probKS  :   The `p`-value of the test, returned as a scalar value in the range `[0,1]`.
    !>                              The output `probKS` is the probability of observing a test statistic as extreme as,
    !>                              or more extreme than, the observed value under the null hypothesis.
    !>                              Small values of `probKS` cast doubt on the validity of the null hypothesis.
    !> \param[out]      Err     :   An object of class [Err_type](@ref err_mod::err_type) containing information
    !>                              about the occurrence of any error during the KS test computation.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    pure subroutine doUniformKS1(np,Point,statKS,probKS,Err)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: doUniformKS1
#endif
        use Sort_mod, only : sortAscending
        use Err_mod, only: Err_type

        implicit none

        integer(IK)     , intent(in)    :: np
        real(RK)        , intent(out)   :: statKS,probKS
        real(RK)        , intent(inout) :: Point(np)
        type(Err_type)  , intent(out)   :: Err

        character(*)    , parameter     :: PROCEDURE_NAME = MODULE_NAME//"@doUniformKS1"
        real(RK)                        :: npSqrt
        real(RK)                        :: cdf,cdfObserved,dt,frac
        integer(IK)                     :: j

        call sortAscending(np,Point,Err)
        if (Err%occurred) then
            ! LCOV_EXCL_START
            Err%msg = PROCEDURE_NAME//Err%msg
            return
            ! LCOV_EXCL_STOP
        end if

        statKS = 0._RK
        cdfObserved = 0._RK
        npSqrt = np
        do j = 1,np
            frac = j/npSqrt
            cdf = Point(j)
            dt = max( abs(cdfObserved-cdf) , abs(frac-cdf) )
            if( dt > statKS ) statKS = dt
            cdfObserved = frac
        end do
        npSqrt = sqrt(npSqrt)
        probKS = getProbKS( (npSqrt+0.12_RK+0.11_RK/npSqrt)*statKS )
    end subroutine doUniformKS1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the two-sample Kolmogorov–Smirnov (KS) test results under the assumption that the
    !> points originate from the same distribution. It is assumed that the two input arrays are sorted ascending on input.
    !>
    !> \param[in]       np1             :   The number of points in the input vector `SortedPoint1`.
    !> \param[in]       np2             :   The number of points in the input vector `SortedPoint2`.
    !> \param[inout]    SortedPoint1    :   The first input sorted sample. On input, it must be sorted in ascending-order.
    !> \param[inout]    SortedPoint2    :   The second input sorted sample. On input, it must be sorted in ascending-order.
    !> \param[out]      statKS          :   The KS test statistic.
    !> \param[out]      probKS          :   The `p`-value of the test, returned as a scalar value in the range `[0,1]`.
    !>                                      The output `probKS` is the probability of observing a test statistic as extreme as,
    !>                                      or more extreme than, the observed value under the null hypothesis.
    !>                                      Small values of `probKS` cast doubt on the validity of the null hypothesis.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    subroutine doSortedKS2(np1,np2,SortedPoint1,SortedPoint2,statKS,probKS)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: doSortedKS2
#endif
        integer(IK) , intent(in)    :: np1, np2
        real(RK)    , intent(in)    :: SortedPoint1(np1), SortedPoint2(np2)
        real(RK)    , intent(out)   :: statKS, probKS
        real(RK)                    :: dummy1,dummy2,dt,np1_RK,np2_RK,npEffective,cdf1,cdf2
        integer(IK)                 :: j1,j2
        np1_RK = np1
        np2_RK = np2
        j1 = 1_IK
        j2 = 1_IK
        cdf1 = 0._RK
        cdf2 = 0._RK
        statKS = 0._RK
        do
            if ( j1<=np1 .and. j2<=np2 )then
                dummy1 = SortedPoint1(j1)
                dummy2 = SortedPoint2(j2)
                if( dummy1 <= dummy2 ) then
                  cdf1 = j1 / np1_RK
                  j1 = j1 + 1
                endif
                if( dummy2 <= dummy1 ) then
                  cdf2 = j2 / np2_RK
                  j2 = j2 + 1_IK
                endif
                dt = abs(cdf2-cdf1)
                if (dt>statKS) statKS = dt
                cycle
            end if
            exit
        end do
        npEffective = sqrt( np1_RK * np2_RK / ( np1_RK + np2_RK ) )
        probKS = getProbKS( ( npEffective + 0.12_RK + 0.11_RK / npEffective ) * statKS )
    end subroutine doSortedKS2

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the Kolmogorov–Smirnov (KS) probability.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    pure function getProbKS(lambda) result(probKS)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS
#endif
        implicit none
        real(RK)   , intent(in) :: lambda
        real(RK)   , parameter  :: EPS1 = 0.001_RK, EPS2 = 1.e-8_RK
        integer(IK), parameter  :: NITER = 100
        integer(IK)             :: j
        real(RK)                :: a2,fac,term,termbf
        real(RK)                :: probKS
        a2 = -2._RK*lambda**2
        fac = 2._RK
        probKS = 0._RK
        termbf = 0._RK
        do j = 1, NITER
            term = fac*exp(a2*j**2)
            probKS = probKS+term
            if (abs(term) <= EPS1*termbf .or. abs(term) <= EPS2*probKS) return
            fac = -fac
            termbf = abs(term)
        end do
        probKS = 1._RK
    end function getProbKS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Returns the uniform CDF on support [0,1). This is rather redundant, aint it? but sometimes, needed.
    !>
    !> \param[in] x : The point at which the CDF must be computed.
    !>
    !> \author
    !> Amir Shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    pure function getUniformCDF(x)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniformCDF
#endif
        implicit none
        real(RK), intent(in) :: x
        real(RK)             :: getUniformCDF
        getUniformCDF = x
    end function getUniformCDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the 1-D histogram (Density plot) of the input vector `X`.
    !> The number of bins in the `X` range (`nxbin`) is determined by the user.
    !> The range of `X`, `[xmin, xmax]`, should be also given by the user.
    !> The program returns two arrays of `Xbin` and `Density(x)` as output.
    !>
    !> \param[in]       method          :   The method by which the hist count should be returned:
    !>                                      + `"pdf"`   : Return the probability density function histogram.
    !>                                      + `"count"` : Return the count histogram.
    !> \param[in]       xmin            :   The minimum of the histogram binning.
    !> \param[in]       xmax            :   The maximum of the histogram binning.
    !> \param[in]       nxbin           :   The number of histogram bins.
    !> \param[in]       np              :   The length of input vector `X`.
    !> \param[in]       X               :   The vector of length `nxbin` of values to be binned.
    !> \param[out]      Xbin            :   The vector of length `nxbin` of values representing the bin left corners.
    !> \param[out]      Density         :   The vector of length `nxbin` of values representing the densities in each bin.
    !> \param[out]      errorOccurred   :   The logical output flag indicating whether error has occurred.
    !>
    !> \author
    !> Amir Shahmoradi, Sep 1, 2017, 12:30 AM, ICES, UT Austin
    subroutine getHist1D(method, xmin, xmax, nxbin, np, X, Xbin, Density, errorOccurred)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHist1D
#endif
        implicit none
        character(*), intent(in)  :: method
        integer(IK) , intent(in)  :: np,nxbin
        real(RK)    , intent(in)  :: xmin,xmax
        real(RK)    , intent(in)  :: X(np)
        real(RK)    , intent(out) :: Xbin(nxbin), Density(nxbin)
        logical     , intent(out) :: errorOccurred
        real(RK)                  :: xbinsize
        integer(IK)               :: i, ip, thisXbin, npEffective

        xbinsize = (xmax-xmin) / real(nxbin,kind=RK)
        Xbin = [ (xmin + real(i-1,kind=RK)*xbinsize,i=1,nxbin) ]

        Density = 0._RK
        npEffective = 0_IK
        do ip = 1, np
            if (X(ip)>=xmin .and. X(ip)<xmax) then
                npEffective = npEffective + 1_IK
                thisXbin = getBin(X(ip),xmin,nxbin,xbinsize)
                Density(thisXbin) = Density(thisXbin) + 1._RK
            end if
        end do

        Xbin = Xbin + 0.5_RK * xbinsize

        if(method=="count") then
            errorOccurred = .false.
            return
        elseif(method=="pdf") then
            Density = Density / real(npEffective,kind=RK)
            errorOccurred = .false.
        else
            errorOccurred = .true.
        end if

    end subroutine getHist1D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the 2-D histogram (Density plot) of a set of data points with (X,Y) coordinates.
    !> The number of bins in the `X` and `Y` directions (`[nxbin, nybin]`) are determined by the user.
    !> The range of `X` and `Y` (`xmin`,`xmax`,`ymin`,`ymax`) should also be given by the user.
    !> The program returns three arrays of `Xbin`, `Ybin`, and `Density(y,x)` as the output.
    !>
    !> \param[in]       histType        :   The method by which the normalization of the histogram counts should be done:
    !>                                      + `"count"`     : Return the count histogram.
    !>                                      + `"pdf"`       : Return the probability density function (PDF) histogram.
    !>                                      + `"pdf(y|x)"`  : Return the conditional PDF of `y` given `x` histogram.
    !>                                      + `"pdf(x|y)"`  : Return the conditional PDF of `x` given `y` histogram.
    !> \param[in]       xmin            :   The minimum of the histogram binning along the x-axis.
    !> \param[in]       xmax            :   The maximum of the histogram binning along the x-axis.
    !> \param[in]       ymin            :   The minimum of the histogram binning along the y-axis.
    !> \param[in]       ymax            :   The maximum of the histogram binning along the y-axis.
    !> \param[in]       nxbin           :   The number of histogram bins along the x-axis.
    !> \param[in]       nybin           :   The number of histogram bins along the y-axis.
    !> \param[in]       np              :   The length of input vector `X`.
    !> \param[in]       X               :   The vector of length `nxbin` of values to be binned.
    !> \param[in]       Y               :   The vector of length `nybin` of values to be binned.
    !> \param[out]      Xbin            :   The vector of length `nxbin` of values representing the bin centers.
    !> \param[out]      Ybin            :   The vector of length `nybin` of values representing the bin centers.
    !> \param[out]      Density         :   The array of shape `(nybin,nxbin)` of values representing the densities per bin.
    !> \param[out]      errorOccurred   :   A logical output value indicating whether an error has occurred.
    !>
    !> \author
    !> Amir Shahmoradi, Sep 1, 2017, 12:00 AM, ICES, UT Austin
    subroutine getHist2D(histType,xmin,xmax,ymin,ymax,nxbin,nybin,np,X,Y,Xbin,Ybin,Density,errorOccurred)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getHist2D
#endif
        !use, intrinsic :: iso_fortran_env, only: output_unit
        use String_mod, only: getLowerCase
        implicit none
        character(*), intent(in)  :: histType
        integer(IK) , intent(in)  :: np,nxbin,nybin
        real(RK)    , intent(in)  :: xmin,xmax,ymin,ymax
        real(RK)    , intent(in)  :: X(np),Y(np)
        real(RK)    , intent(out) :: Xbin(nxbin),Ybin(nybin),Density(nybin,nxbin)
        logical     , intent(out) :: errorOccurred
        character(:), allocatable :: method
        real(RK)                  :: xbinsize,ybinsize
        integer(IK)               :: i, ip, thisXbin, thisYbin, npEffective

        errorOccurred = .false.

        xbinsize = (xmax-xmin) / real(nxbin,kind=RK)
        ybinsize = (ymax-ymin) / real(nybin,kind=RK)
        Xbin = [ (xmin+real(i-1,kind=RK)*xbinsize,i=1,nxbin) ]
        Ybin = [ (ymin+real(i-1,kind=RK)*ybinsize,i=1,nybin) ]

        Density = 0._RK
        npEffective = 0_IK
        do ip = 1,np
            if (X(ip)>=xmin .and. X(ip)<xmax .and. Y(ip)>=ymin .and. Y(ip)<ymax) then
                npEffective = npEffective + 1_IK
                thisXbin = getBin(X(ip),xmin,nxbin,xbinsize)
                thisYbin = getBin(Y(ip),ymin,nybin,ybinsize)
                Density(thisYbin,thisXbin) = Density(thisYbin,thisXbin) + 1._RK
            end if
        end do

        Xbin = Xbin + 0.5_RK * xbinsize
        Ybin = Ybin + 0.5_RK * ybinsize
        method = getLowerCase(trim(adjustl(histType)))
        if(method=="pdf") then
            Density = Density / real(npEffective,kind=RK)
        elseif(method=="pdf(y|x)") then
            do i = 1,nxbin
                Density(1:nybin,i) = Density(1:nybin,i) / sum(Density(1:nybin,i))
            end do
        elseif(method=="pdf(x|y)") then
            do i = 1,nybin
                Density(i,1:nxbin) = Density(i,1:nxbin) / sum(Density(i,1:nxbin))
            end do
        elseif(method=="count") then
            return
        else
            errorOccurred = .true.
        end if

    end subroutine getHist2D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Given the range of the variable `x`, `xmin:xmin+binsize*nbin`, and the number of bins, `nbin`, with which
    !> this range is divided, find which bin the input value `x` falls among `[1:nbin]` bins.
    !> The output `ibin` is the number that identifies the bin.
    !>
    !> \param[in]       x               :   The input value whose bin ID is to be found.
    !> \param[in]       lowerBound      :   The lower limit on the value of `x`.
    !> \param[in]       nbin            :   The number of bins to be considered starting from `lowerBound`.
    !> \param[in]       binsize         :   The size of the bins. It must be exactly `(xmax - xmin) / nbin`.
    !>
    !> \return
    !> `ibin` : The ID of the bin to which the input value `x` belongs.
    !>
    !> \warning
    !> If `x <= xmin` or `x xmin + nbin * binsize`, then `ibin = -1` will be returned to indicate error.
    !>
    !> \remark
    !> If `bmin < x <= bmax` then `x` belongs to this bin.
    !>
    !> \todo
    !> The performance and interface of this routine can be significantly improved.
    !> It is more sensible to pass a contiguous array of bin edges as input instead of `lowerBound` and `binsize`.
    !> If the input point is not within any bins, an index of zero should be returned.
    !>
    !> \author
    !> Version 3.0, Sep 1, 2017, 11:12 AM, Amir Shahmoradi, ICES, The University of Texas at Austin.
    pure function getBin(x,lowerBound,nbin,binsize) result(ibin)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getBin
#endif

        implicit none
        integer(IK), intent(in) :: nbin
        real(RK)   , intent(in) :: x,lowerBound,binsize
        real(RK)                :: xmin,xmid
        integer(IK)             :: ibin,minbin,midbin,maxbin

        if (x<lowerBound .or. x>=lowerBound+nbin*binsize) then
            ! LCOV_EXCL_START
            ibin = -1_IK
            return
            ! LCOV_EXCL_STOP
        end if

        minbin = 1
        maxbin = nbin
        xmin = lowerBound
        loopFindBin: do
            midbin = (minbin+maxbin) / 2
            xmid = xmin + midbin*binsize
            if (x<xmid) then
                if (minbin==midbin) then
                    ibin = minbin
                    exit loopFindBin
                end if
                maxbin = midbin
                cycle loopFindBin
            else
                if (minbin==midbin) then
                    ibin = maxbin
                    exit loopFindBin
                end if
                minbin = midbin
                cycle loopFindBin
            end if
        end do loopFindBin

    end function getBin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the quantiles of an input sample of points, given the input quantile probabilities.
    !>
    !> \param[in]       np                          :   The number of points in the input sample.
    !> \param[in]       nq                          :   The number of output quantiles.
    !> \param[in]       SortedQuantileProbability   :   A sorted ascending-order vector of probabilities at which the quantiles will be returned.
    !> \param[in]       Point                       :   The vector of length `np` representing the input sample.
    !> \param[in]       Weight                      :   The vector of length `np` representing the weights of the points in the input sample.
    !> \param[in]       sumWeight                   :   The sum of the vector weights of the points: `sum(Weight)`.
    !>
    !> \return
    !> `Quantile` : The output vector of length `nq`, representing the quantiles corresponding to the input `SortedQuantileProbability` probabilities.
    !>
    !> \author
    !> Amir Shahmoradi, Sep 1, 2017, 12:00 AM, ICES, UT Austin
    pure function getQuantile(np,nq,SortedQuantileProbability,Point,Weight,sumWeight) result(Quantile)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getQuantile
#endif
        use Constants_mod, only: IK, RK, NEGINF_RK
        use Sort_mod, only: indexArray
        use Err_mod, only: Err_type
        implicit none
        integer(IK) , intent(in)            :: np, nq
        real(RK)    , intent(in)            :: SortedQuantileProbability(nq), Point(np)
        integer(IK) , intent(in), optional  :: Weight(np), sumWeight
        real(RK)                            :: Quantile(nq)
        integer(IK)                         :: ip, iq, iw, weightCounter, Indx(np), SortedQuantileDensity(nq)
        type(Err_type)                      :: Err
        call indexArray(np,Point,Indx,Err)
        if (Err%occurred) then
            ! LCOV_EXCL_START
            Quantile = NEGINF_RK
            return
            ! LCOV_EXCL_STOP
        end if
        iq = 1_IK
        if (present(sumWeight)) then
            weightCounter = 0_IK
            SortedQuantileDensity = nint( SortedQuantileProbability * sumWeight )
            loopWeighted: do ip = 1, np
                do iw = 1, Weight(Indx(ip))
                    weightCounter = weightCounter + 1_IK
                    if (weightCounter>=SortedQuantileDensity(iq)) then
                        Quantile(iq) = Point(Indx(ip))
                        iq = iq + 1_IK
                        if (iq>nq) exit loopWeighted
                    end if
                end do
            end do loopWeighted
        else
            SortedQuantileDensity = nint( SortedQuantileProbability * np )
            loopNonWeighted: do ip = 1, np
                if (ip>=SortedQuantileDensity(iq)) then
                    Quantile(iq) = Point(Indx(ip))
                    iq = iq + 1_IK
                    if (iq>nq) exit loopNonWeighted
                end if
            end do loopNonWeighted
        end if
        if (iq<=nq) Quantile(iq:nq) = Point(Indx(np))
    end function getQuantile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine getClusteredPoint( self & ! LCOV_EXCL_LINE
                                , nd, ndmin, ndmax & ! LCOV_EXCL_LINE
                                , nc, ncmin, ncmax & ! LCOV_EXCL_LINE
                                , Size, sizeMin, sizeMax & ! LCOV_EXCL_LINE
                                , Center, centerMin, centerMax & ! LCOV_EXCL_LINE
                                , Std, stdmin, stdmax & ! LCOV_EXCL_LINE
                                , Eta, etaMin, etaMax & ! LCOV_EXCL_LINE
                                , dist & ! LCOV_EXCL_LINE
                                )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getClusteredPoint
#endif
        use Matrix_mod, only: getInvMatFromCholFac
        use Matrix_mod, only: getCholeskyFactor
        use Constants_mod, only: IK, RK
        use String_mod, only: num2str !, getLowerCase
        use Math_mod, only: getCumSum
        implicit none
        class(ClusteredPoint_type), intent(inout)       :: self
        integer(IK) , intent(in), optional              :: nd, ndmin, ndmax
        integer(IK) , intent(in), optional              :: nc, ncmin, ncmax
        integer(IK) , intent(in), optional              :: sizeMin, sizeMax
        real(RK)    , intent(in), optional              :: stdMin, stdMax
        real(RK)    , intent(in), optional              :: etaMin, etaMax
        real(RK)    , intent(in), optional              :: centerMin, centerMax
        real(RK)    , intent(in), optional              :: Center(:,:), Std(:,:), Eta(:)
        integer(IK) , intent(in), optional              :: Size(:)
        character(*), intent(in), optional              :: dist
        integer(IK) , allocatable                       :: CumSumSize(:)
        real(RK)    , allocatable                       :: NormedPoint(:), InvCovMat(:,:,:)
        real(RK)    , allocatable                       :: CumSumVolNormed(:)
        real(RK)    , allocatable                       :: VolNormed(:)
        real(RK)                                        :: dummy, minLogVol
        logical                                         :: isUniformSuperposed
        logical                                         :: isUniform, isMember
        logical                                         :: isNormal
        integer(IK)                                     :: i, j, ic, ip, membershipCount, minSize

        self%Err%occurred = .false.

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! ndmin, ndmax, nd
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (present(ndmin)) then
            self%ndmin = ndmin
        else
            self%ndmin = 2_IK
        endif

        if (present(ndmax)) then
            self%ndmax = ndmax
        else
            self%ndmax = self%ndmin * 10_IK
        endif

        if (present(nd)) then
            self%nd = nd
        else
            self%nd = getRandInt(self%ndmin, self%ndmax)
        endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! ncmin, ncmax, nc
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (present(ncmin)) then
            self%ncmin = ncmin
        else
            self%ncmin = 1_IK
        endif

        if (present(ncmax)) then
            self%ncmax = ncmax
        else
            self%ncmax = self%ncmin * 10_IK
        endif

        if (present(nc)) then
            self%nc = nc
        else
            self%nc = getRandInt(self%ncmin, self%ncmax)
        endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! sizeMin, sizeMax, Size
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (present(sizeMin)) then
            self%sizeMin = sizeMin
        else
            self%sizeMin = self%nd + 1
        endif

        if (present(sizeMax)) then
            self%sizeMax = sizeMax
        else
            self%sizeMax = self%sizeMin * 1000_IK
        endif

        if (present(Size)) then
            self%Size = Size
        else
            allocate(self%Size(self%nc), CumSumSize(0:self%nc))
            CumSumSize(0) = 0._RK
            do ic = 1, self%nc
                self%Size(ic) = getRandInt(self%sizeMin, self%sizeMax)
                CumSumSize(ic) = CumSumSize(ic-1) + self%Size(ic)
            end do
        endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! centerMin, centerMax, Size
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (present(centerMin)) then
            self%centerMin = centerMin
        else
            self%centerMin = 0._RK
        endif

        if (present(centerMax)) then
            self%centerMax = centerMax
        else
            self%centerMax = 1._RK
        endif

        if (present(Center)) then
            self%Center = Center
        else
            allocate(self%Center(self%nd,self%nc))
            call random_number(self%Center)
            do concurrent(ic = 1:self%nc)
                self%Center(:,ic) = self%centerMin + self%Center(:,ic) * (self%centerMax - self%centerMin)
            end do
        endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! stdMin, stdMax, Size
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (present(stdMin)) then
            self%stdMin = stdMin
        else
            self%stdMin = 0._RK
        endif

        if (present(stdMax)) then
            self%stdMax = stdMax
        else
            self%stdMax = 1._RK
        endif

        if (present(Std)) then
            self%Std = Std
        else
            allocate(self%Std(self%nd,self%nc))
            call random_number(self%Std)
            do concurrent(ic = 1:self%nc)
                self%Std(:,ic) = self%stdMin + self%Std(:,ic) * (self%stdMax - self%stdMin)
            end do
        endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! etaMin, etaMax, eta
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (present(etaMin)) then
            self%etaMin = etaMin
        else
            self%etaMin = 0.1_RK
        endif

        if (present(etaMax)) then
            self%etaMax = etaMax
        else
            self%etaMax = self%etaMin * (1._RK/self%etaMin)**2
        endif

        if (present(Eta)) then
            self%Eta = Eta
        else
            allocate(self%Eta(self%nc))
            call random_number(self%Eta)
            do concurrent(ic = 1:self%nc)
                self%Eta(ic) = self%etaMin + self%Eta(ic) * (self%etaMax - self%etaMin)
            end do
        endif

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Generate the representative matrices of the clusters
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        allocate(self%LogVolume(self%nc))
        allocate(self%ChoDia(self%nd,self%nc))
        allocate(self%ChoLowCovUpp(self%nd,self%nd,self%nc))
        do ic = 1, self%nc

            ! Generate random correlation matrix.

            self%ChoLowCovUpp(:,:,ic) = getRandCorMat(self%nd, self%Eta(ic))
            do j = 1, self%nd
                do i = 1, self%nd
                    self%ChoLowCovUpp(i,j,ic) = self%ChoLowCovUpp(i,j,ic) * self%Std(i,ic) * self%Std(j,ic)
                end do
            end do

            ! Compute the Cholesky factorization.

            call getCholeskyFactor  ( nd = self%nd & ! LCOV_EXCL_LINE
                                    , PosDefMat = self%ChoLowCovUpp(1:self%nd,1:self%nd,ic) & ! LCOV_EXCL_LINE
                                    , Diagonal = self%ChoDia(1:self%nd,ic) & ! LCOV_EXCL_LINE
                                    )
            if (self%ChoDia(1,ic) < 0._RK) then
                self%Err%occurred = .true.
                self%Err%msg = "Singular Covariance matrix detected."
                write(*,"(A)") "ChoLowCovUpp:"
                write(*,"("//num2str(self%nd)//"(F15.8,:,' '))") self%ChoLowCovUpp(1:self%nd,1:self%nd,ic)
                error stop
            end if

            self%LogVolume(ic) = sum(log(self%ChoDia(1:self%nd,ic)))

        end do

        VolNormed = exp( self%LogVolume - maxval(self%LogVolume) )
        CumSumVolNormed = getCumSum(self%nc, VolNormed)
        CumSumVolNormed = CumSumVolNormed / CumSumVolNormed(self%nc)
        

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Set the distribution
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (present(dist)) then
            self%dist = dist
        else
            self%dist = "uniform-mixture"
        end if

        isUniform = self%dist == "uniform"
        isNormal = self%dist == "normal-mixture"
        isUniformSuperposed = self%dist == "uniform-mixture"

        if (.not. (isUniformSuperposed .or. isUniform)) then
            self%Err%occurred = .true.
            self%Err%msg = "No point distribution other than uniform is currently supported."
            return
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Generate the cluster members
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        self%np = sum(self%Size)
        allocate(self%Membership(self%np))
        allocate(self%Point(self%nd,self%np))

        if (isUniformSuperposed) then

            do ic = 1, self%nc
                do ip = CumSumSize(ic-1) + 1, CumSumSize(ic)
                    self%Membership(ip) = ic
                    self%Point(1:self%nd,ip) = getRandMVU   ( nd = self%nd & ! LCOV_EXCL_LINE
                                                            , MeanVec = self%Center(1:self%nd,ic) & ! LCOV_EXCL_LINE
                                                            , CholeskyLower = self%ChoLowCovUpp(1:self%nd,1:self%nd,ic) & ! LCOV_EXCL_LINE
                                                            , Diagonal = self%ChoDia(1:self%nd,ic) & ! LCOV_EXCL_LINE
                                                            )
                end do
            end do

        elseif (isUniform) then

            allocate(NormedPoint(self%nd), InvCovMat(self%nd,self%nd,self%nc))

            do ic = 1, self%nc
                InvCovMat(1:self%nd,1:self%nd,ic) = getInvMatFromCholFac( nd = self%nd & ! LCOV_EXCL_LINE
                                                                        , CholeskyLower = self%ChoLowCovUpp(1:self%nd,1:self%nd,ic) & ! LCOV_EXCL_LINE
                                                                        , CholeskyDiago = self%ChoDia(1:self%nd,ic) & ! LCOV_EXCL_LINE
                                                                        )
            end do

            self%Size = 0_IK ! all cluster sizes will be determined after all point generations.

            do ip = 1, self%np

                ! Check for the distribution type.

                loopAddClusterMember: do

                    call random_number(dummy)
                    ic = minloc(CumSumVolNormed, dim = 1, mask = CumSumVolNormed > dummy)
                    !ic = getRandInt(lowerBound = 1_IK, upperBound = self%nc)

                    self%Point(1:self%nd,ip) = getRandMVU   ( nd = self%nd & ! LCOV_EXCL_LINE
                                                            , MeanVec = self%Center(1:self%nd,ic) & ! LCOV_EXCL_LINE
                                                            , CholeskyLower = self%ChoLowCovUpp(1:self%nd,1:self%nd,ic) & ! LCOV_EXCL_LINE
                                                            , Diagonal = self%ChoDia(1:self%nd,ic) & ! LCOV_EXCL_LINE
                                                            )

                    membershipCount = 0_IK
                    do i = 1, self%nc
                        NormedPoint(1:self%nd) = self%Point(1:self%nd,ip) - self%Center(1:self%nd,i)
                        isMember = isInsideEllipsoid(nd = self%nd, NormedPoint = NormedPoint, InvRepMat = InvCovMat(1:self%nd,1:self%nd,i))
                        if (isMember) membershipCount = membershipCount + 1_IK
                    end do

                    call random_number(dummy)
                    if (dummy < 1._RK / membershipCount) then
                        self%Size(ic) = self%Size(ic) + 1_IK
                        self%Membership(ip) = ic
                        exit loopAddClusterMember
                    end if

                    cycle loopAddClusterMember

                end do loopAddClusterMember

            end do

        end if

    end subroutine getClusteredPoint

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine writeClusteredPoint(self, fileUnit)
        use Constants_mod, only: IK, RK
        implicit none
        class(ClusteredPoint_type)  , intent(in)    :: self
        integer(IK)                 , intent(in)    :: fileUnit
        character(*), parameter                     :: fileFormat = "(*(g0.8,:,','))"
        integer(IK)                                 :: i, j, ic

        write(fileUnit,"(A)") "dist"
        write(fileUnit,fileFormat) self%dist

        write(fileUnit,"(A)") "nd, np, nc"
        write(fileUnit,fileFormat) self%nd, self%np, self%nc

        write(fileUnit,"(A)") "Size"
        write(fileUnit,fileFormat) self%Size

        write(fileUnit,"(A)") "Center"
        write(fileUnit,fileFormat) self%Center

        write(fileUnit,"(A)") "LogVolume"
        write(fileUnit,fileFormat) self%LogVolume

        write(fileUnit,"(A)") "CholeskyLower"
        write(fileUnit,fileFormat) ((self%ChoDia(j,ic), (self%ChoLowCovUpp(i,j,ic), i=j+1,self%nd), j=1,self%nd), ic=1,self%nc)

        write(fileUnit,"(A)") "Point"
        write(fileUnit,fileFormat) self%Point

        write(fileUnit,"(A)") "Membership"
        write(fileUnit,fileFormat) self%Membership

    end subroutine writeClusteredPoint

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Statistics_mod ! LCOV_EXCL_LINE