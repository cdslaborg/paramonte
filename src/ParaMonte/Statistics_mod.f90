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
!!!!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Statistics_mod

    use Constants_mod, only: RK, IK

    implicit none

!logical, save :: paradramPrintEnabled = .false.
!logical, save :: paradisePrintEnabled = .false.

    character(len=*), parameter :: MODULE_NAME = "@Statistics_mod"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface getMean
        module procedure :: getMean_2D
    end interface getMean

    interface getNormData
        module procedure :: getNormData_2D
    end interface getNormData

    interface getVariance
        module procedure :: getVariance_1D, getVariance_2D
    end interface getVariance

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface getLogProbLogNorm
        module procedure :: getLogProbLogNormS, GetLogProbLogNormMP
    end interface getLogProbLogNorm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface GetLogProbNormSP
        module procedure :: getLogProbNormSP_RK, getLogProbNormSP_CK
    end interface GetLogProbNormSP

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface GetLogProbNormMP
        module procedure :: GetLogProbNormMP_RK, GetLogProbNormMP_CK
    end interface GetLogProbNormMP

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface getProbMVNSP
        module procedure :: getProbMVNSP_RK, getProbMVNSP_CK
    end interface getProbMVNSP

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface getProbMVNV
        module procedure :: getProbMVNMP_RK, getProbMVNMP_CK
    end interface getProbMVNV

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface getLogProbMVNSP
        module procedure :: getLogProbMVNSP_RK, getLogProbMVNSP_CK
    end interface getLogProbMVNSP

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface getLogProbMVNMP
        module procedure :: getLogProbMVNMP_RK, getLogProbMVNMP_CK
    end interface getLogProbMVNMP

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface getLogProbNorm
        module procedure :: getLogProbNormSP_RK, GetLogProbNormMP_RK
        module procedure :: getLogProbNormSP_CK, GetLogProbNormMP_CK
    end interface getLogProbNorm

    interface getProbMVN
        module procedure :: getProbMVNSP_RK, getProbMVNMP_RK
        module procedure :: getProbMVNSP_CK, getProbMVNMP_CK
    end interface getProbMVN

    interface getLogProbMVN
        module procedure :: getLogProbMVNSP_RK, getLogProbMVNMP_RK
        module procedure :: getLogProbMVNSP_CK, getLogProbMVNMP_CK
    end interface getLogProbMVN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface getLogProbGausMixSDSP
        module procedure :: getLogProbGausMixSDSP_RK, getLogProbGausMixSDSP_CK
    end interface getLogProbGausMixSDSP

    interface getLogProbGausMixSDMP
        module procedure :: getLogProbGausMixSDMP_RK, getLogProbGausMixSDMP_CK
    end interface getLogProbGausMixSDMP

    interface getLogProbGausMixMDSP
        module procedure :: getLogProbGausMixMDSP_RK, getLogProbGausMixMDSP_CK
    end interface getLogProbGausMixMDSP

    interface getLogProbGausMixMDMP
        module procedure :: getLogProbGausMixMDMP_RK, getLogProbGausMixMDMP_CK
    end interface getLogProbGausMixMDMP

    interface getLogProbGausMix_RK
        module procedure :: getLogProbGausMixSDSP_RK, getLogProbGausMixSDMP_RK, getLogProbGausMixMDSP_RK, getLogProbGausMixMDMP_RK
    end interface getLogProbGausMix_RK

    interface getLogProbGausMix
        module procedure :: getLogProbGausMixSDSP_RK, getLogProbGausMixSDMP_RK, getLogProbGausMixMDSP_RK, getLogProbGausMixMDMP_RK
        module procedure :: getLogProbGausMixSDSP_CK, getLogProbGausMixSDMP_CK, getLogProbGausMixMDSP_CK, getLogProbGausMixMDMP_CK
    end interface getLogProbGausMix

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface getMahalSqSP
        module procedure :: getMahalSqSP_RK, getMahalSqSP_CK
    end interface getMahalSqSP

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface getMahalSqMP
        module procedure :: getMahalSqMP_RK, getMahalSqMP_CK
    end interface getMahalSqMP

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    interface getMahalSq
        module procedure :: getMahalSqSP_RK, getMahalSqMP_RK
        module procedure :: getMahalSqSP_CK, getMahalSqMP_CK
    end interface getMahalSq

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    ! returns square of Mahalanobis distance for a single point. output is a scalar variable.
    ! NOTE: if computation fails, the first element of output will be returned negative.
    pure function getMahalSqSP_RK(nd,MeanVec,InvCovMat,Point)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getMahalSqSP_RK
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: MeanVec(nd)
        real(RK)   , intent(in) :: InvCovMat(nd,nd)        ! Inverse of the covariance matrix
        real(RK)   , intent(in) :: Point(nd)               ! input data points
        real(RK)                :: getMahalSqSP_RK
        real(RK)                :: NormedPoint(nd)
        NormedPoint = Point - MeanVec
        getMahalSqSP_RK = dot_product( NormedPoint , matmul(InvCovMat,NormedPoint) )
    end function getMahalSqSP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns square of the Mahalanobis distance
    ! NOTE: if computation fails, the first element of output (if array) will be returned negative.
    pure function getMahalSqMP_RK(nd,np,MeanVec,InvCovMat,Point)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getMahalSqMP_RK
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK), intent(in) :: nd,np
        real(RK), intent(in)    :: MeanVec(nd)
        real(RK), intent(in)    :: InvCovMat(nd,nd)         ! Inverse of the covariance matrix
        real(RK), intent(in)    :: Point(nd,np)             ! input data points
        real(RK)                :: getMahalSqMP_RK(np)      ! function return
        real(RK)                :: NormedPoint(nd)
        integer(IK)             :: ip
        do ip = 1,np
            NormedPoint = Point(1:nd,ip) - MeanVec
            getMahalSqMP_RK(ip) = dot_product( NormedPoint , matmul(InvCovMat,NormedPoint) )
            if (getMahalSqMP_RK(ip)<0._RK) then
                getMahalSqMP_RK(1) = -1._RK
                return
            end if
        end do
    end function getMahalSqMP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns square of Mahalanobis distance for a single Point_CK. output is a scalar variable.
    ! NOTE: if computation fails, the first element of output will be returned negative.
    pure function getMahalSqSP_CK(nd,MeanVec_CK,InvCovMat_CK,Point_CK)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getMahalSqSP_CK
#endif
        use Constants_mod, only: IK, RK, CK
        implicit none
        integer(IK), intent(in)  :: nd
        complex(CK), intent(in)  :: MeanVec_CK(nd)
        complex(CK), intent(in)  :: InvCovMat_CK(nd,nd)        ! Inverse of the covariance matrix
        complex(CK), intent(in)  :: Point_CK(nd)               ! input data points
        complex(CK)              :: getMahalSqSP_CK            ! function return
        getMahalSqSP_CK = sum( (Point_CK-MeanVec_CK) * matmul(InvCovMat_CK,Point_CK-MeanVec_CK) )
    end function getMahalSqSP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns square of the Mahalanobis distance
    ! NOTE: if computation fails, the first element of output will be returned negative.
    pure function getMahalSqMP_CK(nd,np,MeanVec_CK,InvCovMat_CK,Point_CK)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getMahalSqMP_CK
#endif
        use Constants_mod, only: IK, RK, CK
        implicit none
        integer(IK), intent(in)  :: nd,np
        complex(CK), intent(in)  :: MeanVec_CK(nd)
        complex(CK), intent(in)  :: InvCovMat_CK(nd,nd)       ! Inverse of the covariance matrix
        complex(CK), intent(in)  :: Point_CK(nd,np)           ! input data points
        complex(CK)              :: getMahalSqMP_CK(np)       ! function return
        integer(IK)              :: ip
        do ip = 1,np
            getMahalSqMP_CK(ip) = sum( (Point_CK(1:nd,ip)-MeanVec_CK) * &
            matmul(InvCovMat_CK,Point_CK(1:nd,ip)-MeanVec_CK) )
            if (real(getMahalSqMP_CK(ip))<0._RK) then
                getMahalSqMP_CK(1) = (-1._RK, -1._RK)
                return
            end if
        end do
    end function getMahalSqMP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogProbNormSP_RK(mean,inverseVariance,logSqrtInverseVariance,point)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbNormSP_RK
#endif
        use Constants_mod, only: RK, LOGINVSQRT2PI
        implicit none
        real(RK), intent(in) :: mean,inverseVariance,logSqrtInverseVariance,point
        real(RK)             :: getLogProbNormSP_RK
        getLogProbNormSP_RK = LOGINVSQRT2PI + logSqrtInverseVariance - 0.5_RK * inverseVariance * (point-mean)**2
    end function getLogProbNormSP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function GetLogProbNormMP_RK(np,mean,inverseVariance,logSqrtInverseVariance,Point)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: GetLogProbNormMP_RK
#endif
        use Constants_mod, only: LOGINVSQRT2PI
        implicit none
        integer(IK), intent(in) :: np
        real(RK)   , intent(in) :: mean,inverseVariance,logSqrtInverseVariance,Point(np)
        real(RK)                :: GetLogProbNormMP_RK(np)
        GetLogProbNormMP_RK = LOGINVSQRT2PI + logSqrtInverseVariance - 0.5_RK * inverseVariance * (Point-mean)**2
    end function GetLogProbNormMP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! NOTE: if MahalSq computation fails, output probability will be returned as NullVal%RK from module Constant_mod.
    pure function getProbMVNSP_RK(nd,MeanVec,InvCovMat,sqrtDetInvCovMat,Point)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbMVNSP_RK
#endif
        use Constants_mod, only: INVSQRT2PI, NullVal
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: MeanVec(nd)
        real(RK)   , intent(in) :: InvCovMat(nd,nd)
        real(RK)   , intent(in) :: sqrtDetInvCovMat
        real(RK)   , intent(in) :: Point(nd)
        real(RK)                :: getProbMVNSP_RK, dummy
        dummy = getMahalSqSP_RK(nd,MeanVec,InvCovMat,Point)
        if (dummy<0._RK) then
            getProbMVNSP_RK = NullVal%RK
        else
            getProbMVNSP_RK = INVSQRT2PI**nd * sqrtDetInvCovMat * exp( -0.5_RK * dummy )
        end if
    end function getProbMVNSP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! NOTE: if MahalSq computation fails, output probability will be returned as NullVal%RK from module Constant_mod.
    pure function getProbMVNMP_RK(nd,np,MeanVec,InvCovMat,sqrtDetInvCovMat,Point)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbMVNMP_RK
#endif
        use Constants_mod, only: INVSQRT2PI, NullVal
        implicit none
        integer(IK), intent(in) :: nd,np
        real(RK)   , intent(in) :: MeanVec(nd)
        real(RK)   , intent(in) :: InvCovMat(nd,nd)
        real(RK)   , intent(in) :: sqrtDetInvCovMat
        real(RK)   , intent(in) :: Point(nd,np)
        real(RK)                :: getProbMVNMP_RK(np), Dummy(np)
        Dummy = getMahalSqMP_RK(nd,np,MeanVec,InvCovMat,Point)
        if (Dummy(1)<0._RK) then
            getProbMVNMP_RK = NullVal%RK
        else
            getProbMVNMP_RK = INVSQRT2PI**nd * sqrtDetInvCovMat * exp( -0.5_RK * Dummy )
        end if
    end function getProbMVNMP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! NOTE: if MahalSq computation fails, output probability will be returned as NullVal%RK from module Constant_mod.
    pure function getLogProbMVNSP_RK(nd,MeanVec,InvCovMat,logSqrtDetInvCovMat,Point)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMVNSP_RK
#endif
        use Constants_mod, only: LOGINVSQRT2PI, NullVal
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: MeanVec(nd)
        real(RK)   , intent(in) :: InvCovMat(nd,nd)
        real(RK)   , intent(in) :: logSqrtDetInvCovMat
        real(RK)   , intent(in) :: Point(nd)
        real(RK)                :: getLogProbMVNSP_RK, dummy
        dummy = getMahalSqSP_RK(nd,MeanVec,InvCovMat,Point)
        if (dummy<0._RK) then
            getLogProbMVNSP_RK = NullVal%RK
        else
            getLogProbMVNSP_RK = nd*LOGINVSQRT2PI + logSqrtDetInvCovMat - 0.5_RK * dummy
        end if
    end function getLogProbMVNSP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! NOTE: if MahalSq computation fails, output probability will be returned as NullVal%RK from module Constant_mod.
    pure function getLogProbMVNMP_RK(nd,np,MeanVec,InvCovMat,logSqrtDetInvCovMat,Point)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMVNMP_RK
#endif
        use Constants_mod, only: LOGINVSQRT2PI, NullVal
        implicit none
        integer(IK), intent(in) :: nd,np
        real(RK)   , intent(in) :: MeanVec(nd)
        real(RK)   , intent(in) :: InvCovMat(nd,nd)
        real(RK)   , intent(in) :: logSqrtDetInvCovMat
        real(RK)   , intent(in) :: Point(nd,np)
        real(RK)                :: getLogProbMVNMP_RK(np), Dummy(np)
        Dummy = getMahalSqMP_RK(nd,np,MeanVec,InvCovMat,Point)
        if (Dummy(1)<0._RK) then
            getLogProbMVNMP_RK = NullVal%RK
        else
            getLogProbMVNMP_RK = nd*LOGINVSQRT2PI + logSqrtDetInvCovMat - 0.5_RK * Dummy
        end if
    end function getLogProbMVNMP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getLogProbNormSP_CK(mean_CK,inverseVariance_CK,logSqrtInverseVariance_CK,point_CK)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbNormSP_CK
#endif
        use Constants_mod, only: RK, CK, LOGINVSQRT2PI
        implicit none
        complex(CK), intent(in) :: mean_CK,inverseVariance_CK,logSqrtInverseVariance_CK,point_CK
        complex(CK)             :: getLogProbNormSP_CK
        getLogProbNormSP_CK = LOGINVSQRT2PI + logSqrtInverseVariance_CK - 0.5_RK * inverseVariance_CK * (point_CK-mean_CK)**2
    end function getLogProbNormSP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function GetLogProbNormMP_CK(np,mean_CK,inverseVariance_CK,logSqrtInverseVariance_CK,Point_CK)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: GetLogProbNormMP_CK
#endif
        use Constants_mod, only: IK, RK, CK, LOGINVSQRT2PI
        implicit none
        integer(IK), intent(in) :: np
        complex(CK)   , intent(in) :: mean_CK,inverseVariance_CK,logSqrtInverseVariance_CK,Point_CK(np)
        complex(CK)                :: GetLogProbNormMP_CK(np)
        GetLogProbNormMP_CK = LOGINVSQRT2PI + logSqrtInverseVariance_CK - 0.5_RK * inverseVariance_CK * (Point_CK-mean_CK)**2
    end function GetLogProbNormMP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getProbMVNSP_CK(nd,MeanVec_CK,InvCovMat_CK,sqrtDetInvCovMat_CK,Point_CK)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbMVNSP_CK
#endif
        use Constants_mod, only: IK, RK, CK, INVSQRT2PI, NullVal
        implicit none
        integer(IK), intent(in) :: nd
        complex(CK)   , intent(in) :: MeanVec_CK(nd)
        complex(CK)   , intent(in) :: InvCovMat_CK(nd,nd)
        complex(CK)   , intent(in) :: sqrtDetInvCovMat_CK
        complex(CK)   , intent(in) :: Point_CK(nd)
        complex(CK)                :: getProbMVNSP_CK, dummy
        dummy = getMahalSqSP_CK(nd,MeanVec_CK,InvCovMat_CK,Point_CK)
        if (real(dummy)<0._RK) then
            getProbMVNSP_CK = NullVal%RK
        else
            getProbMVNSP_CK = INVSQRT2PI**nd * sqrtDetInvCovMat_CK * exp( -0.5_RK * dummy )
        end if
    end function getProbMVNSP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getProbMVNMP_CK(nd,np,MeanVec_CK,InvCovMat_CK,sqrtDetInvCovMat_CK,Point_CK)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbMVNMP_CK
#endif
        use Constants_mod, only: IK, RK, CK, INVSQRT2PI, NullVal
        implicit none
        integer(IK), intent(in) :: nd,np
        complex(CK), intent(in) :: MeanVec_CK(nd)
        complex(CK), intent(in) :: InvCovMat_CK(nd,nd)
        complex(CK), intent(in) :: sqrtDetInvCovMat_CK
        complex(CK), intent(in) :: Point_CK(nd,np)
        complex(CK)             :: getProbMVNMP_CK(np), Dummy(np)
        Dummy = getMahalSqMP_CK(nd,np,MeanVec_CK,InvCovMat_CK,Point_CK)
        if (real(Dummy(1))<0._RK) then
            getProbMVNMP_CK = NullVal%RK
        else
            getProbMVNMP_CK = INVSQRT2PI**nd * sqrtDetInvCovMat_CK * exp( -0.5_RK * Dummy )
        end if
    end function getProbMVNMP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getLogProbMVNSP_CK(nd,MeanVec_CK,InvCovMat_CK,logSqrtDetInvCovMat_CK,Point_CK)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMVNSP_CK
#endif
        use Constants_mod, only: IK, RK, CK, LOGINVSQRT2PI, NullVal
        implicit none
        integer(IK), intent(in) :: nd
        complex(CK), intent(in) :: MeanVec_CK(nd)
        complex(CK), intent(in) :: InvCovMat_CK(nd,nd)
        complex(CK), intent(in) :: logSqrtDetInvCovMat_CK
        complex(CK), intent(in) :: Point_CK(nd)
        complex(CK)             :: getLogProbMVNSP_CK, dummy
        dummy = getMahalSqSP_CK(nd,MeanVec_CK,InvCovMat_CK,Point_CK)
        if (real(dummy)<0._RK) then
            getLogProbMVNSP_CK = NullVal%RK
        else
            getLogProbMVNSP_CK  = nd*LOGINVSQRT2PI + logSqrtDetInvCovMat_CK - 0.5_RK * dummy
        end if
    end function getLogProbMVNSP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getLogProbMVNMP_CK(nd,np,MeanVec_CK,InvCovMat_CK,logSqrtDetInvCovMat_CK,Point_CK)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMVNMP_CK
#endif
        use Constants_mod, only: IK, RK, CK, LOGINVSQRT2PI, NullVal
        implicit none
        integer(IK), intent(in) :: nd,np
        complex(CK), intent(in) :: MeanVec_CK(nd)
        complex(CK), intent(in) :: InvCovMat_CK(nd,nd)
        complex(CK), intent(in) :: logSqrtDetInvCovMat_CK
        complex(CK), intent(in) :: Point_CK(nd,np)
        complex(CK)             :: getLogProbMVNMP_CK(np), Dummy(np)
        Dummy = getMahalSqMP_CK(nd,np,MeanVec_CK,InvCovMat_CK,Point_CK)
        if (real(Dummy(1))<0._RK) then
            getLogProbMVNMP_CK = NullVal%RK
        else
            getLogProbMVNMP_CK  = nd*LOGINVSQRT2PI + logSqrtDetInvCovMat_CK - 0.5_RK * Dummy
        end if
    end function getLogProbMVNMP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! SDSP stands for Single-Dimensional Gaussian mixture with Single Point input
    function getLogProbGausMixSDSP_RK(nmode,nd,np,LogAmplitude,MeanVec,InvCovMat,LogSqrtDetInvCovMat,point)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbGausMixSDSP_RK
#endif
        use Constants_mod, only: IK, RK, LOGTINY_RK
        implicit none
        integer(IK), intent(in) :: nmode,nd,np
        real(RK)   , intent(in) :: LogAmplitude(nmode),MeanVec(nmode)
        real(RK)   , intent(in) :: InvCovMat(nmode),LogSqrtDetInvCovMat(nmode)
        real(RK)   , intent(in) :: point
        real(RK)                :: getLogProbGausMixSDSP_RK
        real(RK)                :: normFac,LogProb(nmode)
        integer(IK)             :: imode
        do imode = 1, nmode
            LogProb(imode)  = LogAmplitude(imode) &
                            + getLogProbNormSP_RK(MeanVec(imode),InvCovMat(imode),LogSqrtDetInvCovMat(imode),point)
        end do
        normFac = maxval(LogProb)
        LogProb = LogProb - normFac
        where(LogProb<LOGTINY_RK)
            LogProb = 0._RK
        elsewhere
            LogProb = exp(LogProb)
        end where
        getLogProbGausMixSDSP_RK = normFac + log(sum(LogProb))
    end function getLogProbGausMixSDSP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function getLogProbGausMixSDMP_RK(nmode,nd,np,LogAmplitude,MeanVec,InvCovMat,LogSqrtDetInvCovMat,Point)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbGausMixSDMP_RK
#endif
        use Constants_mod, only: IK, RK, LOGTINY_RK
        implicit none
        integer(IK), intent(in) :: nmode,nd,np
        real(RK)   , intent(in) :: LogAmplitude(nmode),MeanVec(nmode)
        real(RK)   , intent(in) :: InvCovMat(nmode),LogSqrtDetInvCovMat(nmode)
        real(RK)   , intent(in) :: Point(np)
        real(RK)                :: getLogProbGausMixSDMP_RK(np)
        real(RK)                :: NormFac(np),LogProb(nmode,np)
        integer(IK)             :: imode, ip
        do imode = 1, nmode
            LogProb(imode,1:np) = LogAmplitude(imode) &
                                + getLogProbNormMP_RK(np,MeanVec(imode),InvCovMat(imode),LogSqrtDetInvCovMat(imode),Point)
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
            getLogProbGausMixSDMP_RK(ip) = NormFac(ip) + log(sum(LogProb(1:nmode,ip)))
        end do
    end function getLogProbGausMixSDMP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getLogProbGausMixMDSP_RK(nmode,nd,np,LogAmplitude,MeanVec,InvCovMat,LogSqrtDetInvCovMat,Point)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbGausMixMDSP_RK
#endif
        use Constants_mod, only: IK, RK, LOGTINY_RK
        implicit none
        integer(IK), intent(in) :: nmode,nd,np
        real(RK)   , intent(in) :: LogAmplitude(nmode), MeanVec(nd,nmode)
        real(RK)   , intent(in) :: InvCovMat(nd,nd,nmode), LogSqrtDetInvCovMat(nmode)
        real(RK)   , intent(in) :: Point(nd)
        real(RK)                :: getLogProbGausMixMDSP_RK
        real(RK)                :: normFac,LogProb(nmode)
        integer(IK)             :: imode
        do imode = 1, nmode
            LogProb(imode) = LogAmplitude(imode) + &
            getLogProbMVNSP_RK(nd,MeanVec(1:nd,imode),InvCovMat(1:nd,1:nd,imode),LogSqrtDetInvCovMat(imode),Point)
        end do
        normFac = maxval(LogProb)
        LogProb = LogProb - normFac
        where(LogProb<LOGTINY_RK)
            LogProb = 0._RK
        elsewhere
            LogProb = exp(LogProb)
        end where
        getLogProbGausMixMDSP_RK = normFac + log(sum(LogProb))
    end function getLogProbGausMixMDSP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function getLogProbGausMixMDMP_RK(nmode,nd,np,LogAmplitude,MeanVec,InvCovMat,LogSqrtDetInvCovMat,Point)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbGausMixMDMP_RK
#endif
        use Constants_mod, only: IK, RK, LOGTINY_RK
        implicit none
        integer(IK), intent(in) :: nmode,nd,np
        real(RK)   , intent(in) :: LogAmplitude(nmode),MeanVec(nd,nmode)
        real(RK)   , intent(in) :: InvCovMat(nd,nd,nmode), LogSqrtDetInvCovMat(nmode)
        real(RK)   , intent(in) :: Point(nd,np)
        real(RK)                :: getLogProbGausMixMDMP_RK(np)
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
        getLogProbGausMixMDMP_RK(ip) = NormFac(ip) + log(sum(LogProb(1:nmode,ip)))
        end do
    end function getLogProbGausMixMDMP_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! SDSP stands for 1-dimensional Gaussian mixture with scalar input point
    function getLogProbGausMixSDSP_CK(nmode,nd,np,LogAmplitude_CK,MeanVec_CK,InvCovMat_CK,LogSqrtDetInvCovMat_CK,point_CK)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbGausMixSDSP_CK
#endif
        use Constants_mod, only: IK, RK, CK, LOGTINY_RK
        implicit none
        integer(IK), intent(in) :: nmode,nd,np
        complex(CK), intent(in) :: LogAmplitude_CK(nmode),MeanVec_CK(nmode)
        complex(CK), intent(in) :: InvCovMat_CK(nmode),LogSqrtDetInvCovMat_CK(nmode)
        complex(CK), intent(in) :: point_CK
        complex(CK)             :: getLogProbGausMixSDSP_CK
        complex(CK)             :: normFac_CK,LogProb_CK(nmode)
        integer(IK)             :: imode
        do imode = 1, nmode
            LogProb_CK(imode) = LogAmplitude_CK(imode) + &
            getLogProbNormSP_CK(MeanVec_CK(imode),InvCovMat_CK(imode),LogSqrtDetInvCovMat_CK(imode),point_CK)
        end do
        normFac_CK = maxval(real(LogProb_CK))
        LogProb_CK = LogProb_CK - normFac_CK
        where(real(LogProb_CK)<LOGTINY_RK)
            LogProb_CK = 0._RK
        elsewhere
            LogProb_CK = exp(LogProb_CK)
        end where
        getLogProbGausMixSDSP_CK = normFac_CK + log(sum(LogProb_CK))
    end function getLogProbGausMixSDSP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getLogProbGausMixSDMP_CK(nmode,nd,np,LogAmplitude_CK,MeanVec_CK,InvCovMat_CK,LogSqrtDetInvCovMat_CK,Point_CK)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbGausMixSDMP_CK
#endif
        use Constants_mod, only: IK, RK, CK, LOGTINY_RK
        implicit none
        integer(IK), intent(in) :: nmode,nd,np
        complex(CK), intent(in) :: LogAmplitude_CK(nmode),MeanVec_CK(nmode)
        complex(CK), intent(in) :: InvCovMat_CK(nmode),LogSqrtDetInvCovMat_CK(nmode)
        complex(CK), intent(in) :: Point_CK(np)
        complex(CK)             :: getLogProbGausMixSDMP_CK(np)
        complex(CK)             :: normFac_CK(np),LogProb_CK(nmode,np)
        integer(IK)             :: imode, ip
        do imode = 1, nmode
            LogProb_CK(imode,1:np) = LogAmplitude_CK(imode) + &
            getLogProbNormMP_CK(np,MeanVec_CK(imode),InvCovMat_CK(imode),LogSqrtDetInvCovMat_CK(imode),Point_CK)
        end do
        normFac_CK = maxval(real(LogProb_CK),dim=1)
        do ip = 1,np
            LogProb_CK(1:nmode,ip) = LogProb_CK(1:nmode,ip) - normFac_CK(ip)
            do imode = 1,nmode
                if ( real(LogProb_CK(imode,ip)) < LOGTINY_RK ) then
                    LogProb_CK(imode,ip) = 0._RK
                else
                    LogProb_CK(imode,ip) = exp( LogProb_CK(imode,ip) )
                end if
            end do
            getLogProbGausMixSDMP_CK(ip) = normFac_CK(ip) + log(sum(LogProb_CK(1:nmode,ip)))
        end do
    end function getLogProbGausMixSDMP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function getLogProbGausMixMDSP_CK(nmode,nd,np,LogAmplitude_CK,MeanVec_CK,InvCovMat_CK,LogSqrtDetInvCovMat_CK,Point_CK)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbGausMixMDSP_CK
#endif
        use Constants_mod, only: IK, RK, CK, LOGTINY_RK
        implicit none
        integer(IK), intent(in) :: nmode,nd,np
        complex(CK), intent(in) :: LogAmplitude_CK(nmode), MeanVec_CK(nd,nmode)
        complex(CK), intent(in) :: InvCovMat_CK(nd,nd,nmode), LogSqrtDetInvCovMat_CK(nmode)
        complex(CK), intent(in) :: Point_CK(nd)
        complex(CK)             :: getLogProbGausMixMDSP_CK
        complex(CK)             :: normFac_CK,LogProb_CK(nmode)
        integer(IK)             :: imode
        do imode = 1, nmode
            LogProb_CK(imode) = LogAmplitude_CK(imode) + &
            getLogProbMVNSP_CK(nd,MeanVec_CK(1:nd,imode),InvCovMat_CK(1:nd,1:nd,imode)&
                                ,LogSqrtDetInvCovMat_CK(imode),Point_CK)
        end do
        normFac_CK = maxval(real(LogProb_CK))
        LogProb_CK = LogProb_CK - normFac_CK
        where(real(LogProb_CK)<LOGTINY_RK)
            LogProb_CK = 0._RK
        elsewhere
            LogProb_CK = exp(LogProb_CK)
        end where
        getLogProbGausMixMDSP_CK = normFac_CK + log(sum(LogProb_CK))
    end function getLogProbGausMixMDSP_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getLogProbGausMixMDMP_CK(nmode,nd,np,LogAmplitude_CK,MeanVec_CK,InvCovMat_CK,LogSqrtDetInvCovMat_CK,Point_CK)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbGausMixMDMP_CK
#endif
        use Constants_mod, only: IK, RK, CK, LOGTINY_RK
        implicit none
        integer(IK), intent(in) :: nmode,nd,np
        complex(CK), intent(in) :: LogAmplitude_CK(nmode),MeanVec_CK(nd,nmode)
        complex(CK), intent(in) :: InvCovMat_CK(nd,nd,nmode), LogSqrtDetInvCovMat_CK(nmode)
        complex(CK), intent(in) :: Point_CK(nd,np)
        complex(CK)             :: getLogProbGausMixMDMP_CK(np)
        complex(CK)             :: normFac_CK(np),LogProb_CK(nmode,np)
        integer(IK)             :: imode, ip
        do imode = 1, nmode
            LogProb_CK(imode,1:np) = LogAmplitude_CK(imode) + &
            getLogProbMVNMP_CK(nd,np,MeanVec_CK(1:nd,imode),InvCovMat_CK(1:nd,1:nd,imode)&
                                ,LogSqrtDetInvCovMat_CK(imode),Point_CK)
        end do
        normFac_CK = maxval(real(LogProb_CK),dim=1)
        do ip = 1,np
            LogProb_CK(1:nmode,ip) = LogProb_CK(1:nmode,ip) - normFac_CK(ip)
            do imode = 1,nmode
                if ( real(LogProb_CK(imode,ip))<LOGTINY_RK ) then
                    LogProb_CK(imode,ip) = 0._RK
                else
                    LogProb_CK(imode,ip) = exp( LogProb_CK(imode,ip) )
                end if
            end do
            getLogProbGausMixMDMP_CK(ip) = normFac_CK(ip) + log(sum(LogProb_CK(1:nmode,ip)))
        end do
    end function getLogProbGausMixMDMP_CK

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

    pure function getMean_2D(nd,np,Point,Weight) result(Mean)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getMean_2D
#endif
        ! the implementation for one-dimension is very concise and nice: mean = sum(Weight*Point) / sum(Weight)
        implicit none
        integer(IK), intent(in)             :: np,nd            ! np: number of observations, nd: number of variables for each observation
        real(RK)   , intent(in)             :: Point(nd,np)     ! Point is the data matrix
        integer(IK), intent(in), optional   :: Weight(nd,np)    ! sample weight
        real(RK)                            :: Mean(nd)         ! output mean vector
        integer(IK)                         :: ip, SumWeight(nd)
        Mean = 0._RK
        if (present(Weight)) then
            SumWeight = 0._IK
            do ip = 1,np
                SumWeight = SumWeight + Weight(1:nd,ip)
                Mean = Mean + Weight(1:nd,ip) * Point(1:nd,ip)
            end do
            Mean = Mean / real(SumWeight,kind=RK)
        else
            do ip = 1,np
                Mean = Mean + Point(1:nd,ip)
            end do
            Mean = Mean / real(np,kind=RK)
        end if
    end function getMean_2D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getNormData_2D(nd,np,Mean,Point) result(NormData)
#if defined DLL_ENABLED && !defined CFI_ENABLED
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

    pure function getVariance_1D(np,mean,Point,Weight,sumWeight) result(variance)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getVariance_1D
#endif
        ! ATTN: if Weight is present, then SumWeight must be also present
        ! reason for requesting sumWeight: if mean is already given, it means sum weight is also computed a priori
        ! This function is rather redundant as one could also use the concise Fortran syntax to achieve the same goal:
        ! mean = sum(Weight*Point) / sum(Weight)
        ! variance = sum( (Weight*(Point-mean))**2 ) / (sum(Weight)-1)
        ! But the concise version will be slightly slower as it involves three loops instead of two.
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

    pure function getVariance_2D(nd,np,Mean,Point,Weight) result(Variance)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getVariance_2D
#endif
        ! returns the variance of each row in Point weighted by the corresponding Weight if provided
        ! pass the Mean vector by calling getMean() or getMean_2D()
        implicit none
        integer(IK), intent(in)             :: nd, np           ! np is the number of observations (points) whose variance is to be computed
        real(RK)   , intent(in)             :: Mean(nd)         ! the Mean value of the vector Point
        real(RK)   , intent(in)             :: Point(nd,np)     ! Point is the vector of data
        integer(IK), intent(in), optional   :: Weight(nd,np)    ! sample weight
        real(RK)                            :: Variance(nd)     ! output Mean vector
        integer(IK)                         :: ip, SumWeight(nd)
        Variance = 0._RK
        if (present(Weight)) then
            SumWeight = 0._IK
            do ip = 1,np
                SumWeight = SumWeight + Weight(1:nd,ip)
                Variance = Variance + Weight(1:nd,ip) * ( Point(1:nd,ip) - Mean )**2
            end do
            Variance = Variance / real(SumWeight-1_IK,kind=RK)
        else
            do ip = 1,np
                Variance = Variance + ( Point(1:nd,ip) - Mean )**2
            end do
            Variance = Variance / real(np-1_IK,kind=RK)
        end if
    end function getVariance_2D

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Returns the lower triangle Cholesky Factor of the covariance matrix of a set of points in the lower part of CholeskyLower.
    ! The upper part of CholeskyLower, including the diagonal elements of it, contains the Covariance matrix of the sample.
    ! Diagonal, contains on output, the diagonal terms of Cholesky Factor.
    subroutine getSamCholFac(nd,np,Mean,Point,CholeskyLower,Diagonal)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getSamCholFac
#endif
        use Matrix_mod, only: getCholeskyFactor
        implicit none
        integer(IK), intent(in)  :: nd,np                  ! np is the number of observations, nd is the number of parameters for each observation
        real(RK)   , intent(in)  :: Mean(nd)               ! Mean vector
        real(RK)   , intent(in)  :: Point(nd,np)           ! Point is the matrix of the data, CovMat contains the elements of the sample covariance matrix
        real(RK)   , intent(out) :: CholeskyLower(nd,nd)   ! Lower Cholesky Factor of the covariance matrix
        real(RK)   , intent(out) :: Diagonal(nd)           ! Diagonal elements of the Cholesky factorization
        real(RK)                 :: NormedData(np,nd), npMinusOneInverse
        integer(IK)              :: i,j

        do i = 1,np
            NormedData(i,1:nd) = Point(1:nd,i) - Mean
        end do

        ! Only upper half of CovMat is needed
        npMinusOneInverse = 1._RK / real(np-1,kind=RK)
        do j = 1,nd
            do i = 1,j
                ! Get the covariance matrix elements: only the upper half of CovMat is needed
                CholeskyLower(i,j) = dot_product( NormedData(1:np,i) , NormedData(1:np,j) ) * npMinusOneInverse
            end do
        end do

        ! XXX: The efficiency can be improved by merging it with the above loops
        call getCholeskyFactor(nd,CholeskyLower,Diagonal)

  end subroutine getSamCholFac

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This subroutine finds the elements of the symmetric nd*nd sample covariance matrix of a given set of
    ! np observations, each one with nd parameters in the format of a "" np*nd "" matrix.
    ! For a review, refer to Geisser & Cornfield (1963) "Posterior distributions for multivariate normal parameters".
    ! Also refer to Box and Tiao (1973), "Bayesian Inference in Statistical Analysis" Page 421.
    ! Amir Shahmoradi, Oct 16, 2009, 11:14 AM, MTU
    subroutine getSamCovMean(np,nd,Point,CovMat,Mean,MahalSq,InvCovMat,sqrtDetInvCovMat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getSamCovMean
#endif
        use Matrix_mod, only: getInvPosDefMatSqrtDet
        implicit none
        integer(IK), intent(in)            :: np,nd                  ! np is the number of observations, nd is the number of parameters for each observation
        real(RK)   , intent(in)            :: Point(np,nd)           ! Point is the matrix of the data, CovMat contains the elements of the sample covariance matrix
        real(RK)   , intent(out)           :: CovMat(nd,nd)          ! Covariance matrix of the input data
        real(RK)   , intent(out)           :: Mean(nd)               ! Mean vector
        real(RK)   , intent(out), optional :: MahalSq(np)            ! Vector of Mahalanobis Distances Squared, with respect to the mean position of the sample
        real(RK)   , intent(out), optional :: InvCovMat(nd,nd)       ! Inverse Covariance matrix of the input data
        real(RK)   , intent(out), optional :: sqrtDetInvCovMat       ! sqrt determinant of the inverse covariance matrix
        real(RK)   , dimension(nd)         :: DummyVec
        real(RK)   , dimension(np,nd)      :: NormedData
        integer(IK)                        :: i,j

        do j = 1,nd
            Mean(j) = sum(Point(1:np,j)) / real(np,kind=RK)
            NormedData(1:np,j) = Point(1:np,j) - Mean(j)
        end do
        do i = 1,nd
            do j = 1,nd
                CovMat(i,j) = dot_product(NormedData(1:np,i),NormedData(1:np,j))/real(np-1,kind=RK)
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
                DummyVec(j) = dot_product(InvCovMat(1:nd,j),NormedData(i,1:nd))
                end do
                MahalSq(i) = dot_product(NormedData(i,1:nd),DummyVec)
            end do
        end if

    end subroutine getSamCovMean

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This subroutine is exact same thing as getSamCovMean, with the only difference that input data is transposed here on input.
    ! Based on my preliminary benchmarks with Intel 2017 ifort, getSamCovMean() is slightly faster than getSamCovMeanTrans()
    subroutine getSamCovMeanTrans(np,nd,Point,CovMat,Mean,MahalSq,InvCovMat,sqrtDetInvCovMat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
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
        real(RK)   , dimension(nd,np)      :: NormedData
        integer(IK)                        :: i,j

        Mean = 0._RK
        do i = 1,np
            do j = 1,nd
                Mean(j) = Mean(j) + Point(j,i)
            end do
        end do
        Mean = Mean / real(np,kind=RK)

        do i = 1,np
            NormedData(1:nd,i) = Point(1:nd,i) - Mean
        end do

        do i = 1,nd
            do j = 1,nd
                CovMat(i,j) = dot_product(NormedData(i,1:np),NormedData(j,1:np)) / real(np-1,kind=RK)
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
                DummyVec(j) = dot_product(InvCovMat(1:nd,j),NormedData(1:nd,i))
                end do
                MahalSq(i) = dot_product(NormedData(1:nd,i),DummyVec)
                !MahalSq = dot_product(NormedData(1:nd,i),DummyVec)
                !if (maxMahal<MahalSq) maxMahal = MahalSq
            end do
            !maxMahal = maxval(MahalSq)
        end if

    end subroutine getSamCovMeanTrans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This subroutine is exact same thing as getSamCovMeanTrans, with the only differences that only the upper triangle of the covariance matrix is returned.
    ! Also, optional arguments are not available. This subroutine is specifically optimized for use in ParaMCMC.
    subroutine getSamCovUpperMeanTrans(np,nd,Point,CovMatUpper,Mean)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getSamCovUpperMeanTrans
#endif
        use Matrix_mod, only: getInvPosDefMatSqrtDet
        implicit none
        integer(IK), intent(in)            :: np,nd                  ! np is the number of observations, nd is the number of parameters for each observation
        real(RK)   , intent(in)            :: Point(nd,np)           ! Point is the matrix of the data, CovMatUpper contains the elements of the sample covariance matrix
        real(RK)   , intent(out)           :: CovMatUpper(nd,nd)     ! Covariance matrix of the input data
        real(RK)   , intent(out)           :: Mean(nd)               ! Mean vector
        real(RK)                           :: npMinusOneInvReal, NormedData(nd,np)
        integer(IK)                        :: i,j

        Mean = 0._RK
        do i = 1,np
            do j = 1,nd
                Mean(j) = Mean(j) + Point(j,i)
            end do
        end do
        Mean = Mean / real(np,kind=RK)

        do i = 1,np
            NormedData(1:nd,i) = Point(1:nd,i) - Mean
        end do

        npMinusOneInvReal = 1._RK / real(np-1,kind=RK)
        do j = 1,nd
            do i = 1,j
                CovMatUpper(i,j) = dot_product(NormedData(i,1:np),NormedData(j,1:np)) * npMinusOneInvReal
            end do
        end do

  end subroutine getSamCovUpperMeanTrans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This subroutine is exact same thing as getSamCovUpperMeanTrans,
    ! with the only difference that here points can have different weights. This subroutine is specifically optimized for use in ParaMCMC.
    subroutine getWeiSamCovUppMeanTrans(np,sumWeight,nd,Point,Weight,CovMatUpper,Mean)
#if defined DLL_ENABLED && !defined CFI_ENABLED
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
        real(RK)                            :: NormedData(nd,np)
        integer(IK)                         :: i,j,ip

        Mean = 0._RK
        do i = 1,np
            do j = 1,nd
                Mean(j) = Mean(j) + Weight(i)*Point(j,i)
            end do
        end do
        Mean = Mean / real(sumWeight,kind=RK)

        do i = 1,np
            NormedData(1:nd,i) = Point(1:nd,i) - Mean
        end do

        sumWeightMinusOneInvReal = 1._RK / real(sumWeight-1,kind=RK)
        do j = 1,nd
            do i = 1,j
                CovMatUpper(i,j) = 0
                do ip = 1,np
                    CovMatUpper(i,j) = CovMatUpper(i,j) + Weight(ip)*NormedData(i,ip)*NormedData(j,ip)
                end do
                CovMatUpper(i,j) = CovMatUpper(i,j) * sumWeightMinusOneInvReal
            end do
        end do

    end subroutine getWeiSamCovUppMeanTrans

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This subroutine combines the mean and covariance of two samples given as input.
    ! It uses the recursion equation similar to http://stats.stackexchange.com/questions/97642/how-to-combine-sample-means-and-sample-variances
    ! See also https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Covariance to update the covariance:
    ! See also https://stats.stackexchange.com/questions/43159/how-to-calculate-pooled-variance-of-two-groups-given-known-group-variances-mean
    subroutine combineCovMean(nd,npA,MeanA,CovMatA,npB,MeanB,CovMatB,Mean,CovMat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: combineCovMean
#endif

        implicit none
        integer(IK), intent(in)       :: nd
        integer(IK), intent(in)       :: npA,npB
        real(RK)   , intent(in)       :: MeanA(nd),CovMatA(nd,nd)
        real(RK)   , intent(in)       :: MeanB(nd),CovMatB(nd,nd)
        real(RK)   , intent(out)      :: Mean(nd),CovMat(nd,nd)
        real(RK)   , dimension(nd,1)  :: MatMeanA,MatMeanB,MatMean
        real(RK)                      :: npAB

        npAB = real(npA + npB,kind=RK)
        MatMeanA(1:nd,1) = MeanA
        MatMeanB(1:nd,1) = MeanB

        ! First find Mean:
        Mean = ( real(npA,kind=RK)*MeanA + real(npB,kind=RK)*MeanB ) / npAB
        MatMean(1:nd,1) = Mean

        ! Now find new Covariance matrix:
        CovMat  = real(npA,kind=RK) * ( CovMatA + matmul(MatMeanA,transpose(MatMeanA)) ) &
                + real(npB,kind=RK) * ( CovMatB + matmul(MatMeanB,transpose(MatMeanB)) )
        CovMat  = CovMat/npAB - matmul(MatMean,transpose(MatMean))

    end subroutine combineCovMean

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This subroutine is the same as combineCovMean, with the IMPORTANT difference that
    ! only the upper triangles and diagonals of the input covariance matrices need to be given by the user: CovMatUpper, CovMatUpperB
    subroutine combineMeanCovUpper(nd,npA,MeanVecA,CovMatUpperA,npB,MeanVecB,CovMatUpperB,MeanVec,CovMatUpper)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: combineMeanCovUpper
#endif

        implicit none
        integer(IK), intent(in)  :: nd
        integer(IK), intent(in)  :: npA,npB
        real(RK)   , intent(in)  :: MeanVecA(nd),CovMatUpperA(nd,nd)
        real(RK)   , intent(in)  :: MeanVecB(nd),CovMatUpperB(nd,nd)
        real(RK)   , intent(out) :: MeanVec(nd),CovMatUpper(nd,nd)
        real(RK)                 :: npAreal, npBreal, npABinverse
        integer(IK)              :: i,j

        npAreal = real(npA,kind=RK)
        npBreal = real(npB,kind=RK)
        npABinverse = 1._RK / real(npA + npB,kind=RK)

        ! First find MeanVec:
        MeanVec = npABinverse * ( npAreal*MeanVecA + npBreal*MeanVecB )

        ! Now find new Covariance matrix:
        do j = 1,nd
          do i = 1,j
            CovMatUpper(i,j) = ( npAreal * ( CovMatUpperA(i,j) + MeanVecA(i) * MeanVecA(j) ) &
                               + npBreal * ( CovMatUpperB(i,j) + MeanVecB(i) * MeanVecB(j) ) &
                               ) * npABinverse - MeanVec(i) * MeanVec(j)
            end do
        end do

    end subroutine combineMeanCovUpper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This subroutine is the same as combineMeanCovUpper, with the IMPORTANT difference that the resulting Mean and CovMat
    ! are now returned as CovMatB and MeanB.
    subroutine mergeMeanCovUpper(nd,npA,MeanVecA,CovMatUpperA,npB,MeanVecB,CovMatUpperB)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: mergeMeanCovUpper
#endif

        implicit none
        integer(IK), intent(in)     :: nd
        integer(IK), intent(in)     :: npA,npB
        real(RK)   , intent(in)     :: MeanVecA(nd),CovMatUpperA(nd,nd)
        real(RK)   , intent(inout)  :: MeanVecB(nd),CovMatUpperB(nd,nd)
        real(RK)                    :: MeanVec(nd)
        real(RK)                    :: npABinverse, npAnpABinverseProd, npBnpABinverseProd
        integer(IK)                 :: i,j

        npABinverse = 1._RK / real(npA + npB,kind=RK)
        npAnpABinverseProd = real(npA,kind=RK) * npABinverse
        npBnpABinverseProd = real(npB,kind=RK) * npABinverse

        do j = 1,nd
            MeanVec(j) = npAnpABinverseProd*MeanVecA(j) + npBnpABinverseProd*MeanVecB(j)
            do i = 1,j
                CovMatUpperB(i,j)   = npAnpABinverseProd * ( CovMatUpperA(i,j) + MeanVecA(i) * MeanVecA(j) ) &
                                    + npBnpABinverseProd * ( CovMatUpperB(i,j) + MeanVecB(i) * MeanVecB(j) ) &
                                    - MeanVec(i) * MeanVec(j)
            end do
        end do
        MeanVecB = MeanVec

    end subroutine mergeMeanCovUpper

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This code returns a normally distributed deviate with zero mean and unit variance.
    function getRandGaus()
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandGaus
#endif

        implicit none
        integer(IK), save :: iset=0
        real(RK)   , save :: gset
        real(RK)          :: fac,rsq,vec(2)
        real(RK)          :: getRandGaus

        if (iset == 0) then
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
            getRandGaus = vec(2)*fac
            iset = 1
        else
            getRandGaus = gset
            iset = 0
        endif

    end function getRandGaus

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This code returns a normally distributed deviate with the given mean and standard deviation.
    function getRandNorm(mean,std)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandNorm
#endif
        implicit none
        real(RK), intent(in) :: mean, std
        real(RK)       :: getRandNorm
        getRandNorm = mean + std*getRandGaus()
    end function getRandNorm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This code returns a log-normally distributed deviate with the given mean and standard deviation.
    function getRandLogn(mean,std)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandLogn
#endif
        implicit none
        real(RK), intent(in) :: mean, std
        real(RK)       :: getRandLogn
        getRandLogn = exp( mean + std*getRandGaus() )
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
    ! Amir Shahmoradi, March 22, 2012, 2:21 PM, IFS, UTEXAS
    subroutine getMVNDev(nd,MeanVec,CovMatIn,X)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getMVNDev
#endif

        use Matrix_mod, only: getCholeskyFactor

        implicit none
        integer(IK), intent(in)  :: nd
        real(RK)   , intent(in)  :: MeanVec(nd), CovMatIn(nd,nd)
        real(RK)   , intent(out) :: X(nd)
        real(RK)                 :: CovMat(nd,nd), Diagonal(nd), DummyVec(nd)
        integer(IK)              :: i

        CovMat = CovMatIn
        call getCholeskyFactor(nd,CovMat,Diagonal)
        if (Diagonal(1)<0._RK) then
            write(*,*) 'getCholeskyFactor() failed in getMVNDev()'
            stop
        end if
        do i=1,nd
            DummyVec(i) = getRandGaus()
            x(i) = DummyVec(i) * Diagonal(i)
        end do
        do i=2,nd
            x(i) = x(i) + dot_product(CovMat(i,1:i-1),DummyVec(1:i-1))
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
    ! Amir Shahmoradi, April 25, 2016, 2:21 PM, IFS, UTEXAS
    subroutine getMVUDev(nd,MeanVec,CovMatIn,X)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getMVUDev
#endif

        use Matrix_mod, only: getCholeskyFactor

        implicit none
        integer(IK), intent(in)  :: nd
        real(RK)   , intent(in)  :: MeanVec(nd)
        real(RK)   , intent(in)  :: CovMatIn(nd,nd)
        real(RK)   , intent(out) :: X(nd)
        real(RK)                 :: Diagonal(nd), DummyVec(nd), CovMat(nd,nd), dummy
        integer(IK)              :: i

        CovMat = CovMatIn
        call getCholeskyFactor(nd,CovMat,Diagonal)
        if (Diagonal(1)<0._RK) then
            error stop
            !call abortProgram( output_unit , 1 , 1 , 'Statitistics@getMVUDev()@getCholeskyFactor() failed.' )
        end if
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
            X(i) = X(i) + dot_product(CovMat(i,1:i-1),DummyVec(1:i-1))
        end do

        X = X + MeanVec

    end subroutine getMVUDev

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getRandMVN(nd,MeanVec,CholeskyLower,Diagonal) result(RandMVN)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandMVN
#endif
        ! Amir Shahmoradi, April 23, 2017, 12:36 AM, ICES, UTEXAS
        ! returns a multivariate Normal random vector.
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: MeanVec(nd)
        real(RK)   , intent(in) :: CholeskyLower(nd,nd), Diagonal(nd)   ! Cholesky lower triangle and its diagonal terms, calculated from the input CovMat.
        real(RK)                :: RandMVN(nd), dummy
        integer(IK)             :: j,i
        RandMVN = 0._RK
        do j = 1,nd
            dummy = getRandGaus()
            RandMVN(j) = RandMVN(j) + Diagonal(j) * dummy
            do i = j+1,nd
                RandMVN(i) = RandMVN(i) + CholeskyLower(i,j) * dummy
            end do
        end do
        RandMVN = RandMVN + MeanVec
    end function getRandMVN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!    ! Given an input Mean vector of length nd, Covariance Matrix of dimension (nd,nd), as well as a vector of integers representing
!    ! the variables (corresponding to CovMat columns) that are given
!    ! This subroutine gives out a conditional Multivariate Normal Random deviate.
!    ! random p-tivariate normal deviate, given that the first pg variables x1 are given (i.e. fixed).
!    ! For a review of Multivariate Normal distribution: Applied Multivariate Statistical Analysis, Johnson, Wichern, 1998, 4th ed.
!    ! Amir Shahmoradi, Oct 20, 2009, 9:12 PM, MTU
!    function getCondRandMVN(nd,MeanVec,CovMat,nIndIndx,IndIndx) result(CondRandMVN)
!        use Matrix_mod, only: getRegresCoef
!        implicit none
!        integer(IK), intent(in) :: nd, nIndIndx, IndIndx(nIndIndx)
!        real(RK)   , intent(in) :: MeanVec(nd), CovMat(nd,nd)
!        real(RK)                :: CondRandMVN(nd),
!        integer(IK)             :: j, i
!    end function getCondRandMVN

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Amir Shahmoradi, April 23, 2017, 1:36 AM, ICES, UTEXAS
    ! Given the Cholesky Lower triangle and diagonals of a given covariance matrix, this function return one point
    ! uniformly randomly drawn from inside of an nd-ellipsoid, whose nd elements getRandMVU(i), i=1,nd are guaranteed to be in range
    ! MeanVec(i) - sqrt(CovMat(i,i)) < getRandMVU(i) < MeanVec(i) + sqrt(CovMat(i,i))
    function getRandMVU(nd,MeanVec,CholeskyLower,Diagonal)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandMVU
#endif
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: MeanVec(nd)
        real(RK)   , intent(in) :: CholeskyLower(nd,nd) ! Cholesky lower triangle, calculated from the input CovMat.
        real(RK)   , intent(in) :: Diagonal(nd)         ! Cholesky diagonal terms, calculated from the input CovMat.
        real(RK)                :: getRandMVU(nd), dummy, DummyVec(nd), sumSqDummyVec
        integer(IK)             :: i,j
        sumSqDummyVec = 0._RK
        do j=1,nd
            DummyVec(j) = getRandGaus()
            sumSqDummyVec = sumSqDummyVec + DummyVec(j)**2
        end do
        call random_number(dummy)
        dummy = (dummy**(1._RK/dble(nd)))/sqrt(sumSqDummyVec)
        DummyVec = DummyVec * dummy  ! DummyVec(j) * dummy is a uniform random point from inside of nd-sphere
        getRandMVU = 0._RK
        do j = 1,nd
            getRandMVU(j) = getRandMVU(j) + Diagonal(j) * DummyVec(j)
            do i = j+1,nd
                getRandMVU(i) = getRandMVU(i) + CholeskyLower(i,j) * DummyVec(j)
            end do
        end do
        getRandMVU = getRandMVU + MeanVec
    end function getRandMVU

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! determines if the input NormedPoint (normalized with respect to the mean of the target ellipsoid) is within or on the boundary 
    ! of the ellipsoid that is centered at origin and its boundary is described by the representative matrix Sigma (RepMat), 
    ! such that x^T Sigma^(-1) x = 1 for all x on the boundary. Note that the input matrix is the inverse of RepMat: InvRepMat
    pure function isInsideEllipsoid(nd,NormedPoint,InvRepMat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: isInsideEllipsoid
#endif
        use Math_mod, only: getLogVolEllipsoid
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: NormedPoint(nd)
        real(RK)   , intent(in) :: InvRepMat(nd,nd)
        logical                 :: isInsideEllipsoid
        isInsideEllipsoid = dot_product(NormedPoint,matmul(InvRepMat,NormedPoint)) <= 1._RK
    end function isInsideEllipsoid

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! NOTE: if MahalSq computation fails, output probability will be returned as NullVal%RK from module Constant_mod.
    pure function getLogProbMVU(nd,logSqrtDetInvCovMat) result(logProbMVU)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbMVU
#endif
        use Math_mod, only: getLogVolEllipsoid
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: logSqrtDetInvCovMat
        real(RK)                :: logProbMVU
        logProbMVU = -getLogVolEllipsoid(nd,logSqrtDetInvCovMat)
    end function getLogProbMVU

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Amir Shahmoradi, April 23, 2017, 1:36 AM, ICES, UTEXAS
    ! This is algorithm is similar to getRandMVU, with the only difference that points are drawn randomly from the surface of the ellipsoid instead of inside of its interior.
    ! Note that the distribution of points on the surface of the ellipsoid is NOT uniform.
    ! Regions of high curvature will have more points randomly sampled from them.
    ! generating uniform random points on arbitrary-dimension ellipsoids is not a task with trivial solution!
    function getRandPointOnEllipsoid(nd,MeanVec,CholeskyLower,Diagonal)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandPointOnEllipsoid
#endif
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: MeanVec(nd)
        real(RK)   , intent(in) :: CholeskyLower(nd,nd) ! Cholesky lower triangle, calculated from the MVN CovMat.
        real(RK)   , intent(in) :: Diagonal(nd)         ! Cholesky diagonal terms, calculated from the MVN CovMat.
        real(RK)                :: getRandPointOnEllipsoid(nd), DummyVec(nd), sumSqDummyVec
        integer(IK)             :: i,j
        sumSqDummyVec = 0._RK
        do j=1,nd
            DummyVec(j) = getRandGaus()
            sumSqDummyVec = sumSqDummyVec + DummyVec(j)**2
        end do
        DummyVec = DummyVec / sqrt(sumSqDummyVec)    ! DummyVec is a random point on the surface of nd-sphere.
        getRandPointOnEllipsoid = 0._RK
        do j = 1,nd
            getRandPointOnEllipsoid(j) = getRandPointOnEllipsoid(j) + Diagonal(j) * DummyVec(j)
            do i = j+1,nd
                getRandPointOnEllipsoid(i) = getRandPointOnEllipsoid(i) + CholeskyLower(i,j) * DummyVec(j)
            end do
        end do
        getRandPointOnEllipsoid = getRandPointOnEllipsoid + MeanVec
    end function getRandPointOnEllipsoid

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns the complementary error function with fractional error everywhere less than 1.2*10^-7.
    ! amir shahmoradi, Monday March 6, 2017, 3:22 pm, ICES, The University of Texas at Austin.
    pure function erfcc(x)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: erfcc
#endif
        implicit none
        real(RK), intent(in) :: x
        real(RK)             :: erfcc
        real(RK)             :: t,z
        z = abs(x)
        t = 1._RK/(1._RK+0.5_RK*z)
        erfcc = t*exp(-z*z-1.26551223_RK+t*(1.00002368_RK+t*(.37409196_RK+t*&
                (.09678418_RK+t*(-.18628806_RK+t*(.27886807_RK+t*(-1.13520398_RK+t*&
                (1.48851587_RK+t*(-.82215223_RK+t*.17087277_RK)))))))))
        if (x < 0._RK) erfcc = 2._RK - erfcc
    end function erfcc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogProbLogNormS(logMean,inverseVariance,logSqrtInverseVariance,logPoint)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbLogNormS
#endif
        use Constants_mod, only: LOGINVSQRT2PI
        implicit none
        real(RK), intent(in) :: logMean,inverseVariance,logSqrtInverseVariance,logPoint
        real(RK)             :: getLogProbLogNormS
        getLogProbLogNormS = LOGINVSQRT2PI + logSqrtInverseVariance - logPoint - 0.5_RK * inverseVariance * (logPoint-logMean)**2
    end function getLogProbLogNormS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function GetLogProbLogNormMP(np,logMean,inverseVariance,logSqrtInverseVariance,LogPoint)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: GetLogProbLogNormMP
#endif
        use Constants_mod, only: LOGINVSQRT2PI
        implicit none
        integer(IK), intent(in) :: np
        real(RK)   , intent(in) :: logMean,inverseVariance,logSqrtInverseVariance,LogPoint(np)
        real(RK)                :: GetLogProbLogNormMP(np)
        GetLogProbLogNormMP = LOGINVSQRT2PI + logSqrtInverseVariance - LogPoint - 0.5_RK * inverseVariance * (LogPoint-logMean)**2
    end function GetLogProbLogNormMP

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ! Returns random numbers with priodicity larger than ~2*10**18 random numbers.
  ! For more info see Press et al. (1990) Numerical Recipes.
    function getRandRealLecuyer(idum)  ! do not change idum between calls
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandRealLecuyer
#endif
        implicit none
        integer(IK), intent(inout) :: idum
        integer(IK), parameter     :: im1=2147483563, im2=2147483399, imm1=im1-1, ia1=40014, ia2=40692
        integer(IK), parameter     :: iq1=53668, iq2=52774, ir1=12211, ir2=3791, ntab=32, ndiv=1+imm1/ntab
        real(RK)   , parameter     :: am=1._RK/real(im1,kind=RK), eps=1.2e-7_RK, rnmx=1._RK-eps
        real(RK)                   :: getRandRealLecuyer
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
        getRandRealLecuyer = min(am*iy,rnmx)
    end function getRandRealLecuyer

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns integer random number in the range [lowerBound,upperBound],
    ! using the real random number generator of Lecuyer: Statistics@getRandRealLecuyer()
    function getRandIntLecuyer(lowerBound,upperBound,idum)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandIntLecuyer
#endif
        implicit none
        integer(IK), intent(in)    :: lowerBound,upperBound
        integer(IK), intent(inout) :: idum
        integer(IK)                :: getRandIntLecuyer
        getRandIntLecuyer = lowerBound + nint( getRandRealLecuyer(idum)*real(upperBound-lowerBound,kind=RK) )
    end function getRandIntLecuyer

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns integer random number in the range [lowerBound,upperBound], using built-in random number generator of Fortran.
    function getRandInt(lowerBound,upperBound)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandInt
#endif
        implicit none
        integer(IK), intent(in) :: lowerBound,upperBound
        real(RK)                :: dummy
        integer(IK)             :: getRandInt
        call random_number(dummy)
        getRandInt = lowerBound + nint( dummy*real(upperBound-lowerBound,kind=RK) )
    end function getRandInt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns real random number in the range [lowerBound,upperBound], using built-in random number generator of Fortran.
    function getRandUniform(lowerBound,upperBound) result(unifrnd)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandUniform
#endif
        implicit none
        real(RK), intent(in)    :: lowerBound,upperBound
        real(RK)                :: unifrnd
        call random_number(unifrnd)
        unifrnd = lowerBound + unifrnd * real(upperBound-lowerBound,kind=RK)
    end function getRandUniform

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ATTN: alpha must be > 0, else negative random Gamma variable (-1) will be output to indicate error occurrence.
    ! returns a Gamma distributed random number, following the prescription in GSL library.
    function getRandGamma(alpha)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandGamma
#endif
        implicit none
        real(RK), intent(in) :: alpha
        real(RK)             :: getRandGamma
        real(RK)             :: c,u,v,z
        if (alpha<=0._RK) then  ! illegal value of alpha
            getRandGamma = -1._RK
            return
        else
            getRandGamma = alpha
            if (getRandGamma<1._RK) getRandGamma = getRandGamma + 1._RK
            getRandGamma = getRandGamma - 0.3333333333333333
            c = 3._RK*sqrt(getRandGamma)
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
                if ( log(u)>=0.5_RK*z**2+getRandGamma*(1._RK-v+log(v)) ) cycle
                getRandGamma = getRandGamma*v
                exit
            end do
            if (alpha<1._RK) then
                call random_number(u)
                getRandGamma = getRandGamma * u**(1._RK/alpha)
            end if
        end if
    end function getRandGamma

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ATTN: alpha must be > 1, else a negative random Gamma value (-1) will be output, to indicate error occurrence.
    ! returns a Gamma-distributed random number whose shape parameter is integer
    function getRandGammaIntShape(alpha)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandGammaIntShape
#endif
        implicit none
        integer(IK), intent(in) :: alpha
        real(RK)                :: getRandGammaIntShape
        real(RK)                :: am,e,h,s,x,y,Vector(2),Array(5)
        if (alpha < 1) then  ! illegal value of alpha
            getRandGammaIntShape = -1._RK
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
        getRandGammaIntShape = x
    end function getRandGammaIntShape

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ATTN: both alpha, beta must be > 1, else a negative random Beta variable will be output to indicate error.
    ! returns a random BEta-distributed variable.
    function getRandBeta(alpha,beta)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandBeta
#endif
        implicit none
        real(RK), intent(in) :: alpha,beta
        real(RK)             :: getRandBeta
        real(RK)             :: x
        if ( alpha>0._RK .and. beta>0._RK ) then  ! illegal value of alpha
            x = getRandGamma(alpha)
            getRandBeta = x / ( x + getRandGamma(beta) )
        else
            getRandBeta = -1._RK
        end if
    end function getRandBeta

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns an exponentially-distributed random variable with mean 1/invMean
    function getRandExpWithInvMean(invMean) result(randExp)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandExpWithInvMean
#endif
        implicit none
        real(RK), intent(in)    :: invMean
        real(RK)                :: randExp
        call random_number(randExp)
        randExp = -log(randExp) * invMean
    end function getRandExpWithInvMean

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns a standard (mean=1) exponential-distributed random variable
    function getRandExp()
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandExp
#endif
        implicit none
        real(RK) :: getRandExp
        call random_number(getRandExp)
        getRandExp = -log(getRandExp)
    end function getRandExp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ATTN: nd>1 and eta>0.0, otherwise the first element of output, getRandCorMat(1,1), will be set to -1 to indicate error.
    ! returns a random correlation matrix
    function getRandCorMat(nd,eta)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandCorMat
#endif
        use Matrix_mod, only: getCholeskyFactor
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: eta
        real(RK)                :: getRandCorMat(nd,nd), dummy
        real(RK)                :: beta,sumSqDummyVec,DummyVec(nd-1),W(nd-1),Diagonal(nd-1)
        integer(IK)             :: m,j,i

        if (nd<2 .or. eta<=0._RK) then  ! illegal value for eta. set getRandCorMat=0, return
            getRandCorMat(1,1) = -1._RK
            return
        end if

        do m = 1,nd
            getRandCorMat(m,m) = 1._RK
        end do
        beta = eta + 0.5_RK*(nd-2._RK)
        dummy = getRandBeta(beta,beta)
        if (dummy<=0._RK .or. dummy>=1._RK) then
            error stop
            !call abortProgram( output_unit , 1 , 1 , 'Statitistics@getRandCorMat() failed. Random Beta variable out of bound: ' // num2str(dummy) )
        end if
        getRandCorMat(1,2) = 2._RK * dummy - 1._RK    ! for the moment, only the upper half of getRandCorMat is needed, the lower half will contain cholesky lower triangle.

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
            call getCholeskyFactor(m,getRandCorMat(1:m,1:m),Diagonal(1:m))
            if (Diagonal(1)<0._RK) then
                error stop
                !call abortProgram( output_unit , 1 , 1 , 'Statitistics@getRandCorMat()@getCholeskyFactor() failed.' )
            end if
            DummyVec(1:m) = 0._RK
            do j = 1,m
                DummyVec(j) = DummyVec(j) + Diagonal(j) * W(j)
                do i = j+1,m
                    DummyVec(i) = DummyVec(i) + getRandCorMat(i,j) * DummyVec(j)
                end do
            end do
            getRandCorMat(1:m,m+1) = DummyVec(1:m)
        end do
        do i=1,nd-1
            getRandCorMat(i+1:nd,i) = getRandCorMat(i,i+1:nd)
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

    ! Returns a random correlation matrix using Monte Carlo rejection method.
    ! Only the upper half of getRandCorMatRejection is the correlatrion matrix, lower half is giberish.
    ! Therefore this subroutine is very slow for high matrix dimensions ( nd>10 ).
    function getRandCorMatRejection(nd,minRho,maxRho)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandCorMatRejection
#endif
        use Matrix_mod, only: isPosDef
        implicit none
        integer(IK), intent(in) :: nd
        real(RK)   , intent(in) :: minRho,maxRho
        real(RK)                :: getRandCorMatRejection(nd,nd),RhoVec(nd*(nd-1))
        integer(IK)             :: i,j,irho
        if (maxRho<minRho .or. nd<1) then
            error stop
            !call abortProgram( output_unit , 1 , 1 , 'Statitistics@getRandCorMatRejection() failed: Invalid input values.' )
        end if
        if (nd==1) then
          getRandCorMatRejection = 1._RK
        else
            rejection: do
                call random_number(RhoVec)
                RhoVec = minRho + RhoVec*(maxRho-minRho)
                irho = 0
                do j=1,nd
                    getRandCorMatRejection(j,j) = 1._RK
                    do i=1,j-1
                        irho = irho + 1
                        getRandCorMatRejection(i,j) = RhoVec(irho)
                    end do
                end do
                if (isPosDef(nd,getRandCorMatRejection)) exit rejection
                cycle rejection
            end do rejection
        end if
        do j=1,nd-1
            getRandCorMatRejection(j+1:nd,j) = getRandCorMatRejection(j,j+1:nd)
        end do
  end function getRandCorMatRejection

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getUpperCorMatFromUpperCovMat(nd,CovMatUpper) result(CorMat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getUpperCorMatFromUpperCovMat
#endif
        implicit none
        integer(IK)  , intent(in) :: nd
        real(RK)     , intent(in) :: CovMatUpper(nd,nd)   ! only upper half needed
        real(RK)                  :: CorMat(nd,nd)
        real(RK)                  :: StdVec(nd)
        integer(IK)               :: i,j
        do j=1,nd
            StdVec(j) = sqrt(CovMatUpper(j,j))
            do i=1,j
                CorMat(i,j) = CovMatUpper(i,j) / ( StdVec(j) * StdVec(i) )
            end do
        end do
    end function getUpperCorMatFromUpperCovMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getUpperCovMatFromUpperCorMat(nd,StdVec,CorMat) result(CovMat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getUpperCovMatFromUpperCorMat
#endif
        implicit none
        integer(IK)  , intent(in) :: nd
        real(RK)     , intent(in) :: StdVec(nd), CorMat(nd,nd)   ! only upper half needed
        real(RK)                  :: CovMat(nd,nd)
        integer(IK)               :: i,j
        do j=1,nd
            CovMat(j,j) = StdVec(j)**2
            do i=1,j-1
                CovMat(i,j) = CorMat(i,j) * StdVec(j) * StdVec(i)
            end do
        end do
    end function getUpperCovMatFromUpperCorMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getUpperCovMatFromLowerCorMat(nd,StdVec,CorMat) result(CovMat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getUpperCovMatFromLowerCorMat
#endif
        implicit none
        integer(IK)  , intent(in) :: nd
        real(RK)     , intent(in) :: StdVec(nd), CorMat(nd,nd)   ! only upper half needed
        real(RK)                  :: CovMat(nd,nd)
        integer(IK)               :: i,j
        do j=1,nd
            CovMat(j,j) = StdVec(j)**2
            do i=1,j-1
                CovMat(i,j) = CorMat(j,i) * StdVec(j) * StdVec(i)
            end do
        end do
    end function getUpperCovMatFromLowerCorMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLowerCovMatFromUpperCorMat(nd,StdVec,CorMat) result(CovMat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLowerCovMatFromUpperCorMat
#endif
        implicit none
        integer(IK)  , intent(in) :: nd
        real(RK)     , intent(in) :: StdVec(nd), CorMat(nd,nd)   ! only upper half needed
        real(RK)                  :: CovMat(nd,nd)
        integer(IK)               :: i,j
        do j=1,nd
            CovMat(j,j) = StdVec(j)**2
            do i=1,j-1
                CovMat(j,i) = CorMat(i,j) * StdVec(j) * StdVec(i)
            end do
        end do
    end function getLowerCovMatFromUpperCorMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLowerCovMatFromLowerCorMat(nd,StdVec,CorMat) result(CovMat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLowerCovMatFromLowerCorMat
#endif
        implicit none
        integer(IK)  , intent(in) :: nd
        real(RK)     , intent(in) :: StdVec(nd), CorMat(nd,nd)   ! only upper half needed
        real(RK)                  :: CovMat(nd,nd)
        integer(IK)               :: i,j
        do j=1,nd
            CovMat(j,j) = StdVec(j)**2
            do i=1,j-1
                CovMat(j,i) = CorMat(j,i) * StdVec(j) * StdVec(i)
            end do
        end do
    end function getLowerCovMatFromLowerCorMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getCovMatFromCorMat(nd,StdVec,CorMat) result(CovMat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getCovMatFromCorMat
#endif
        implicit none
        integer(IK)  , intent(in) :: nd
        real(RK)     , intent(in) :: StdVec(nd), CorMat(nd,nd)   ! only upper half needed
        real(RK)                  :: CovMat(nd,nd)
        integer(IK)               :: i,j
        do j=1,nd
            CovMat(j,j) = StdVec(j)**2
            do i=1,j-1
                CovMat(i,j) = CorMat(i,j) * StdVec(j) * StdVec(i)
                CovMat(j,i) = CovMat(i,j)
            end do
        end do
    end function getCovMatFromCorMat

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getGeoLogPDF_old(successProb,logPdfPrecision,minSeqLen,seqLen) result(GeoLogPDF)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getGeoLogPDF_old
#endif
        !   return the Geometric distribution PDF values for a range of trials, starting at index 1.
        !   If the probability of success on each trial is successProb,
        !   then the probability that the kth trial (out of k trials) is
        !   the first success is GeoLogPDF(k).
        !
        !   Input
        !
        !       successProb
        !
        !           The probability of success.
        !
        !       logPdfPrecision
        !
        !           The precision value below which the PDF is practically considered to be zero.
        !
        !       minSeqLen
        !
        !           The minimum length of the range of k values for which the PDF will be computed.
        !
        !       seqLen
        !
        !           The length of the range of k values for which the PDF will be computed.
        !           If provided, it will overwrite the the output sequence length as inferred from 
        !           the combination of minSeqLen and logPdfPrecision.
        !
        !   Output
        !
        !       GeoLogPDF
        !
        !           An allocatable representing the geometric PDF over a range of k values, whose length is 
        !           seqLen, or if not provided, is determined from the values of logPdfPrecision and minSeqLen.
        !
        use Constants_mod, only: IK, RK
        implicit none
        real(RK)    , intent(in)            :: successProb
        real(RK)    , intent(in), optional  :: logPdfPrecision
        integer(IK) , intent(in), optional  :: minSeqLen
        integer(IK) , intent(in), optional  :: seqLen
        real(RK)    , allocatable           :: GeoLogPDF(:)
        real(RK)    , parameter             :: LOG_PDF_PRECISION = log(0.001_RK)
        real(RK)                            :: logProbFailure
        integer(IK)                         :: lenGeoLogPDF, i
        logProbFailure = log(1._RK - successProb)
        if (present(seqLen)) then
            lenGeoLogPDF = seqLen
        else
            if (present(logPdfPrecision)) then
                lenGeoLogPDF = ceiling(  logPdfPrecision / logProbFailure)
            else
                lenGeoLogPDF = ceiling(LOG_PDF_PRECISION / logProbFailure)
            end if
            if (present(minSeqLen)) lenGeoLogPDF = max(minSeqLen,lenGeoLogPDF)
        end if
        allocate(GeoLogPDF(lenGeoLogPDF))
        GeoLogPDF(1) = log(successProb)
        do i = 2, lenGeoLogPDF
            GeoLogPDF(i) = GeoLogPDF(i-1) + logProbFailure
        end do
    end function getGeoLogPDF_old

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogProbGeo(numTrial, SuccessStep, successProb) result(LogProbGeo)
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

    function fitGeoLogPDF_old(numTrial, SuccessStep, LogCount) result(PowellMinimum)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: fitGeoLogPDF_old
#endif
        !
        !   return a fit of the Geometric distribution PDF to the input natural logarithm of a sequence of Counts, 
        !   the ith element of which represents the number of successes after SuccessStep(i) tries in a Bernoulli trail.
        !   It is assumed that the input LogCount represents a sequence of counts at
        !   the first success is GeoLogPDF(k).
        !
        !   Input
        !
        !       numTrial
        !
        !           The length of the input vector LogCount
        !
        !       LogCount
        !
        !           A vector of real values representing the natural logarithms of the counts of success at each Bernoulli trial, 
        !           sequentially, from 1 to numTrial.
        !
        !   Output
        !
        !       PowellMinimum
        !
        !           An object of type PowellMinimum_type from Optimiziation_mod, containing the best-fit successProb and the 
        !           normalization constant of the fit in the vector component xmin.
        !
        use Optimization_mod, only: PowellMinimum_type
        use Constants_mod, only: IK, RK, POSINF_RK, NEGINF_RK
        implicit none
        integer(IK) , intent(in)    :: numTrial
        integer(IK) , intent(in)    :: SuccessStep(numTrial)
        real(RK)    , intent(in)    :: LogCount(numTrial)
        type(PowellMinimum_type)    :: PowellMinimum

        real(RK)                    :: BestFitSuccessProbNormFac(2) ! vector of the two parameters
        real(RK)    , parameter     :: SUCCESS_PROB_INIT_GUESS = 0.23_RK
        real(RK)    , parameter     :: FISHER_TRANS_SUCCESS_PROB_INIT_GUESS = atanh(2*(SUCCESS_PROB_INIT_GUESS - 0.5_RK))

        ! do Fisher transformation to make the limits infinity
        BestFitSuccessProbNormFac = [FISHER_TRANS_SUCCESS_PROB_INIT_GUESS, LogCount(1)]

        PowellMinimum = PowellMinimum_type  ( ndim = 2_IK &
                                            , getFuncMD = getSumDistSq &
                                            , StartVec = BestFitSuccessProbNormFac &
                                            )
        if (PowellMinimum%Err%occurred) return
        PowellMinimum%xmin(1) = 0.5_RK * tanh(PowellMinimum%xmin(1)) + 0.5_RK ! reverse Fisher-transform

    contains

        function getSumDistSq(ndim,successProbFisherTransNormFac) result(sumDistSq)
            !use Constants_mod, only: IK, RK
            implicit none
            integer(IK) , intent(in)    :: ndim
            real(RK)    , intent(in)    :: successProbFisherTransNormFac(ndim)
            real(RK)                    :: sumDistSq, successProb
            successProb = 0.5_RK*tanh(successProbFisherTransNormFac(1)) + 0.5_RK ! reverse Fisher-transform
            !sumDistSq = sum( (LogCount - getGeoLogPDF(successProb=successProb,seqLen=numTrial) - successProbFisherTransNormFac(2) )**2 )
            sumDistSq = sum(    ( LogCount &
                                - numTrial * successProbFisherTransNormFac(2) &
                                - getLogProbGeo(numTrial = numTrial, SuccessStep = SuccessStep, successProb = successProb) &
                                )**2 &
                            )
        end function getSumDistSq

    end function fitGeoLogPDF_old

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogProbGeoCyclic(successProb,maxNumTrial,numTrial,SuccessStep) result(LogProbGeoCyclic)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbGeoCyclic
#endif
        !   Compute the natural logarithm of the Geometric distribution PDF of a limited range of Bernoulli trials, 
        !   starting at index 1 up to maxNumTrial. In other words, upon reaching the trial maxNumTrial, 
        !   the Bernoulli trials count restart from index 1. This Cyclic Geometric distribution is 
        !   particularly useful in the parallelization studies of Monte Carlo simulation.
        !
        !   Input
        !
        !       successProb
        !
        !           The probability of success.
        !
        !       maxNumTrial
        !
        !           The maximum number of trails possible. After maxNumTrial tries, 
        !           the Geometric distribution restarts from index 1.
        !
        !       numTrial < maxNumTrial
        !
        !           The length of the array SuccessStep. 
        !
        !       SuccessStep(1:numTrial)
        !
        !           An array of integers that represent the steps at which the Bernoulli successes occur. 
        !
        !           WARNING: Any value of SuccessStep must be an integer numbers between 1 and maxNumTrial.
        !           WARNING: The onus is on the user to ensure this condition holds.
        !
        !   Output
        !
        !       LogProbGeoCyclic(1:numTrial)
        !
        !           A real-valued array of length numTrial whose values represent the probabilities of having 
        !           Bernoulli successes at the corresponding SuccessStep values.
        !
        use Constants_mod, only: IK, RK, NEGLOGINF_RK
        implicit none
        real(RK)    , intent(in)    :: successProb
        integer(IK) , intent(in)    :: maxNumTrial
        integer(IK) , intent(in)    :: numTrial
        integer(IK) , intent(in)    :: SuccessStep(numTrial)
        real(RK)                    :: LogProbGeoCyclic(numTrial)
        real(RK)                    :: failureProb, logProbSuccess, logProbFailure, logDenominator, exponentiation
        integer(IK)                 :: i
        if (successProb>0._RK .and. successProb<1._RK) then
            failureProb = 1._RK - successProb
            logProbSuccess = log(successProb)
            logProbFailure = log(failureProb)
            exponentiation = maxNumTrial * logProbFailure
            if (exponentiation<NEGLOGINF_RK) then
                logDenominator = 0._RK
            else
                logDenominator = log(1._RK-exp(exponentiation))
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

    function fitGeoCyclicLogPDF(maxNumTrial, numTrial, SuccessStep, LogCount) result(PowellMinimum)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: fitGeoLogPDF
#endif
        !
        !   Return a fit of the Cyclic Geometric distribution PDF to the input natural logarithm of a sequence of Counts, 
        !   the ith element of which represents the number of successes after SuccessStep(i) tries in a Bernoulli trail.
        !
        !   Input
        !
        !       maxNumTrial
        !
        !           The maximum number of trails possible. After maxNumTrial tries, 
        !           the Geometric distribution restarts from index 1.
        !
        !       numTrial
        !
        !           The length of the input vector SuccessStep and LogCount.
        !
        !       SuccessStep(1:numTrial)
        !
        !           An array of integers that represent the steps at which the Bernoulli successes occur. 
        !
        !           WARNING: THESE VALUES MUST BE INTEGER NUMBERS BETWEEN 1 AND MAXNUMTRIAL.
        !           WARNING: THE ONUS IS ON THE USER TO ENSURE THIS CONDITION HOLDS.
        !
        !       LogCount(1:numTrial)
        !
        !           A vector of real values representing the natural logarithms of the counts of success at the corresponding 
        !           Bernoulli trials specified by elements of SuccessStep. 
        !
        !   Output
        !
        !       PowellMinimum
        !
        !           An object of type PowellMinimum_type from Optimiziation_mod, containing the best-fit successProb and the 
        !           normalization constant of the fit in the vector component xmin.
        !
        use Optimization_mod, only: PowellMinimum_type
        use Constants_mod, only: IK, RK, POSINF_RK, NEGINF_RK
        implicit none
        integer(IK) , intent(in)    :: maxNumTrial
        integer(IK) , intent(in)    :: numTrial
        integer(IK) , intent(in)    :: SuccessStep(numTrial)
        real(RK)    , intent(in)    :: LogCount(numTrial)
        type(PowellMinimum_type)    :: PowellMinimum

        real(RK)                    :: BestFitSuccessProbNormFac(2) ! vector of the two parameters
        real(RK)    , parameter     :: SUCCESS_PROB_INIT_GUESS = 0.23_RK
        real(RK)    , parameter     :: FISHER_TRANS_SUCCESS_PROB_INIT_GUESS = atanh(2*(SUCCESS_PROB_INIT_GUESS - 0.5_RK))

        ! do Fisher transformation to make the limits infinity
        BestFitSuccessProbNormFac = [FISHER_TRANS_SUCCESS_PROB_INIT_GUESS, 0._RK] ! LogCount(1)]

        PowellMinimum = PowellMinimum_type  ( ndim = 2_IK &
                                            , getFuncMD = getSumDistSq &
                                            , StartVec = BestFitSuccessProbNormFac &
                                            )
        if (PowellMinimum%Err%occurred) return
        PowellMinimum%xmin(1) = 0.5_RK * tanh(PowellMinimum%xmin(1)) + 0.5_RK ! reverse Fisher-transform

    contains

        pure function getSumDistSq(ndim,successProbFisherTransNormFac) result(sumDistSq)
            !use Constants_mod, only: IK, RK
            implicit none
            integer(IK) , intent(in)    :: ndim
            real(RK)    , intent(in)    :: successProbFisherTransNormFac(ndim)
            real(RK)                    :: sumDistSq, successProb
            successProb = 0.5_RK * tanh(successProbFisherTransNormFac(1)) + 0.5_RK ! reverse Fisher-transform
            !sumDistSq = sum( (LogCount - getGeoLogPDF(successProb=successProb,seqLen=numTrial) - successProbFisherTransNormFac(2) )**2 )
            sumDistSq = sum(    ( LogCount &
                                - getLogProbGeoCyclic(successProb=successProb, maxNumTrial=maxNumTrial, numTrial=numTrial, SuccessStep=SuccessStep) &
                                - successProbFisherTransNormFac(2) &
                                )**2 &
                            )
        end function getSumDistSq

    end function fitGeoCyclicLogPDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns standard normal distribution PDF value.
    function getSNormPDF(z)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getSNormPDF
#endif
        use Constants_mod, only: INVSQRT2PI
        implicit none
        real(RK), intent(in) :: z
        real(RK)             :: getSNormPDF
        getSNormPDF = INVSQRT2PI * exp( -0.5_RK*z**2 )
    end function getSNormPDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns the non-standard normal distribution PDF value.
    function getNormPDF(mean,stdev,variance,x)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormPDF
#endif
        use Constants_mod, only: INVSQRT2PI
        implicit none
        real(RK), intent(in) :: mean,stdev,variance,x
        real(RK)             :: getNormPDF
        getNormPDF = INVSQRT2PI * exp( -(x-mean)**2/(2._RK*variance) ) / stdev
    end function getNormPDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getNormCDF(mean,stdev,x)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getNormCDF
#endif
        use Constants_mod, only: SPR,SQRT2
        implicit none
        real(RK), intent(in) :: mean,stdev,x
        real(RK)             :: getNormCDF
        getNormCDF = 0.5_RK * ( 1._RK + erf( real((x-mean)/(SQRT2*stdev),kind=SPR) ) )
    end function getNormCDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getSNormCDF(x)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getSNormCDF
#endif
        use Constants_mod, only: SPR,SQRT2
        implicit none
        real(RK), intent(in) :: x
        real(RK)             :: getSNormCDF
        getSNormCDF = 0.5_RK * ( 1._RK + erf( real(x/SQRT2,kind=SPR) ) )
    end function getSNormCDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ATTN: if x is not in [0,1], a negative value for getBetaCDF will be return to indicate error.
    ! The efficiency of this code can be improved by making x a vector on input.
    function getBetaCDF(a,b,x)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaCDF
#endif
        use Constants_mod, only : SPR
        implicit none
        real(RK), intent(in) :: a,b,x
        real(RK)             :: bt
        real(RK)             :: getBetaCDF
        if (x < 0._RK .or. x > 1._RK) then
            getBetaCDF = -1._RK
            return
        end if
        if (x == 0._RK .or. x == 1._RK) then
            bt = 0._RK
        else
            bt = exp( log_gamma(real(a+b,kind=SPR)) - log_gamma(real(a,kind=SPR)) - log_gamma(real(b,kind=SPR)) &
               + a*log(x) + b*log(1._RK-x) )
        end if
        if ( x < (a+1._RK) / (a+b+2._RK) ) then
            getBetaCDF = bt * getBetaContinuedFraction(a,b,x) / a
        else
            getBetaCDF = 1._RK - bt * getBetaContinuedFraction(b,a,1._RK-x) / b
        end if
    end function getBetaCDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getBetaContinuedFraction(a,b,x)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getBetaContinuedFraction
#endif
        implicit none
        real(RK)   , intent(in) :: a,b,x
        real(RK)   , parameter  :: eps = epsilon(x), fpmin = tiny(x)/eps
        integer(IK), parameter  :: maxit = 100
        real(RK)                :: aa,c,d,del,qab,qam,qap
        real(RK)                :: getBetaContinuedFraction
        integer(IK)             :: m,m2
        qab = a+b
        qap = a+1._RK
        qam = a-1._RK
        c = 1._RK
        d = 1._RK-qab*x/qap
        if (abs(d) < fpmin) d = fpmin
        d = 1._RK/d
        getBetaContinuedFraction = d
        do m = 1,maxit
            m2 = 2*m
            aa = m*(b-m)*x/((qam+m2)*(a+m2))
            d = 1._RK+aa*d
            if (abs(d) < fpmin) d = fpmin
            c = 1._RK+aa/c
            if (abs(c) < fpmin) c = fpmin
            d = 1._RK/d
            getBetaContinuedFraction = getBetaContinuedFraction*d*c
            aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
            d = 1._RK+aa*d
            if (abs(d) < fpmin) d = fpmin
            c = 1._RK+aa/c
            if (abs(c) < fpmin) c = fpmin
            d = 1._RK/d
            del = d*c
            getBetaContinuedFraction = getBetaContinuedFraction*del
            if (abs(del-1._RK) <=  eps) exit
        end do
        if (m > maxit) then
            error stop
            !call abortProgram( output_unit , 1 , 1 , &
            !'Statitistics@getBetaContinuedFraction() failed: a or b too big, or maxit too small.' )
        end if
    end function getBetaContinuedFraction

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine doKS1(np,Point,getCDF,statKS,probKS,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
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
            Err%msg = PROCEDURE_NAME//Err%msg
            return
        end if

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

    ! This assumes that points are coming from a uniform distribution in [0,1]. So, all Point must be in [0,1] on input.
    pure subroutine doUniformKS1(np,Point,statKS,probKS,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
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
            Err%msg = PROCEDURE_NAME//Err%msg
            return
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

    ! it is assumed that the input data1 and data2 arrays are sorted in ascending-order on input
    subroutine doSortedKS2(np1,np2,SortedPoint1,SortedPoint2,statKS,probKS)
#if defined DLL_ENABLED && !defined CFI_ENABLED
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

    pure function getProbKS(lambda) result(probKS)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getProbKS
#endif
        implicit none
        real(RK)   , intent(in) :: lambda
        real(RK)   , parameter  :: EPS1 = 0.001_RK, EPS2 = 1.e-8_RK
        integer(IK), parameter  :: NITER=100
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

    ! returns the uniform CDF on support [0,1). rather redundant, aint it? but sometimes, needed
    pure function getUniformCDF(x)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getUniformCDF
#endif
        implicit none
        real(RK), intent(in) :: x
        real(RK)             :: getUniformCDF
        getUniformCDF = x
    end function getUniformCDF

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Returns the 1-D histogram (Density plot) of a vector X.
    ! The number of bins in the X range (nxbin) is determined by the user
    ! The range of X (xmin,xmax) should also be given by the user
    ! The program returns 2 arrays of Xbin & Density(x) as output.
    ! Amir Shahmoradi, Sep 1, 2017, 12:30 AM, ICES, UT Austin
    subroutine getHist1d(method,xmin,xmax,nxbin,np,X,Xbin,Density,errorOccurred)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getHist1d
#endif
        implicit none
        character(*), intent(in)  :: method
        integer(IK) , intent(in)  :: np,nxbin
        real(RK)    , intent(in)  :: xmin,xmax
        real(RK)    , intent(in)  :: X(np)
        real(RK)    , intent(out) :: Xbin(nxbin),Density(nxbin)
        logical     , intent(out) :: errorOccurred
        real(RK)                  :: xbinsize
        integer(IK)               :: i,ip,thisXbin

        errorOccurred = .false.
        Density = 0._RK
        xbinsize = (xmax-xmin) / real(nxbin,kind=RK)
        Xbin = [ (xmin+real(i-1,kind=RK)*xbinsize,i=1,nxbin) ]

        do ip = 1,np
            thisXbin = getBin(X(ip),xmin,nxbin,xbinsize)
            Density(thisXbin) = Density(thisXbin) + 1._RK
        end do

        Xbin = Xbin + 0.5_RK * xbinsize
       !method = getLowerCase(trim(adjustl(histType)))
        if(method=='pdf') then
            Density = Density / real(np,kind=RK)
        elseif(method=='count') then
            return
        else
            errorOccurred = .true.
        end if
    end subroutine getHist1d

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Returns the 2-D histogram (Density plot) of a set of data points with (X,Y) coordinates.
    ! The number of bins in the X and Y directions (nxbin, nybin) are determined by the user
    ! The range of X & Y (xmin,xmax,ymin,ymax) should also be given by the user
    ! The program returns 3 arrays of Xbin, Ybin & Density(y,x) as the output.
    ! Amir Shahmoradi, Sep 1, 2017, 12:00 AM, ICES, UT Austin
    subroutine getHist2d(histType,xmin,xmax,ymin,ymax,nxbin,nybin,np,X,Y,Xbin,Ybin,Density)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getHist2d
#endif
        !use, intrinsic :: iso_fortran_env, only: output_unit
        use String_mod, only: getLowerCase
        implicit none
        character(*), intent(in)  :: histType
        integer(IK) , intent(in)  :: np,nxbin,nybin
        real(RK)    , intent(in)  :: xmin,xmax,ymin,ymax
        real(RK)    , intent(in)  :: X(np),Y(np)
        real(RK)    , intent(out) :: Xbin(nxbin),Ybin(nybin),Density(nybin,nxbin)
        character(:), allocatable :: method
        real(RK)                  :: xbinsize,ybinsize
        integer(IK)               :: i,ip,thisXbin,thisYbin

        Density = 0._RK
        xbinsize = (xmax-xmin) / real(nxbin,kind=RK)
        ybinsize = (ymax-ymin) / real(nybin,kind=RK)
        Xbin = [ (xmin+real(i-1,kind=RK)*xbinsize,i=1,nxbin) ]
        Ybin = [ (ymin+real(i-1,kind=RK)*ybinsize,i=1,nybin) ]

        do ip = 1,np
            thisXbin = getBin(X(ip),xmin,nxbin,xbinsize)
            thisYbin = getBin(Y(ip),ymin,nybin,ybinsize)
            Density(thisYbin,thisXbin) = Density(thisYbin,thisXbin) + 1._RK
        end do

        Xbin = Xbin + 0.5_RK * xbinsize
        Ybin = Ybin + 0.5_RK * ybinsize
        method = getLowerCase(trim(adjustl(histType)))
        if(method=='pdf') then
            Density = Density / real(np,kind=RK)
        elseif(method=='pdf(y|x)') then
            do i = 1,nxbin
                Density(1:nybin,i) = Density(1:nybin,i) / sum(Density(1:nybin,i))
            end do
        elseif(method=='pdf(x|y)') then
            do i = 1,nybin
                Density(i,1:nxbin) = Density(i,1:nxbin) / sum(Density(i,1:nxbin))
            end do
        elseif(method=='count') then
            return
        else
            error stop
            !call abortProgram( output_unit , 1 , 1 , 'Statistics@getHist2d() failed. The requested histType does not exist.' )
        end if

    end subroutine getHist2d

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ATTN: if x <= xmin or x xmin+nbin*binsize then getBin=-1 will be returned to indicate error.
    ! Given the range of the variable x, (xmin:xmin+binsize*nbin), and the number of bins, nbin,
    ! which this range is going to be divided, this FUNCTION finds which bin\in[1,nbin], the input variable X belongs to.
    ! The output is the number that identifies the bin.
    ! On input, binsize must be exactly (xmax-xmin)/nbin.
    ! if bmin < x <= bmax then x belongs to this bin.
    ! Version 3.0, Sep 1, 2017, 11:12 AM, Amir Shahmoradi, ICES, The University of Texas at Austin.
    pure function getBin(x,lowerBound,nbin,binsize)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getBin
#endif

        implicit none
        integer(IK), intent(in) :: nbin
        real(RK)   , intent(in) :: x,lowerBound,binsize
        real(RK)                :: xmin,xmid
        integer(IK)             :: getBin,minbin,midbin,maxbin

        if (x<lowerBound .or. x>=lowerBound+nbin*binsize) then
            getBin = -1_IK
            return
        end if

        minbin = 1
        maxbin = nbin
        xmin = lowerBound
        loopFindBin: do
            midbin = (minbin+maxbin) / 2
            xmid = xmin + midbin*binsize
            if (x<xmid) then
                if (minbin==midbin) then
                    getBin = minbin
                    exit loopFindBin
                end if
                maxbin = midbin
                cycle loopFindBin
            else
                if (minbin==midbin) then
                    getBin = maxbin
                    exit loopFindBin
                end if
                minbin = midbin
                cycle loopFindBin
            end if
        end do loopFindBin

    end function getBin

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getQuantile(np,nq,SortedQuantileProbability,Point,Weight,sumWeight) result(Quantile)
#if defined DLL_ENABLED && !defined CFI_ENABLED
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
        real(RK)                            :: probability
        integer(IK)                         :: ip, iq, iw, weightCounter, Indx(np), SortedQuantileDensity(nq)
        type(Err_type)                      :: Err
        iq = 1_IK
        Quantile = 0._RK
        probability = 0._RK
        weightCounter = 0_IK
        call indexArray(np,Point,Indx,Err)
        if (Err%occurred) then
            Quantile = NEGINF_RK
            return
        end if
        if (present(sumWeight)) then
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
    end function getQuantile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Statistics_mod