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

!>  \brief This module contains mathematical procedures.
!>  \author Amir Shahmoradi

module Knn_mod

    use Err_mod, only: Err_type
    use Constants_mod, only: IK, RK
    use JaggedArray_mod, only: IV => IntVec_type
    use JaggedArray_mod, only: RV => RealVec_type
    implicit none

    character(*), parameter :: MODULE_NAME = "@Knn_mod"

    type :: MinDistEdge_type
        integer(IK)                 :: np
        integer(IK) , allocatable   :: Index(:) ! np
        real(RK)    , allocatable   :: LenSq(:) ! np
    end type MinDistEdge_type

    interface MinDistEdge_type
        module procedure :: getMinDistEdge
    end interface MinDistEdge_type

    type :: Hub_type
        integer(IK)                 :: nh
        integer(IK) , allocatable   :: NodeIndex(:) ! nh
        integer(IK) , allocatable   :: EdgeCount(:) ! nh
        type(IV)    , allocatable   :: EdgeIndex(:) ! nh
        type(RV)    , allocatable   :: EdgeLenSq(:) ! nh
        type(MinDistEdge_type)      :: MinDistEdge
        type(Err_type)              :: Err
    end type Hub_type

    interface Hub_type
        module procedure :: hubify
    end interface Hub_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the probability density function of the `k`th nearest neighbor distance being `distSq` in a Poisson process of
    !> dimension `nd` with the rate density `lambda`.
    !>
    pure function getLogProbKnnPoisson(nd, k, lambda, distance) result(logProbKnnPoisson)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogProbKnnPoisson
#endif
        use Constants_mod, only: IK, RK
        use Math_mod, only: getLogVolUnitBall
        use Math_mod, only: getLogFactorial
        implicit none
        integer(IK) , intent(in)    :: nd, k
        real(RK)    , intent(in)    :: lambda, distance
        real(RK)                    :: logProbKnnPoisson
        real(RK)                    :: logVolUnitBall
        real(RK)                    :: logDistance
        logDistance = log(distance)
        logVolUnitBall = getLogVolUnitBall(nd)
        logProbKnnPoisson   = log(real(nd,RK)) + k * (log(lambda) + logVolUnitBall + nd * logDistance) &
                            - lambda * exp(logVolUnitBall) * distance**nd &
                            - logDistance - getLogFactorial(k-1_IK)
    end function getLogProbKnnPoisson

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the average minimum distance in dimension `nd` with the rate density `logDensity`.
    pure function getMeanMinDist(nd, logDensity, logVolUnitBall) result(logMeanMinDist)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMeanMinDist
#endif
        use Constants_mod, only: IK, RK ! LCOV_EXCL_LINE
        !use Math_mod, only: getLogVolUnitBall
        implicit none
        integer(IK) , intent(in)    :: nd
        real(RK)    , intent(in)    :: logDensity, logVolUnitBall
        real(RK)                    :: logMeanMinDist
        real(RK)                    :: ndInverse
        !logVolUnitBall = getLogVolUnitBall(nd)
        ndInverse = 1._RK / nd
        logMeanMinDist = log_gamma((nd + 1_IK) * ndInverse) - ndInverse * (logDensity + logVolUnitBall)
    end function getMeanMinDist

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the log rate density `logDensity` form the average minimum distance in dimension `nd`.
    pure function getLogDensity(nd, logMeanMinDist, logVolUnitBall) result(logDensity)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogDensity
#endif
        use Constants_mod, only: IK, RK ! LCOV_EXCL_LINE
        implicit none
        integer(IK) , intent(in)    :: nd
        real(RK)    , intent(in)    :: logMeanMinDist
        real(RK)    , intent(in)    :: logVolUnitBall
        real(RK)                    :: logDensity
        logDensity = nd * (log_gamma((nd + 1_IK) / real(nd,RK)) - logMeanMinDist) - logVolUnitBall
    end function getLogDensity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getShortestSumDistSqPointIndex(nd, np, Point) result(innestIndex)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getShortestSumDistSqPointIndex
#endif
        use Constants_mod, only: IK, RK, HUGE_RK ! LCOV_EXCL_LINE
        implicit none
        integer(IK) , intent(in)    :: nd, np
        real(RK)    , intent(in)    :: Point(nd,np)
       !real(RK)                    :: DistanceSq(np,np)
        integer(IK)                 :: innestIndex
        integer(IK)                 :: ip, jp
        real(RK)                    :: sumDistSq
        real(RK)                    :: sumDistSqOld
        innestIndex = 0_IK
        sumDistSqOld = HUGE_RK
        do jp = 1, np
            sumDistSq = 0._RK
            do ip = 1, np
                !DistanceSq(ip,jp) = sum((Point(1:nd,ip) - Point(1:nd,jp))**2)
                sumDistSq = sumDistSq + sum((Point(1:nd,jp) - Point(1:nd,ip))**2)
            end do
            if (sumDistSq < sumDistSqOld) then
                sumDistSqOld = sumDistSq
                innestIndex = jp
            end if
        end do
    end function getShortestSumDistSqPointIndex

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getDistSq(nd, np, Point) result(DistSq)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDistSq
#endif
        use Constants_mod, only: IK, RK ! LCOV_EXCL_LINE
        implicit none
        integer(IK) , intent(in)    :: nd, np
        real(RK)    , intent(in)    :: Point(nd,np)
        real(RK)                    :: DistSq(np,np)
        integer(IK)                 :: jp, ip
        do concurrent(jp = 1:np)
            DistSq(1:jp-1,jp) = DistSq(jp,1:jp-1)
            DistSq(jp,jp) = 0._RK
            do ip = jp + 1, np
                DistSq(ip,jp) = sum( (Point(1:nd,ip) - Point(1:nd,jp))**2 )
            end do
        end do
    end function getDistSq

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \warning
    !> The input value `np` must be larger than 1 at all times. A value of `np = 1` is meaningless.
    !> In addition, a value of `np = 1` is also meaningless, although this procedure can handle `np = 2` gracefully.
    pure function getMinDistEdge(np, DistSq) result(MinDistEdge)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getMinDistEdge
#endif
        use Constants_mod, only: IK, RK, HUGE_RK ! LCOV_EXCL_LINE
        implicit none
        integer(IK) , intent(in)    :: np
        real(RK)    , intent(in)    :: DistSq(np,np)
        type(MinDistEdge_type)      :: MinDistEdge
        integer(IK)                 :: pindex1
        integer(IK)                 :: pindex2
        integer(IK)                 :: jp

        allocate(MinDistEdge%Index(np), MinDistEdge%LenSq(np))

        MinDistEdge%np = np

        jp = 1_IK
        pindex2 = minloc(DistSq(2:np,jp), dim = 1) + 1
        MinDistEdge%Index(jp) = pindex2
        MinDistEdge%LenSq(jp) = DistSq(pindex2,jp)

        jp = np
        pindex1 = minloc(DistSq(1:jp-1,jp), dim = 1)
        MinDistEdge%Index(jp) = pindex1
        MinDistEdge%LenSq(jp) = DistSq(pindex1,jp)

        do concurrent(jp = 2:np-1)
            pindex1 = minloc(DistSq(1:jp-1,jp), dim = 1)
            pindex2 = minloc(DistSq(jp+1:np,jp), dim = 1) + jp
            if (DistSq(pindex1,jp) < DistSq(pindex2,jp)) then
                MinDistEdge%Index(jp) = pindex1
                MinDistEdge%LenSq(jp) = DistSq(pindex1,jp)
            else
                MinDistEdge%Index(jp) = pindex2
                MinDistEdge%LenSq(jp) = DistSq(pindex2,jp)
            end if
        end do

    end function getMinDistEdge

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function hubify(np, DistSq) result(Hub)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: hubify
#endif
        use Constants_mod, only: IK, RK, HUGE_RK
        use Unique_mod, only: findUnique

        implicit none

        integer(IK) , intent(in)    :: np
        real(RK)    , intent(in)    :: DistSq(np,np)
        type(Hub_type)              :: Hub

        character(*), parameter     :: PROCEDURE_NAME = MODULE_NAME//"@hubify()"

        integer(IK)                 :: ih

        Hub%MinDistEdge = MinDistEdge_type(np, DistSq)
        call findUnique ( lenVector = np & ! LCOV_EXCL_LINE
                        , Vector = Hub%MinDistEdge%Index & ! LCOV_EXCL_LINE
                        , lenUnique = Hub%nh & ! LCOV_EXCL_LINE
                        , UniqueValue = Hub%NodeIndex & ! LCOV_EXCL_LINE
                        , UniqueCount = Hub%EdgeCount & ! LCOV_EXCL_LINE
                        , UniqueIndex = Hub%EdgeIndex & ! LCOV_EXCL_LINE
                        , Err = Hub%Err & ! LCOV_EXCL_LINE
                        )
        if (Hub%Err%occurred) then
            ! LCOV_EXCL_START
            Hub%Err%msg = PROCEDURE_NAME//Hub%Err%msg
            return
            ! LCOV_EXCL_STOP
        end if

        allocate(Hub%EdgeLenSq(Hub%nh))
        do concurrent(ih = 1:Hub%nh)
            Hub%EdgeLenSq(ih)%Vector = Hub%MinDistEdge%LenSq(Hub%EdgeIndex(ih)%Vector)
        end do

    end function hubify

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Knn_mod ! LCOV_EXCL_LINE