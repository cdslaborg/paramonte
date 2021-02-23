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

!>  \brief This module contains mathematical/computational geometry procedures.
!>  \author Amir Shahmoradi

module Geometry_mod

    implicit none

    character(*), parameter :: MODULE_NAME = "@Geometry_mod"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Compute the *effective* total volume of all bounded ellipsoidal regions together while considering the possible overlaps.
    !> This computation is done aggressively via Monte Carlo simulation.
    !>
    !> \param[in]       nd              :   The number of dimensions.
    !> \param[in]       nc              :   The number of input clusters.
    !> \param[in]       nc              :   The number of input clusters.
    !> \param[in]       nsim            :   The number of simulated points for computing the effective sum of the volumes of all clusters `logSumVolNormed`
    !>                                      (**optional**, if not provided, the sum of the cluster bounding volumes will be computed as if no bounding ellipsoids overlap).
    !> \param[in]       Center          :   An array of size `(nd,nc)` representing the centers of the ellipsoids.
    !> \param[in]       InvCovMat       :   An array of size `(nd,nd,nc)` representing the inverse covariance matrices of the ellipsoids.
    !> \param[in]       LogVolNormed    :   An array of size `(nc)` representing the logarithms of the volumes of the ellipsoids normalized to the volume of an nd-ball.
    !> \param[in]       ScaleFactorSq   :   An array of size `(nc)` representing the maximum of Mahalanobis-distance-squared of the members of each cluster from its center.
    !> \param[in]       ChoLowCovUpp    :   An array of size `(nd,nd,nc)` whose lower-triangle represents the Cholesky lower triangle of the covariance matrices of the ellipsoids.
    !> \param[in]       ChoDia          :   An array of size `(nd,nc)` representing the diagonal elements of the Cholesky lower triangle of the covariance matrices of the ellipsoids.
    !> \param[out]      Overlap         :   An array of size `(nc,nc)` representing the diagonal elements of the Cholesky lower triangle of the covariance matrices of the ellipsoids.
    !> \param[inout]    logSumVolNormed :   The logarithm of the sum of `LogVolNormed`. On output, it will be rewritten with the effective log-sum of the (normalized) volumes.
    subroutine getEffLogVol ( nd & ! LCOV_EXCL_LINE
                            , nc & ! LCOV_EXCL_LINE
                            , nsim & ! LCOV_EXCL_LINE
                            , Center & ! LCOV_EXCL_LINE
                            , InvCovMat & ! LCOV_EXCL_LINE
                            , LogVolNormed & ! LCOV_EXCL_LINE
                            , ScaleFactorSq & ! LCOV_EXCL_LINE
                            , ChoLowCovUpp & ! LCOV_EXCL_LINE
                            , ChoDia & ! LCOV_EXCL_LINE
                            , logSumVolNormed & ! LCOV_EXCL_LINE
                            , Overlap & ! LCOV_EXCL_LINE
                            )
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getEffLogVol
#endif
        use Statistics_mod, only: getRandMVU
        use Constants_mod, only: IK, RK
        implicit none

        integer(IK) , intent(in)    :: nd, nc, nsim
        real(RK)    , intent(in)    :: Center(nd,nc)
        real(RK)    , intent(in)    :: ChoDia(nd,nc)
        real(RK)    , intent(in)    :: LogVolNormed(nc)
        real(RK)    , intent(in)    :: ScaleFactorSq(nc)
        real(RK)    , intent(in)    :: InvCovMat(nd,nd,nc)
        real(RK)    , intent(in)    :: ChoLowCovUpp(nd,nd,nc)
        real(RK)    , intent(inout) :: logSumVolNormed
        integer(IK) , intent(out)   :: Overlap(nc,nc)

        integer(IK)                 :: ic, jc, isim
        integer(IK)                 :: membershipCount
        integer(IK) , allocatable   :: NumSim(:)
        real(RK)    , allocatable   :: VolNormed(:)
        real(RK)    , allocatable   :: SimPointVol(:)
        real(RK)                    :: SimPoint(nd)
        real(RK)                    :: SimPointNormed(nd)
        real(RK)                    :: volumeContribution
        real(RK)                    :: sumVolNormed
        real(RK)                    :: mahalSq

        ! Compute the expected number of simulated points in each bounded region.

        VolNormed = exp(LogVolNormed - logSumVolNormed)
        NumSim = ceiling(nsim * VolNormed)
        SimPointVol = VolNormed / NumSim

        ! Simulate points within each cluster and compute the effective sum of volumes by excluding overlaps.

        sumVolNormed = 0._RK
        do ic = 1, nc

            do isim = 1, NumSim(ic)

                SimPoint(1:nd) = getRandMVU ( nd = nd & ! LCOV_EXCL_LINE
                                            , MeanVec = Center(1:nd,ic) & ! LCOV_EXCL_LINE
                                            , CholeskyLower = ChoLowCovUpp(1:nd,1:nd,ic) & ! LCOV_EXCL_LINE
                                            , Diagonal = ChoDia(1:nd,ic) & ! LCOV_EXCL_LINE
                                            )

                membershipCount = 1_IK
                volumeContribution = SimPointVol(ic)
                do jc = 1, nc
                    if (jc/=ic) then
                        SimPointNormed(1:nd) = SimPoint(1:nd) - Center(1:nd,jc)
                        mahalSq = dot_product( SimPointNormed(1:nd) , matmul(InvCovMat(1:nd,1:nd,jc), SimPointNormed(1:nd)) )
                        if (mahalSq < ScaleFactorSq(jc)) then
                            membershipCount = membershipCount + 1_IK
                            volumeContribution = volumeContribution + SimPointVol(jc)
                        end if
                    end if
                end do

                if (membershipCount==1_IK) then
                    sumVolNormed = sumVolNormed + volumeContribution
                else
                    sumVolNormed = sumVolNormed + volumeContribution / membershipCount
                end if

            end do

        end do

        logSumVolNormed = log(sumVolNormed) + logSumVolNormed

    end subroutine getEffLogVol

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Geometry_mod ! LCOV_EXCL_LINE