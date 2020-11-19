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

!>  \brief This module contains procedures for computing the cross-correlation of time series data.
!>  @author Amir Shahmoradi

module CrossCorr_mod

    implicit none

    character(len=*), parameter :: MODULE_NAME = "@CrossCorr_mod"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the integrated autocorrelation (IAC) via the BatchMeans method.
    !>
    !> @param[in]   np              : The number of data points in the input time series data.
    !> @param[in]   Point           : The input data series data vector.
    !> @param[in]   Weight          : The vector of weights of the input data points (optional, default = array of ones).
    !> @param[in]   batchSize       : The batch size (optional, default = computed from the input parameters).
    !>
    !> \return
    !> `iac` : The integrated autocorrelation (IAC) via the BatchMeans method.
    !>
    !> \remark
    !> Note that np must be large enough to get a meaningful answer.
    function getBatchMeansIAC(np,Point,Weight,batchSize) result(iac)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getBatchMeansIAC
#endif
        use Constants_mod, only: IK, RK
        use Statistics_mod, only: getVariance
        use Math_mod, only: getCumSum
        !use Stirng_mod, only: num2str
        !use Err_mod, only: Err_type

        implicit none

        character(len=*), parameter         :: PROCEDURE_NAME = MODULE_NAME // "@getBatchMeansIAC()"

        integer(IK), intent(in)             :: np
        real(RK), intent(in)                :: Point(np)
        integer(IK), intent(in), optional   :: Weight(np), batchSize
        !type(Err_type), intent(out)         :: Err
        real(RK)                            :: iac

        integer(IK)                         :: CumSumWeight(np), currentSampleEndLoc, sampleSize, sampleCount, isample
        integer(IK)                         :: ip, ipVerbose, ipStart, ipEnd, npEffective
        real(RK)                            :: avgPoint, varPoint, avgBatchMean, varBatchMean
        real(RK)                            :: diffSquared, sampleSizeInverse
        real(RK), allocatable               :: BatchMean(:)

        !Err%occurred = .false.

        if (present(Weight)) then
            CumSumWeight = getCumSum(np,Weight)
        else
            CumSumWeight(np) = np
        end if

        ! compute batch size and count
        if (present(batchSize)) then
            sampleSize = batchSize
        else
            sampleSize = int( real(CumSumWeight(np),kind=RK)**(0.666666666666666_RK) )
        end if
        sampleSizeInverse = 1._RK / real(sampleSize,kind=RK)
        sampleCount = CumSumWeight(np) / sampleSize
        npEffective = sampleCount * sampleSize

        if (sampleCount<2) then
            !Err%occurred = .true.
            !Err%msg = PROCEDURE_NAME // ": sampleCount<10: " // num2str(sampleCount)
            !iac = -huge(iac)
            iac = 1
            return
        end if

        ! xxx: here goes another GFortran 7.3 bug: EndOfLineLoc is assumed already allocated, despite the first appearance here.
        if (allocated(BatchMean)) deallocate(BatchMean)
        allocate(BatchMean(sampleCount))

        ! improvement: iterate from the end to the beginning of the chain to ignore initial points instead of the last points.
        ! this would be beneficial for MCMC samples

        ! compute the Batch-Means avergage and variance, also average of Point
        avgPoint = 0._RK
        ipStart = 1
        ipEnd   = 0
        if (present(Weight)) then
            ip = 1
            isample = 1
            ipVerbose = 0
            currentSampleEndLoc = sampleSize
            BatchMean(isample) = 0._RK
            loopOverWeight: do
                ipVerbose = ipVerbose + 1
                if (ipVerbose>CumSumWeight(ip)) ip = ip + 1
                if (ipVerbose>currentSampleEndLoc) then ! we are done with the current batch
                    avgPoint = avgPoint + BatchMean(isample)
                    BatchMean(isample) = BatchMean(isample) * sampleSizeInverse
                    if (ipVerbose>npEffective) exit loopOverWeight  ! condition equivalent to currentSampleEndLoc==npEffective
                    currentSampleEndLoc = currentSampleEndLoc + sampleSize
                    isample = isample + 1
                    BatchMean(isample) = 0._RK
                end if
                BatchMean(isample)  = BatchMean(isample) + Point(ip)
            end do loopOverWeight
        else    ! there is no weight
            do isample = 1, sampleCount
                BatchMean(isample) = 0._RK
                ipEnd = ipEnd + sampleSize
                do ip = ipStart, ipEnd
                    BatchMean(isample)  = BatchMean(isample) + Point(ip)
                end do
                ipStart = ipEnd + 1_IK
                avgPoint = avgPoint + BatchMean(isample)
                BatchMean(isample) = BatchMean(isample) * sampleSizeInverse
            end do
        end if
        avgBatchMean = sum( BatchMean ) / real(sampleCount,kind=RK)
        varBatchMean = sum( (BatchMean - avgBatchMean)**2 ) / real(sampleCount-1,kind=RK)
        avgPoint = avgPoint / real(npEffective,kind=RK)

        ! compute the variance of Point

        varPoint = 0._RK
        if (present(Weight)) then
            ip = 1
            ipVerbose = 0
            diffSquared = ( Point(ip) - avgPoint )**2
            loopComputeVarPoint: do
                ipVerbose = ipVerbose + 1
                if (ipVerbose>npEffective) exit loopComputeVarPoint
                if (ipVerbose>CumSumWeight(ip)) then
                    ip = ip + 1 ! by definition, ip never become > np, otherwise it leads to disastrous errors
                    diffSquared = ( Point(ip) - avgPoint )**2
                end if
                varPoint = varPoint + diffSquared
            end do loopComputeVarPoint
        else
            do ip = 1, npEffective
                varPoint = varPoint + ( Point(ip) - avgPoint )**2
            end do
        end if
        varPoint = varPoint / real(npEffective-1,kind=RK)

        ! compute the IAC

        iac = sampleSize * varBatchMean / varPoint

    end function getBatchMeansIAC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the integrated autocorrelation (IAC) based on the cumulative autocorrelation.
    !>
    !> @param[in]   np              :   The number of data points in the input time series data.
    !> @param[in]   Point           :   The input data series data vector.
    !> @param[in]   Weight          :   The vector of weights of the input data points (optional, default = array of ones).
    !> @param[in]   significance    :   The significance in units of standard deviation below which the autocorrelation is 
    !>                                  considered noise (optional, default = 2).
    !>
    !> \return
    !> `iac` : The integrated autocorrelation (IAC).
    function getCumSumIAC(np,Point,Weight,significance) result(iac)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSumIAC
#endif
        use Constants_mod, only: IK, RK
        use Math_mod, only: getCumSum
        implicit none
        integer(IK) , intent(in)            :: np
        real(RK)    , intent(in)            :: Point(np)
        integer(IK) , intent(in), optional  :: Weight(np), significance ! in units of sigma
        real(RK)                            :: iac, meanPoint, normFac, NormedData(np), cutoff
        real(RK)    , allocatable           :: AutoCorr(:)
        integer(IK)                         :: i, paddedLen, sumWeight, significanceDefault, cutoffIndex
        significanceDefault = 2_IK
        if (present(significance)) significanceDefault = significance
        if (present(Weight)) then
            sumWeight = sum(Weight)
            meanPoint = sum(Point*Weight) / real(sumWeight,kind=RK)
        else
            sumWeight = np
            meanPoint = sum(Point) / real(np,kind=RK)
        end if
        NormedData = Point - meanPoint
        paddedLen = getPaddedLen(sumWeight)
        AutoCorr = getCrossCorrFFTweighted  ( nData1    = np            &
                                            , nData2    = np            &
                                            , paddedLen = paddedLen     &
                                            , Data1     = NormedData    &
                                            , Data2     = NormedData    &
                                            , Weight1   = Weight        &
                                            , Weight2   = Weight        &
                                            )
        normFac = 1._RK / AutoCorr(1)
        AutoCorr = AutoCorr * normFac
        ! For autocorrelation, under the assumption of a completely random series, the ACF standard error reduces to sqrt(1/ndata)
        cutoff = significanceDefault * sqrt(1._RK/sumWeight) ! standardErrorAutoCorr
        cutoffIndex = 1_IK
        do i = 1, paddedLen
            if (AutoCorr(i)<cutoff) then
                cutoffIndex = i
                exit
            end if
        end do
        iac = 2_IK * sum(AutoCorr(1:cutoffIndex)) - 1_IK
    end function getCumSumIAC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the integrated autocorrelation (IAC) based on the maximum cumulative autocorrelation.
    !>
    !> @param[in]   np              : The number of data points in the input time series data.
    !> @param[in]   Point           : The input data series data vector.
    !> @param[in]   Weight          : The vector of weights of the input data points (optional, default = array of ones).
    !>
    !> \return
    !> `maxIAC` : The integrated autocorrelation (IAC).
    function getMaxCumSumIAC(np,Point,Weight) result(maxIAC)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getMaxCumSumIAC
#endif
        use Constants_mod, only: IK, RK
        use Math_mod, only: getCumSum
        implicit none
        integer(IK), intent(in)             :: np
        real(RK), intent(in)                :: Point(np)
        integer(IK), intent(in), optional   :: Weight(np)
        real(RK)                            :: maxIAC, meanPoint, normFac, NormedData(np)
        real(RK), allocatable               :: AutoCorr(:)
        integer(IK)                         :: paddedLen, sumWeight
        if (present(Weight)) then
            sumWeight = sum(Weight)
            meanPoint = sum(Point*Weight) / real(sumWeight,kind=RK)
        else
            sumWeight = np
            meanPoint = sum(Point) / real(np,kind=RK)
        end if
        NormedData = Point - meanPoint
        paddedLen = getPaddedLen(sumWeight)
        AutoCorr = getCrossCorrFFTweighted  ( nData1    = np            &
                                            , nData2    = np            &
                                            , paddedLen = paddedLen     &
                                            , Data1     = NormedData    &
                                            , Data2     = NormedData    &
                                            , Weight1   = Weight        &
                                            , Weight2   = Weight        &
                                            )
        normFac = 1._RK / AutoCorr(1)
        AutoCorr = AutoCorr * normFac
        maxIAC = 2_IK * maxval(getCumSum(paddedLen,AutoCorr)) - 1_IK
    end function getMaxCumSumIAC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the smallest length of vector that is a power of `base` and larger than the input vector length `actualLen`.
    !>
    !> @param[in]   actualLen   : The input vector length.
    !> @param[in]   base        : The base of the exponentiation.
    !>
    !> \return
    !> `paddedLen` : The minimum power-of-`base` length given `actualLen`.
    !>
    !> \remark
    !> This method is used to compute the cross-correlation. Usage:
    !> `actualLen = max( actualLen1, actualLen2 )`
    !>
    !> \remark
    !> For weighted-data cross-correlation computation, try,
    !> `actualLen = max( sum(Weight1(1:actualLen1)), sum(Weight2(1:actualLen2)) )`.
    pure function getPaddedLen(actualLen,base) result(paddedLen)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getPaddedLen
#endif
        use Constants_mod, only: IK, RK
        integer(IK), intent(in)         :: actualLen
        real(RK), intent(in), optional  :: base
        integer(IK)                     :: paddedLen
        paddedLen = 2 ** ( getNextExponent(real(actualLen,kind=RK),base) + 1 )
    end function getPaddedLen

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> Return the exponent that yields the smallest real number larger than the input number `absoluteValue`.
    !>
    !> @param[in]   absoluteValue   : The input real number.
    !> @param[in]   base            : The base of the exponentiation.
    !>
    !> \return
    !> `nextExponent` : The output minimum exponent.
    pure function getNextExponent(absoluteValue,base) result(nextExponent)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getNextExponent
#endif
        use Constants_mod, only: IK, RK, INVLN2
        real(RK), intent(in)            :: absoluteValue
        real(RK), intent(in), optional  :: base
        integer(IK)                     :: nextExponent
        if (present(base)) then
            nextExponent = ceiling( log(absoluteValue) / log(base) )
        else    ! assume the base is 2
            nextExponent = ceiling( log(absoluteValue) * INVLN2 )
        end if
    end function getNextExponent

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return an array that is extended and padded with zeros for the requested length `paddedLen`.
    !>
    !> @param[in]   currentLen  : The length of the input array.
    !> @param[in]   Array       : The array to be extended.
    !> @param[in]   paddedLen   : The requested new length of the array.
    !>
    !> \return
    !> `ArrayPadded` : The output extended array, padded with zeros.
    pure function padZero(currentLen,Array,paddedLen) result(ArrayPadded)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: padZero
#endif
        use Constants_mod, only: IK, RK
        integer(IK) , intent(in)            :: currentLen
        real(RK)    , intent(in)            :: Array(currentLen)
        integer(IK) , intent(in), optional  :: paddedLen
        real(RK)    , allocatable           :: ArrayPadded(:)
        integer(IK)                         :: i, paddedSize
        if (present(paddedLen)) then
            paddedSize = paddedLen
        else
            paddedSize = 2 ** ( getNextExponent(real(currentLen,kind=RK)) + 1 )
        end if
        allocate(ArrayPadded(paddedSize))
        do concurrent(i=1:currentLen)
            ArrayPadded(i) = Array(i)
        end do
        do concurrent(i=currentLen+1:paddedSize)
            ArrayPadded(i) = 0._RK
        end do
    end function padZero

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the cross-correlation of the two input data vectors, (including any user-supplied zero padding), computed via Fast-Fourier Transform.
    !>
    !> @param[in]   ndata   : The lengths of the input arrays. It MUST be an integer power of two.
    !> @param[in]   Data1   : The first array.
    !> @param[in]   Data2   : The second array.
    !>
    !> \return
    !> `CrossCorrFFT` : A vector of same length as the input arrays containing the cross-correlation.
    !> The answer is returned as the first `ndata` points in ans stored in wrap-around order, i.e., cross-correlation
    !> at increasingly negative lags are in `CrossCorrFFT(ndata)` on down to `CrossCorrFFT(ndata/2+1)`, while cross-correlation
    !> increasingly positive lags are in `CrossCorrFFT(1)` (zero lag) on up to `CrossCorrFFT(ndata/2)`.
    !>
    !> \remark
    !> Sign convention of this routine: If `Data1` lags Data2, i.e., is shifted to the
    !> right of it, then ans will show a peak at positive lags.
    !>
    !> \remark
    !> For autocorrelation, under the assumption of a completely random series,
    !> the ACF standard error reduces to: \f$\sqrt{ 1 / \texttt{ndata} }\f$
    function getCrossCorrFFT(ndata,Data1,Data2) result(CrossCorrFFT)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getCrossCorrFFT
#endif
        use Constants_mod, only: IK, RK, CK
        implicit none
        integer(IK), intent(in)                 :: ndata
        real(RK), intent(inout)                 :: Data1(ndata),Data2(ndata)
        real(RK)                                :: CrossCorrFFT(ndata)
        character(len=*), parameter             :: PROCEDURE_NAME = MODULE_NAME//"@getCrossCorrFFT()"
        complex(CK), dimension(ndata/2)         :: Cdat1,Cdat2
        integer(IK)                             :: no2
        if (iand(ndata,ndata-1)/=0) then
            write(*,*) PROCEDURE_NAME//": paddedLen must be a power of 2."
            error stop
        end if
        no2=ndata/2
        call realft(ndata,Data1,1,Cdat1)
        call realft(ndata,Data2,1,Cdat2)
        Cdat1(1)=cmplx(real(Cdat1(1))*real(Cdat2(1))/no2, aimag(Cdat1(1))*aimag(Cdat2(1))/no2, kind=CK)
        Cdat1(2:)=Cdat1(2:)*conjg(Cdat2(2:))/no2
        call realft(ndata,CrossCorrFFT,-1,Cdat1)
    end function getCrossCorrFFT

    !> \brief
    !> Return the cross-correlation of the two input *weighted* data vectors,
    !> (including any user-supplied zero padding), computed via Fast-Fourier Transform.
    !>
    !> @param[in]   ndata1      : The length of the first input array `Data1`.
    !> @param[in]   ndata2      : The length of the second input array `Data2`.
    !> @param[in]   paddedLen   : The length by which the input data vectors must be extended and padded.
    !> @param[in]   Data1       : The first array.
    !> @param[in]   Data2       : The second array.
    !> @param[in]   Weight1     : The weights of the elements in the first array.
    !> @param[in]   Weight2     : The weights of the elements in the second array.
    !>
    !> \return
    !> `CrossCorrFFT` : A vector of length `paddedLen` containing the cross-correlation.
    function getCrossCorrFFTweighted(nData1,nData2,paddedLen,Data1,Data2,Weight1,Weight2) result(CrossCorrFFT)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getCrossCorrFFTweighted
#endif
        use Constants_mod, only: IK, RK, CK
        implicit none
        integer(IK), intent(in)                 :: nData1, nData2, paddedLen
        real(RK), intent(inout)                 :: Data1(nData1),Data2(nData2)
        integer(IK), intent(in), optional       :: Weight1(nData1), Weight2(nData2)
        real(RK)                                :: CrossCorrFFT(paddedLen)
        character(len=*), parameter             :: PROCEDURE_NAME = MODULE_NAME//"@getCrossCorrFFTweighted()"
        complex(CK), dimension(paddedLen/2)     :: Cdat1,Cdat2
        integer(IK)                             :: no2
        if (iand(paddedLen,paddedLen-1)/=0) then
            write(*,*) PROCEDURE_NAME//": paddedLen must be a power of 2."
            error stop
        end if
        no2 = paddedLen / 2
        call realftWeighted(ndata1,paddedLen/4,Data1,Cdat1,Weight1)
        call realftWeighted(ndata2,paddedLen/4,Data2,Cdat2,Weight2)
        Cdat1(1)=cmplx(real(Cdat1(1))*real(Cdat2(1))/no2, aimag(Cdat1(1))*aimag(Cdat2(1))/no2, kind=CK)
        Cdat1(2:)=Cdat1(2:)*conjg(Cdat2(2:))/no2
        call realft(paddedLen,CrossCorrFFT,-1,Cdat1)
    end function getCrossCorrFFTweighted

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine realft(ndata,Data,isign,Zdata)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: realft
#endif
        use Constants_mod, only: IK, RK, CK
        use Misc_mod, only: zroots_unity
        implicit none
        integer(IK), intent(in)             :: ndata, isign
        real(RK), intent(inout)             :: Data(ndata)
        complex(CK), optional, target       :: Zdata(ndata/2)
        integer(IK)                         :: nh,nq
        complex(CK), dimension(ndata/4)     :: w
        complex(CK), dimension(ndata/4-1)   :: h1,h2
        complex(CK), dimension(:), pointer  :: Cdata
        complex(CK)                         :: z
        real(RK)                            :: c1=0.5_RK,c2
        nh = ndata / 2
        nq = ndata / 4
        if (present(Zdata)) then
            Cdata=>Zdata
            if (isign == 1) Cdata=cmplx(Data(1:ndata-1:2),Data(2:ndata:2),kind=CK)
        else
            allocate(Cdata(nh))
            Cdata=cmplx(Data(1:ndata-1:2),Data(2:ndata:2),kind=CK)
        end if
        if (isign == 1) then
            c2=-0.5_RK
            call four1(nh,Cdata,+1)
        else
            c2=0.5_RK
        end if
        w=zroots_unity(sign(ndata,isign),nq)
        w=cmplx(-aimag(w),real(w),kind=CK)
        h1=c1*(Cdata(2:nq)+conjg(Cdata(nh:nq+2:-1)))
        h2=c2*(Cdata(2:nq)-conjg(Cdata(nh:nq+2:-1)))
        Cdata(2:nq)=h1+w(2:nq)*h2
        Cdata(nh:nq+2:-1)=conjg(h1-w(2:nq)*h2)
        z=Cdata(1)
        if (isign == 1) then
            Cdata(1)=cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=CK)
        else
            Cdata(1)=cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=CK)
            call four1(nh,Cdata,-1)
        end if
        if (present(Zdata)) then
            if (isign /= 1) then
                Data(1:ndata-1:2)=real(Cdata)
                Data(2:ndata:2)=aimag(Cdata)
            end if
        else
            Data(1:ndata-1:2)=real(Cdata)
            Data(2:ndata:2)=aimag(Cdata)
            nullify(Cdata)
        end if
    end subroutine realft

    subroutine realftWeighted(ndata,paddedLenQuarter,Data,Zdata,Weight)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: realftWeighted
#endif
        use Constants_mod, only: IK, RK, CK
        use Misc_mod, only: zroots_unity
        implicit none
        integer(IK), intent(in)                     :: ndata, paddedLenQuarter
        integer(IK), intent(in), optional           :: Weight(ndata)
        real(RK), intent(inout)                     :: Data(ndata)
        complex(CK)                                 :: Zdata(2*paddedLenQuarter)
        integer(IK)                                 :: paddedLenHalf, id, idEnd, iw, counter
        complex(CK), dimension(paddedLenQuarter)    :: w
        complex(CK), dimension(paddedLenQuarter-1)  :: h1,h2
        real(RK)                                    :: zReal,zImag
        real(RK)                                    :: c1 = 0.5_RK, c2
        paddedLenHalf = 2 * paddedLenQuarter
        if (present(Weight)) then
            iw = 1
            counter = 0
            loopOverData: do id = 1, ndata
                loopOverWeight: do
                    if (iw>Weight(id)) then
                        iw = 1
                        cycle loopOverData
                    elseif (iw==Weight(id)) then
                        counter = counter + 1
                        if (id==ndata) then
                            Zdata(counter) = cmplx( Data(id) , 0._RK , kind=CK )
                            exit loopOverData
                        else
                            Zdata(counter) = cmplx( Data(id) , Data(id+1) , kind=CK )
                            iw = 2
                            cycle loopOverData
                        end if
                    else
                        counter = counter + 1
                        Zdata(counter) = cmplx( Data(id) , Data(id) , kind=CK )
                        iw = iw + 2
                        cycle loopOverWeight
                    end if
                end do loopOverWeight
            end do loopOverData
            Zdata(counter+1:paddedLenHalf) = cmplx( 0._RK , 0._RK , kind=CK )
        else
            if (mod(ndata,2)==0) then
                idEnd = ndata / 2
            else
                idEnd = (ndata-1) / 2
            end if
            do concurrent(id=1:idEnd)
                Zdata(id) = cmplx( Data(2*id-1) , Data(2*id) , kind=CK )
            end do
            Zdata(idEnd+1:paddedLenHalf) = cmplx( 0._RK , 0._RK , kind=CK )
        end if
        c2=-0.5_RK
        call four1(paddedLenHalf,Zdata,1_IK)
        w=zroots_unity(sign(2*paddedLenHalf,1_IK),paddedLenQuarter)
        w=cmplx(-aimag(w),real(w),kind=CK)
        h1=c1*(Zdata(2:paddedLenQuarter)+conjg(Zdata(paddedLenHalf:paddedLenQuarter+2:-1)))
        h2=c2*(Zdata(2:paddedLenQuarter)-conjg(Zdata(paddedLenHalf:paddedLenQuarter+2:-1)))
        Zdata(2:paddedLenQuarter)=h1+w(2:paddedLenQuarter)*h2
        Zdata(paddedLenHalf:paddedLenQuarter+2:-1)=conjg(h1-w(2:paddedLenQuarter)*h2)
        zReal = real(Zdata(1)); zImag = aimag(Zdata(1))
        Zdata(1) = cmplx(zReal+zImag,zReal-zImag,kind=CK)
    end subroutine realftWeighted

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine four1(ndata,Data,isign)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: four1
#endif
        use Constants_mod, only: IK, RK, CK, TWOPI
        use Misc_mod, only: arth
        implicit none
        integer(IK), intent(in)    :: ndata, isign
        complex(CK), intent(inout) :: Data(ndata)
        complex(CK), allocatable   :: dat(:,:), temp(:,:)
        complex(CK), allocatable   :: w(:),wp(:)
        real(RK)   , allocatable   :: theta(:)
        integer(IK) :: m1,m2,j
        m1=2**ceiling(0.5_RK*log(real(ndata,RK))/0.693147_RK)
        m2=ndata/m1
        allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
        dat=reshape(Data,shape(dat))
        call fourrow(dat,isign)
        theta=arth(0,isign,m1)*TWOPI/ndata
        wp=cmplx(-2.0_RK*sin(0.5_RK*theta)**2,sin(theta),kind=CK)
        w=cmplx(1.0_RK,0.0_RK,kind=CK)
        do j=2,m2
            w=w*wp+w
            dat(:,j)=dat(:,j)*w
        end do
        temp=transpose(dat)
        call fourrow(temp,isign)
        Data=reshape(temp,shape(Data))
        deallocate(dat,w,wp,theta,temp)
    end subroutine four1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine fourrow(Data,isign)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: fourrow
#endif
        use Constants_mod, only: IK, RK, CK, PI
        use Misc_mod, only: swap
        implicit none
        complex(CK), dimension(:,:), intent(inout) :: Data
        integer(IK), intent(in)                    :: isign
        integer(IK)                                :: n,i,istep,j,m,mmax,n2
        real(RK)                                   :: theta
        complex(CK), dimension(size(Data,1))       :: temp
        complex(CK)                                :: w,wp
        complex(CK)                                :: ws
        n=size(Data,2)
        n2=n/2
        j=n2
        do i=1,n-2
            if (j > i) call swap(Data(:,j+1),Data(:,i+1))
            m=n2
            do
                if (m < 2 .or. j < m) exit
                j=j-m
                m=m/2
            end do
            j=j+m
        end do
        mmax=1
        do
            if (n <= mmax) exit
            istep=2*mmax
            theta=PI/(isign*mmax)
            wp=cmplx(-2.0_RK*sin(0.5_RK*theta)**2,sin(theta),kind=CK)
            w=cmplx(1.0_RK,0.0_RK,kind=CK)
            do m=1,mmax
                ws=w
                do i=m,n,istep
                    j=i+mmax
                    temp=ws*Data(:,j)
                    Data(:,j)=Data(:,i)-temp
                    Data(:,i)=Data(:,i)+temp
                end do
                w=w*wp+w
            end do
            mmax=istep
        end do
    end subroutine fourrow

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! The following are for the AutoCorrelation computation using the conventional definition of correlation.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getInverseSumNormedDataSq(nd,np,NormedData) result(InverseSumNormedDataSq)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getInverseSumNormedDataSq
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK), intent(in)     :: nd,np
        real(RK)   , intent(in)     :: NormedData(nd,np)
        real(RK)                    :: InverseSumNormedDataSq(nd)
        integer(IK)                 :: ip
        InverseSumNormedDataSq = 0._RK
        do ip = 1, np
            InverseSumNormedDataSq = InverseSumNormedDataSq + NormedData(1:nd,ip)**2
        end do
        InverseSumNormedDataSq = 1._RK / InverseSumNormedDataSq
    end function getInverseSumNormedDataSq

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure subroutine getAutoCorrSlow(nd,np,NormedData,nlag,Lag,AutoCorr,InverseSumNormedDataSq)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getAutoCorrSlow
#endif

        use Constants_mod, only: IK, RK
        implicit none
        integer(IK), intent(in)             :: nd,np,nlag,Lag(nlag)
        real(RK)   , intent(in)             :: NormedData(nd,np)
        real(RK)   , intent(in), optional   :: InverseSumNormedDataSq(nd)
        real(RK)   , intent(out)            :: AutoCorr(nd,nlag)
        real(RK)                            :: InverseSumNormedDataSqComputed(nd)
        integer(IK)                         :: ip, ilag

        if (any(Lag>np-1_IK)) then
            AutoCorr    = -huge(1._RK)
            return
        end if

        if (present(InverseSumNormedDataSq)) then
            InverseSumNormedDataSqComputed = InverseSumNormedDataSq
        else
            InverseSumNormedDataSqComputed = 0._RK
            do ip = 1, np
                InverseSumNormedDataSqComputed = InverseSumNormedDataSqComputed + NormedData(1:nd,ip)**2
            end do
            InverseSumNormedDataSqComputed = 1._RK / InverseSumNormedDataSqComputed
        end if

        ! Now compute the non-normalized covariances

        do ilag = 1, nlag
            AutoCorr(1:nd,ilag) = 0._RK
            do ip = 1, np-Lag(ilag)
                AutoCorr(1:nd,ilag) = AutoCorr(1:nd,ilag) + NormedData(1:nd,ip) * NormedData(1:nd,ip+Lag(ilag))
            end do
            AutoCorr(1:nd,ilag) = AutoCorr(1:nd,ilag) * InverseSumNormedDataSqComputed
        end do

    end subroutine getAutoCorrSlow

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module CrossCorr_mod