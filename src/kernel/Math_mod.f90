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

!>  \brief This module contains mathematical procedures.
!>  @author Amir Shahmoradi

module Math_mod

    implicit none

    character(*), parameter :: MODULE_NAME = "@Math_mod"

    interface getCumSum
        module procedure :: getCumSum_IK, getCumSum_RK
    end interface getCumSum

    interface getCumSumReverse
        module procedure :: getCumSumReverse_IK, getCumSumReverse_RK
    end interface getCumSumReverse

    interface getLogSumExp
        module procedure :: getLogSumExp_RK, getLogSumExp_CK
    end interface getLogSumExp

    interface getLogSubExp
        module procedure :: getLogSubExp_RK
    end interface getLogSubExp

    interface getLogEggBox
        module procedure :: getLogEggBoxSD_RK, getLogEggBoxMD_RK
        module procedure :: getLogEggBoxSD_CK, getLogEggBoxMD_CK
    end interface getLogEggBox

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the distance squared between the two input points.
    !>
    !> @param[in]   nd      : The size of the input vectors.
    !> @param[in]   Point1  : The first point.
    !> @param[in]   Point2  : The second point.
    !>
    !> \return
    !> `distanceSq` : The distance squared.
    pure function getDistanceSq(nd,Point1,Point2) result(distanceSq)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getDistanceSq
#endif
        use Constants_mod, only: IK, RK
        integer(IK) , intent(in)    :: nd
        real(RK)    , intent(in)    :: Point1(nd),Point2(nd)
        real(RK)                    :: distanceSq
        integer(IK)                 :: i
        distanceSq = 0._RK
        do i = 1, nd
            distanceSq = distanceSq + (Point2(i)-Point1(i))**2
        end do
    end function getDistanceSq


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the correlation coefficient (`-1 < corCoef < 1`) corresponding to the input Fisher z-transformation.
    !>
    !> @param[in]   fisherTrans : The Fisher transformation.
    !>
    !> \return
    !> `corCoef` : The correlation coefficient corresponding to the input Fisher transformation.
    pure elemental function getCorCeofFromFisherTrans(fisherTrans) result(corCoef)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getCorCeofFromFisherTrans
#endif
        use Constants_mod, only: RK
        real(RK), intent(in) :: fisherTrans
        real(RK)             :: corCoef !, twiceFisherTrans
        corCoef = tanh(fisherTrans)
        !twiceFisherTrans = 2._RK * corCoef
        !fisherTrans = (exp(twiceFisherTrans)-1._RK) / (exp(twiceFisherTrans)+1._RK)
    end function getCorCeofFromFisherTrans


    !> \brief
    !> Return Fisher z-transformation of an input correlation coefficient (`-1 < corCoef < 1`).
    !>
    !> @param[in]   corCoef : The input correlation coefficient.
    !>
    !> \return
    !> `fisherTrans` : The output Fisher transformation.
    pure elemental function getFisherTransFromCorCoef(corCoef) result(fisherTrans)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getFisherTransFromCorCoef
#endif
        use Constants_mod, only: RK
        real(RK), intent(in) :: corCoef
        real(RK)             :: fisherTrans
        fisherTrans = atanh(corCoef)
        !fisherTrans = 0.5_RK * log(1._RK+corCoef) / (1._RK-corCoef)
    end function getFisherTransFromCorCoef

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the cumulative sum of the input integer array.
    !>
    !> @param[in]   vecLen  : The length of the input vector.
    !> @param[in]   Vec     : The input vector.
    !>
    !> \return
    !> `CumSum` : An integer array of length `vecLen` representing the cumulative sum.
    !>
    !> \remark
    !> The first element of `CumSum` is the same as the first element of `Vec`.
    pure function getCumSum_IK(vecLen,Vec) result(CumSum)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSum_IK
#endif
        use Constants_mod, only: IK
        integer(IK), intent(in) :: vecLen
        integer(IK), intent(in) :: Vec(vecLen)
        integer(IK)             :: CumSum(vecLen)
        integer(IK)             :: i
        CumSum(1) = Vec(1)
        do i = 2, vecLen
            CumSum(i) = CumSum(i-1) + Vec(i)
        end do
    end function getCumSum_IK

    !> \brief
    !> Return the cumulative sum of the input real array.
    !> @param[in]   vecLen  : The length of the input vector.
    !> @param[in]   Vec     : The input vector.
    !>
    !> \return
    !> `CumSum` : A real array of length `vecLen` representing the cumulative sum.
    !>
    !> \remark
    !> The first element of `CumSum` is the same as the first element of `Vec`.
    pure function getCumSum_RK(vecLen,Vec) result(CumSum)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSum_RK
#endif
        use Constants_mod, only: IK, RK
        integer(IK), intent(in) :: vecLen
        real(RK)   , intent(in) :: Vec(vecLen)
        real(RK)                :: CumSum(vecLen)
        integer(IK)             :: i
        CumSum(1) = Vec(1)
        do i = 2, vecLen
            CumSum(i) = CumSum(i-1) + Vec(i)
        end do
    end function getCumSum_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the cumulative sum of the input integer array, in the backward direction.
    !> @param[in]   vecLen  : The length of the input vector.
    !> @param[in]   Vec     : The input vector.
    !>
    !> \return
    !> `CumSumReverse` : An integer array of length `vecLen` representing the cumulative sum.
    !>
    !> \remark
    !> The last element of `CumSumReverse` is the same as the last element of `Vec`.
    pure function getCumSumReverse_IK(vecLen,Vec) result(CumSumReverse)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSumReverse_IK
#endif
        use Constants_mod, only: IK
        integer(IK), intent(in) :: vecLen
        integer(IK), intent(in) :: Vec(vecLen)
        integer(IK)             :: CumSumReverse(vecLen)
        integer(IK)             :: i, indx
        CumSumReverse(1) = Vec(vecLen)
        do i = vecLen-1,1,-1
            indx = vecLen - i
            CumSumReverse(indx+1) = CumSumReverse(indx) + Vec(i)
        end do
    end function getCumSumReverse_IK

    !> \brief
    !> Return the cumulative sum of the input real array, in the backward direction.
    !>
    !> @param[in]   vecLen  : The length of the input vector.
    !> @param[in]   Vec     : The input vector.
    !>
    !> \return
    !> `CumSumReverse` : A real array of length `vecLen` representing the cumulative sum.
    !>
    !> \remark
    !> The last element of `CumSumReverse` is the same as the last element of `Vec`.
    pure function getCumSumReverse_RK(vecLen,Vec) result(CumSumReverse)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSumReverse_RK
#endif
        use Constants_mod, only: IK, RK
        integer(IK), intent(in) :: vecLen
        real(RK)   , intent(in) :: Vec(vecLen)
        real(RK)                :: CumSumReverse(vecLen)
        integer(IK)             :: i, indx
        CumSumReverse(1) = Vec(vecLen)
        do i = vecLen-1,1,-1
            indx = vecLen - i
            CumSumReverse(indx+1) = CumSumReverse(indx) + Vec(i)
        end do
    end function getCumSumReverse_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return `log( exp(logValueLarger) - exp(logValueSamller) )` robustly (without overflow or underflow).
    !>
    !> @param[in]   logValueLarger  : A real number.
    !> @param[in]   logValueSamller : A real number.
    !>
    !> \return
    !> `logSubExp` : A real number.
    !>
    !> \warning
    !> The onus is on the user to ensure `logValueLarger > logValueSamller`.
    pure function getLogSubExp_RK(logValueLarger,logValueSamller) result(logSubExp)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExp_RK
#endif
        use Constants_mod, only: RK
        real(RK)   , intent(in) :: logValueLarger, logValueSamller
        real(RK)                :: logSubExp
        logSubExp = logValueLarger + log( 1._RK - exp(logValueSamller-logValueLarger) )
    end function getLogSubExp_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the logarithm of the sum of the exponential of the input real vector robustly (without overflow or underflow).
    !>
    !> @param[in]   lenLogValue : The length of the input vector.
    !> @param[in]   LogValue    : The input vector of log-values whose log-sum-exp must be computed.
    !> @param[in]   maxLogValue : The maximum of the input `LogValue` argument (optional).
    !>
    !> \return
    !> `logSumExp` : A real number.
    pure function getLogSumExp_RK(lenLogValue, LogValue, maxLogValue) result(logSumExp)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExp_RK
#endif
        use Constants_mod, only: IK, RK, LOGTINY_RK
        integer(IK) , intent(in)            :: lenLogValue
        real(RK)    , intent(in)            :: LogValue(lenLogValue)
        real(RK)    , intent(in), optional  :: maxLogValue
        real(RK)                            :: logSumExp
        real(RK)                            :: maxLogValueDefault
        real(RK)                            :: LogValueCopy(lenLogValue)
        integer(IK)                         :: i
        if (present(maxLogValue)) then
            maxLogValueDefault = maxLogValue
        else
            maxLogValueDefault = maxval(LogValue)
        end if
        LogValueCopy = LogValue - maxLogValueDefault
        do concurrent(i=1:lenLogValue)
            if (LogValueCopy(i)<LOGTINY_RK) then
                LogValueCopy(i) = 0._RK
            else
                LogValueCopy(i) = exp(LogValueCopy(i))
            end if
        end do
        logSumExp = maxLogValueDefault + log(sum(LogValueCopy))
    end function getLogSumExp_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the logarithm of the sum of the exponential of the input complex vector robustly (without overflow or underflow).
    !>
    !> @param[in]   lenLogValue : The length of the input vector.
    !> @param[in]   LogValue    : The input vector of log-values whose log-sum-exp must be computed.
    !> @param[in]   maxLogValue : The maximum of the real component of the input `LogValue` argument (optional).
    !>
    !> \return
    !> `logSumExp` : A complex number.
    pure function getLogSumExp_CK(lenLogValue, LogValue, maxLogValue) result(logSumExp)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExp_CK
#endif
        use Constants_mod, only: IK, RK, LOGTINY_RK
        integer(IK) , intent(in)            :: lenLogValue
        complex(RK) , intent(in)            :: LogValue(lenLogValue)
        complex(RK) , intent(in), optional  :: maxLogValue
        complex(RK)                         :: logSumExp
        complex(RK)                         :: LogValueCopy(lenLogValue)
        complex(RK)                         :: maxLogValueDefault
        integer(IK)                         :: i
        if (present(maxLogValue)) then
            maxLogValueDefault = maxLogValue
        else
            maxLogValueDefault = maxval(real(LogValue))
        end if
        LogValueCopy = LogValue - maxLogValueDefault
        do concurrent(i=1:lenLogValue)
            if (real(LogValueCopy(i))<LOGTINY_RK) then
                LogValueCopy(i) = 0._RK
            else
                LogValueCopy(i) = exp(LogValueCopy(i))
            end if
        end do
        logSumExp = maxLogValueDefault + log(sum(LogValueCopy))
    end function getLogSumExp_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the logarithm of the egg-box probability density function in one dimension.
    pure function getLogEggBoxSD_RK(constant,exponent,coef,point) result(logEggBox)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogEggBoxSD_RK
#endif
        use Constants_mod, only: IK, RK
        implicit none
        real(RK), intent(in) :: constant
        real(RK), intent(in) :: exponent
        real(RK), intent(in) :: coef
        real(RK), intent(in) :: point
        real(RK)             :: logEggBox
        logEggBox = exponent * log( constant + cos(coef*point) )
    end function getLogEggBoxSD_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the logarithm of the egg-box probability density function in one dimension, as a complex number.
    pure function getLogEggBoxSD_CK(constant,exponent,coef,point) result(logEggBox)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogEggBoxSD_CK
#endif
        use Constants_mod, only: IK, RK
        implicit none
        complex(RK), intent(in) :: constant
        complex(RK), intent(in) :: exponent
        complex(RK), intent(in) :: coef
        complex(RK), intent(in) :: Point
        complex(RK)             :: logEggBox
        logEggBox = exponent * log( constant + cos(coef*point) )
    end function getLogEggBoxSD_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the logarithm of the egg-box probability density function in multiple dimensions, as a real number.
    pure function getLogEggBoxMD_RK(nd,constant,exponent,Coef,Point) result(logEggBox)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogEggBoxMD_RK
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK), intent(in) :: nd
        real(RK), intent(in) :: constant
        real(RK), intent(in) :: exponent
        real(RK), intent(in) :: Coef(nd)
        real(RK), intent(in) :: Point(nd)
        real(RK)             :: logEggBox
        integer(IK)          :: i
        logEggBox = 0._RK
        do i = 1,nd
            logEggBox = logEggBox * cos(Coef(i)*Point(i))
        end do
        logEggBox = exponent * log( constant + logEggBox )
    end function getLogEggBoxMD_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the logarithm of the egg-box probability density function in multiple dimensions, as a complex number.
    pure function getLogEggBoxMD_CK(nd,constant,exponent,Coef,Point) result(logEggBox)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogEggBoxMD_CK
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK), intent(in) :: nd
        complex(RK), intent(in) :: constant
        complex(RK), intent(in) :: exponent
        complex(RK), intent(in) :: Coef(nd)
        complex(RK), intent(in) :: Point(nd)
        complex(RK)             :: logEggBox
        integer(IK)          :: i
        logEggBox = 0._RK
        do i = 1,nd
            logEggBox = logEggBox * cos(Coef(i)*Point(i))
        end do
        logEggBox = exponent * log( constant + logEggBox )
    end function getLogEggBoxMD_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the log(factorial) for a whole integer input. This is basically `logGamma( positiveInteger + 1 )`.
    !>
    !> \param[in]   positiveInteger :   The input positive integer whose log factorial must be computed.
    !>
    !> \return
    !> `logFactorial`  :   The natural logarithm of the factorial.
    !>
    !> \remark
    !> This function is mostly useful for large input integers.
    pure function getLogFactorial(positiveInteger) result(logFactorial)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogFactorial
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK), intent(in) :: positiveInteger
        integer(IK)             :: i
        real(RK)                :: logFactorial
        logFactorial = 0._RK
        do i = 2, positiveInteger
            logFactorial = logFactorial + log(real(i,kind=RK))
        end do
    end function getLogFactorial

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the factorial for a whole integer input. This is basically `Gamma( intNum + 1 )`.
    !>
    !> \param[in]   positiveInteger :   The input positive integer whose log factorial must be computed.
    !>
    !> \return
    !> `factorial`  :   The factorial.
    pure function getFactorial(positiveInteger) result(factorial)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getFactorial
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK), intent(in) :: positiveInteger
        integer(IK)             :: i
        real(RK)                :: factorial
        factorial = 1._RK
        do i = 2, positiveInteger
            factorial = factorial * i
        end do
    end function getFactorial

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the Gamma function for a half integer input.
    pure function getGammaHalfInt(halfIntNum)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaHalfInt
#endif
        use Constants_mod, only: IK, RK, SQRTPI
        implicit none
        real(RK), intent(in) :: halfIntNum
        real(RK)             :: getGammaHalfInt
        integer(IK)          :: i,k
        getGammaHalfInt = SQRTPI
        k = nint(halfIntNum-0.5_RK,kind=IK)    ! halfIntNum = k + 1/2
        do i = k+1, 2*k
            getGammaHalfInt = getGammaHalfInt * i / 4._RK
        end do
    end function getGammaHalfInt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! see for example: http://math.stackexchange.com/questions/606184/volume-of-n-dimensional-ellipsoid
    ! see for example: https://en.wikipedia.org/wiki/Volume_of_an_n-ball
    ! see for example: https://en.wikipedia.org/wiki/Particular_values_of_the_Gamma_function
    ! getEllVolCoef = PI^(nd/2) / gamma(nd/2+1) where n is just a positive integer
    !
    !> \brief
    !> Return the coefficient of the volume of an `nd`-dimensional ellipsoid.
    !>
    !> param[in]    nd : The number of dimensions.
    !>
    !> \return
    !> `ellVolCoef` : The coefficient of the volume of an `nd`-dimensional ellipsoid.
    pure function getEllVolCoef(nd) result(ellVolCoef)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getEllVolCoef
#endif
        use Constants_mod, only: IK, RK, PI
        implicit none
        integer(IK), intent(in) :: nd
        integer(IK)             :: i,k
        real(RK), parameter     :: FOUR_PI = PI * 4._RK
        real(RK)                :: ellVolCoef
        if (mod(nd,2_IK)==0_IK) then ! nd is even

            ellVolCoef = PI
            do i = 2_IK, nd / 2_IK
                ellVolCoef = ellVolCoef * PI / i    ! nd = 2k ; ellVolCoef = PI^(k) / Factorial(k)
            end do

        else ! nd is an odd integer

            ! nd = 2k-1 ; gamma(nd/2 + 1) = gamma(k + 1/2) ; gamma(k+1/2) = sqrt(PI) * (2k)! / (4^k * k!)
            k = (nd + 1_IK) / 2_IK

            ! This is to avoid an extra unnecessary division of ellVolCoef by PI
            ellVolCoef = 4._RK / (k + 1_IK)

            do i = k+2_IK, 2_IK*k
                ! ellVolCoef = PI^(k-1/2) / gamma(k+1/2) = PI^(k+1) * 4^k * k! / (2k)!
                ellVolCoef = ellVolCoef * FOUR_PI / i
            end do

        end if
    end function getEllVolCoef

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! see for example: http://math.stackexchange.com/questions/606184/volume-of-n-dimensional-ellipsoid
    ! see for example: https://en.wikipedia.org/wiki/Volume_of_an_n-ball
    ! see for example: https://en.wikipedia.org/wiki/Particular_values_of_the_Gamma_function
    ! getEllVolCoef = PI^(nd/2) / gamma(nd/2+1) where n is just a positive integer
    !
    !> \brief
    !> Return the coefficient of the volume of an `nd`-dimensional ellipsoid.
    !>
    !> param[in]    nd : The number of dimensions.
    !>
    !> \return
    !> `logEllVolCoef` : The natural logarithm of the volume of an `nd`-dimensional unit-ball.
    pure function getLogEllVolCoef(nd) result(logEllVolCoef)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getEllVolCoef
#endif
        use Constants_mod, only: IK, RK, PI
        implicit none
        integer(IK), intent(in) :: nd
        integer(IK)             :: i,k
        real(RK), parameter     :: FOUR_PI = PI * 4._RK
        real(RK)                :: logEllVolCoef
        if (mod(nd,2_IK)==0_IK) then ! nd is even

            logEllVolCoef = PI
            do i = 2_IK, nd / 2_IK
                logEllVolCoef = logEllVolCoef * PI / i    ! nd = 2k ; ellVolCoef = PI^(k) / Factorial(k)
            end do

        else ! nd is an odd integer

            ! nd = 2k-1 ; gamma(nd/2 + 1) = gamma(k + 1/2) ; gamma(k+1/2) = sqrt(PI) * (2k)! / (4^k * k!)
            k = (nd + 1_IK) / 2_IK

            ! This is to avoid an extra unnecessary division of ellVolCoef by PI
            logEllVolCoef = 4._RK / (k + 1_IK)

            do i = k+2_IK, 2_IK*k
                ! ellVolCoef = PI^(k-1/2) / gamma(k+1/2) = PI^(k+1) * 4^k * k! / (2k)!
                logEllVolCoef = logEllVolCoef * FOUR_PI / i
            end do

        end if
        logEllVolCoef = log(logEllVolCoef)
    end function getLogEllVolCoef

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the logarithm of the volume of an `nd`-dimensional ball of unit-radius.
    !>
    !> param[in]    nd  : The number of dimensions.
    !>
    !> \return
    !> `logVolUnitBall` : The logarithm of the volume of an `nd`-dimensional ball of unit-radius.
    !>
    !> \brief
    !> This routine is an exact replacement for [getLogEllVolCoef](@ref getlogellvolcoef).
    !> However, it is, on average, 10%-15% faster than [getLogEllVolCoef](@ref getlogellvolcoef) based
    !> on benchmarks with Intel ifort 19.4 with all optimization flags on.
    !
    ! see for example: http://math.stackexchange.com/questions/606184/volume-of-n-dimensional-ellipsoid
    ! see for example: https://en.wikipedia.org/wiki/Volume_of_an_n-ball
    ! see for example: https://en.wikipedia.org/wiki/Particular_values_of_the_Gamma_function
    ! getEllVolCoef = PI^(nd/2) / gamma(nd/2+1) where n is just a positive integer
    pure function getLogVolUnitBall(nd) result(logVolUnitBall)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogVolUnitBall
#endif
        use Constants_mod, only: IK, RK, PI
        implicit none
        integer(IK) , intent(in)    :: nd
        real(RK)    , parameter     :: LOG_PI = log(PI)
        integer(IK)                 :: ndHalfInteger
        real(RK)                    :: logVolUnitBall
        real(RK)                    :: ndHalfReal
        if (mod(nd,2_IK)==0_IK) then ! nd is even
            ndHalfInteger = nd / 2_IK
            logVolUnitBall = ndHalfInteger*LOG_PI - getLogFactorial(ndHalfInteger) ! nd = 2k ; logVolUnitBall = PI^(k) / Factorial(k)
        else ! nd is an odd integer
            ndHalfReal = 0.5_RK * real(nd,kind=RK)
            logVolUnitBall = ndHalfReal * LOG_PI - log_gamma(ndHalfReal+1._RK)
        end if
    end function getLogVolUnitBall

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the logarithm of the volume of an `nd`-dimensional hyper-ellipsoid.
    !>
    !> param[in]    nd                  :   The number of dimensions.
    !> param[in]    logSqrtDetCovMat    :   The logarithm of the square root of the determinant of the representative covariance matrix of the ellipsoid.
    !>
    !> \return
    !> `logVolEllipsoid` : The logarithm of the volume of an `nd`-dimensional hyper-ellipsoid.
    !
    ! see for example: https://math.stackexchange.com/questions/332391/volume-of-hyperellipsoid/332434
    ! see for example: https://math.stackexchange.com/questions/2854930/volume-of-an-n-dimensional-ellipsoid
    pure function getLogVolEllipsoid(nd,logSqrtDetCovMat) result(logVolEllipsoid)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogVolEllipsoid
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK) , intent(in)    :: nd
        real(RK)    , intent(in)    :: logSqrtDetCovMat
        real(RK)                    :: logVolEllipsoid
        logVolEllipsoid = getLogVolUnitBall(nd) + logSqrtDetCovMat
    end function getLogVolEllipsoid

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the volume of a list of `nEllipsoid` `nd`-dimensional hyper-ellipsoids.
    !>
    !> @param[in]   nd                  :   The number of dimensions of the ellipsoids.
    !> @param[in]   nEllipsoid          :   The number of ellipsoids.
    !> @param[in]   LogSqrtDetCovMat    :   The vector of natural logarithms of the square roots of the determinants of the covariance matrices of the ellipsoids.
    !>
    !> \return
    !> `logVolEllipsoid` : The logarithm of the volume of an `nd`-dimensional hyper-ellipsoid.
    !>
    !> \remark
    !> There is really no reason to use this function for HPC solutions, as it can be likely done more efficiently.
    pure function getLogVolEllipsoids(nd,nEllipsoid,LogSqrtDetCovMat) result(LogVolEllipsoids)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogVolEllipsoids
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK)                 :: i
        integer(IK) , intent(in)    :: nd
        integer(IK) , intent(in)    :: nEllipsoid
        real(RK)    , intent(in)    :: LogSqrtDetCovMat(nEllipsoid)
        real(RK)                    :: LogVolEllipsoids(nEllipsoid)
        real(RK)                    :: logVolUnitBall
        logVolUnitBall = getLogVolUnitBall(nd)
        do i = 1, nEllipsoid
            LogVolEllipsoids(i) = logVolUnitBall + LogSqrtDetCovMat(i)
        end do
    end function getLogVolEllipsoids

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the lower incomplete gamma function for the specified `exponent` and upper limit.
    !> `tolerance` represents the relative accuracy.
    !>
    !> \warning
    !> Do not set tolerance to a number larger than 1, or else the code will crash spectacularly.
    !>
    !> \remark
    !> `logGammaExponent = log_gamma(exponent)`
    !>
    !> \warning
    !> On input, both `exponent` and `upperLim` must be positive, otherwise, a `lowerGamma = -infinity` will be returned to signal error.
    ! This algorithm borrows from the `gammp` implementation of the Numerical Recipes by Press et al 1992.
    pure function getLowerGamma(exponent,logGammaExponent,upperLim,tolerance) result(lowerGamma)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLowerGamma
#endif
        use Constants_mod, only: IK, RK, HUGE_RK
        implicit none
        real(RK), intent(in)            :: exponent, logGammaExponent, upperLim
        real(RK), intent(in), optional  :: tolerance
        real(RK)                        :: lowerGamma
        if ( upperLim < 0._RK .or. exponent <= 0._RK ) then
            lowerGamma = -HUGE_RK
            return
        elseif ( upperLim < exponent + 1._RK ) then
            lowerGamma = getGammaSeries(exponent,log_gamma(exponent),upperLim,tolerance)
        else
            lowerGamma = 1._RK - getGammaContFrac(exponent,logGammaExponent,upperLim,tolerance)
        end if
    end function getLowerGamma

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the lower incomplete gamma function for the specified exponent and upper limit.
    !>
    !> `tolerance` represents the relative accuracy.
    !>
    !> \warning
    !> Do not set `tolerance` to a number larger than 1, or else the code will crash spectacularly.
    !>
    !> \warning
    !> On input, both `exponent` and `lowerLim` must be positive, otherwise, a `upperGamma = -infinity` will be returned to signal error.
    ! This algorithm borrows from the gammq implementation of the Numerical Recipes by Press et al 1992.
    pure function getUpperGamma(exponent,logGammaExponent,lowerLim,tolerance) result(upperGamma)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getUpperGamma
#endif
        use Constants_mod, only: IK, RK, HUGE_RK
        implicit none
        real(RK), intent(in)            :: exponent, logGammaExponent, lowerLim
        real(RK), intent(in), optional  :: tolerance
        real(RK)                        :: upperGamma
        if ( lowerLim < 0._RK .or. exponent <= 0._RK ) then
            upperGamma = -HUGE_RK
            return
        elseif (lowerLim < exponent + 1._RK) then
            upperGamma = 1._RK - getGammaSeries(exponent,logGammaExponent,lowerLim,tolerance)
        else
            upperGamma = getGammaContFrac(exponent,logGammaExponent,lowerLim,tolerance)
        end if
    end function getUpperGamma

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the lower incomplete gamma function `P(exponent, upperLim)` evaluated by its series representation as `gammaSeries`.
    !>
    !> \param[in]   exponent            :   The exponent in lower incomplete gamma function `P(exponent, upperLim)`.
    !> \param[in]   logGammaExponent    :   The input `log( gamma(exponent) )`.
    !> \param[in]   upperLim            :   The upper limit in lower incomplete gamma function `P(exponent, upperLim)`.
    !> \param[in]   tolerance           :   The input relative accuracy.
    !>
    !> \warning
    !> Do not set `tolerance` to a number larger than 1, or else the code will crash spectacularly.
    !>
    !> \warning
    !> If the algorithm fails to converge, `gammaSeries` will be set to negative infinity on output to signal error.
    ! This algorithm borrows from the gser implementation of the Numerical Recipes by Press et al 1992.
    pure function getGammaSeries(exponent,logGammaExponent,upperLim,tolerance) result(gammaSeries)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaSeries
#endif
        use Constants_mod, only: IK, RK, HUGE_RK
        implicit none
        real(RK), intent(in)            :: exponent, logGammaExponent, upperLim
        real(RK), intent(in), optional  :: tolerance
        real(RK)                        :: gammaSeries
        integer(IK), parameter          :: ITMAX = 100
        real(RK), parameter             :: EPS_DEFAULT = epsilon(upperLim)
        real(RK)                        :: ap,del,summ,eps
        integer(IK)                     :: iter
        if (present(tolerance)) then
            eps = tolerance
        else
            eps = EPS_DEFAULT
        end if
        if (upperLim == 0._RK) then
            gammaSeries = 0._RK
            return
        end if
        ap = exponent
        summ = 1.0_RK / exponent
        del = summ
        do iter = 1, ITMAX
            ap = ap + 1.0_RK
            del = del*upperLim / ap
            summ = summ + del
            if (abs(del) < abs(summ)*eps) exit
        end do
        if (iter>ITMAX) then
            gammaSeries = -HUGE_RK
        else
            gammaSeries = summ * exp( exponent*log(upperLim) - upperLim - logGammaExponent )
        end if
    end function getGammaSeries

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the incomplete gamma function `Q(exponent, lowerLim)` evaluated by its continued fraction representation as `gammaContFrac`.
    !>
    !> @param[in]   logGammaExponent    :   This is the `log( gamma(exponent) )`.
    !> @param[in]   tolerance           :   Represents the relative accuracy (optional).
    !> @param[in]   exponent            :   The exponent.
    !> @param[in]   lowerLim            :   The lower limit.
    !>
    !> \warning
    !> If the algorithm fails to converge, `gammaContFrac` will be set to negative infinity on output to signal error.
    ! This algorithm borrows from the gcf implementation of the Numerical Recipes by Press et al 1992.
    pure function getGammaContFrac(exponent,logGammaExponent,lowerLim,tolerance) result(gammaContFrac)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaContFrac
#endif
        use Constants_mod, only: IK, RK, HUGE_RK
        implicit none
        real(RK), intent(in)            :: exponent, logGammaExponent, lowerLim
        real(RK), optional, intent(in)  :: tolerance
        real(RK)                        :: gammaContFrac
        integer(IK), parameter          :: ITMAX = 100_IK                               !< The maximum allowed number of iterations.
        real(RK), parameter             :: EPS_DEFAULT = epsilon(lowerLim)              !< The default relative accuracy.
        real(RK), parameter             :: FPMIN_DEFAULT = tiny(lowerLim) / EPS_DEFAULT !< A number near the smallest representable floating-point number.
        real(RK)                        :: eps, fpmin
        real(RK)                        :: an,b,c,d,del,h
        integer(IK)                     :: iter
        if (lowerLim == 0._RK) then
            gammaContFrac = 1._RK
            return
        end if
        if (present(tolerance)) then
            eps = tolerance
            fpmin = tiny(lowerLim) / eps
        else
            eps = EPS_DEFAULT
            fpmin = FPMIN_DEFAULT
        end if
        b = lowerLim + 1._RK - exponent
        c = 1._RK / fpmin
        d = 1._RK / b
        h = d
        do iter = 1, ITMAX
            an = -iter * (iter-exponent)
            b = b + 2.0_RK
            d = an*d + b
            if (abs(d) < fpmin) d = fpmin
            c = b + an/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_RK / d
            del = d*c
            h = h*del
            if (abs(del-1._RK) <= eps) exit
        end do
       !if (iter > ITMAX) call nrerror('exponent too large, ITMAX too small in getGammaContFrac')
        if (iter>ITMAX) then
            gammaContFrac = -HUGE_RK
        else
            gammaContFrac = exp(exponent*log(lowerLim) - lowerLim - logGammaExponent)*h
        end if
    end function getGammaContFrac

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Math_mod