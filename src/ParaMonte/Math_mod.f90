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

    ! returns the correlation coefficient (-1<corCoef<1) corresponding to the input Fisher z-transformation.
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

    ! returns the correlation coefficient (-1<corCoef<1) corresponding to the input Fisher z-transformation.
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

    ! returns Fisher z-transformation of an input correlation coefficient (-1<corCoef<1).
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

    pure function getCumSumReverse_IK(vecLen,Vec) result(CumSum)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSumReverse_IK
#endif
        use Constants_mod, only: IK
        integer(IK), intent(in) :: vecLen
        integer(IK), intent(in) :: Vec(vecLen)
        integer(IK)             :: CumSum(vecLen)
        integer(IK)             :: i, indx
        CumSum(1) = Vec(vecLen)
        do i = vecLen-1,1,-1
            indx = vecLen - i
            CumSum(indx+1) = CumSum(indx) + Vec(i)
        end do
    end function getCumSumReverse_IK

    pure function getCumSumReverse_RK(vecLen,Vec) result(CumSum)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getCumSumReverse_RK
#endif
        use Constants_mod, only: IK, RK
        integer(IK), intent(in) :: vecLen
        real(RK)   , intent(in) :: Vec(vecLen)
        real(RK)                :: CumSum(vecLen)
        integer(IK)             :: i, indx
        CumSum(1) = Vec(vecLen)
        do i = vecLen-1,1,-1
            indx = vecLen - i
            CumSum(indx+1) = CumSum(indx) + Vec(i)
        end do
    end function getCumSumReverse_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogSubExp_RK(logFunc1,logFunc2) result(logSubExp)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSubExp_RK
#endif
        ! returns log( exp(logFunc1) - exp(logFunc2) ) robustly
        ! onus is on the user to ensure logFunc1 > logFunc2. To ensure wrong input does not go unnoticed
        use Constants_mod, only: RK
        real(RK)   , intent(in) :: logFunc1, logFunc2
        real(RK)                :: logSubExp
        logSubExp = logFunc1 + log( 1._RK - exp(logFunc2-logFunc1) )
    end function getLogSubExp_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogSumExp_RK(nFunc,LogFunc)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExp_RK
#endif
        use Constants_mod, only: IK, RK, LOGTINY_RK
        integer(IK), intent(in) :: nFunc
        real(RK)   , intent(in) :: LogFunc(nFunc)
        real(RK)                :: getLogSumExp_RK
        real(RK)                :: LogFuncCopy(nFunc)
        real(RK)                :: normFac
        integer(IK)             :: i
        normFac = maxval(LogFunc)
        LogFuncCopy = LogFunc - normFac
        do concurrent(i=1:nFunc)
            if (LogFuncCopy(i)<LOGTINY_RK) then
                LogFuncCopy(i) = 0._RK
            else
                LogFuncCopy(i) = exp(LogFuncCopy(i))
            end if
        end do
        getLogSumExp_RK = normFac + log(sum(LogFuncCopy))
    end function getLogSumExp_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogSumExp_CK(nFunc,LogFunc)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogSumExp_CK
#endif
        use Constants_mod, only: IK, RK, LOGTINY_RK
        integer(IK), intent(in) :: nFunc
        complex(RK), intent(in) :: LogFunc(nFunc)
        complex(RK)             :: getLogSumExp_CK
        complex(RK)             :: LogFuncCopy(nFunc)
        complex(RK)             :: normFac
        integer(IK)             :: i
        normFac = maxval(real(LogFunc))
        LogFuncCopy = LogFunc - normFac
        do concurrent(i=1:nFunc)
            if (real(LogFuncCopy(i))<LOGTINY_RK) then
                LogFuncCopy(i) = 0._RK
            else
                LogFuncCopy(i) = exp(LogFuncCopy(i))
            end if
        end do
        getLogSumExp_CK = normFac + log(sum(LogFuncCopy))
    end function getLogSumExp_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogEggBoxSD_RK(constant,exponent,coef,point) result(logEggBox)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogEggBoxSD_RK
#endif
        ! This simple function returns the log (in natural Neper base) of the EggBox function for nd > 1.
        use Constants_mod, only: IK, RK
        implicit none
        real(RK), intent(in) :: constant
        real(RK), intent(in) :: exponent
        real(RK), intent(in) :: Coef
        real(RK), intent(in) :: Point
        real(RK)             :: logEggBox
        logEggBox = exponent * log( constant + cos(coef*point) )
    end function getLogEggBoxSD_RK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogEggBoxSD_CK(constant,exponent,coef,point) result(logEggBox)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogEggBoxSD_CK
#endif
        ! This simple function returns the log (in natural Neper base) of the EggBox function for nd > 1.
        use Constants_mod, only: IK, RK
        implicit none
        complex(RK), intent(in) :: constant
        complex(RK), intent(in) :: exponent
        complex(RK), intent(in) :: Coef
        complex(RK), intent(in) :: Point
        complex(RK)             :: logEggBox
        logEggBox = exponent * log( constant + cos(coef*point) )
    end function getLogEggBoxSD_CK

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLogEggBoxMD_RK(nd,constant,exponent,Coef,Point) result(logEggBox)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogEggBoxMD_RK
#endif
        ! This simple function returns the log (in natural Neper base) of the EggBox function for nd > 1.
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

    pure function getLogEggBoxMD_CK(nd,constant,exponent,Coef,Point) result(logEggBox)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogEggBoxMD_CK
#endif
        ! This simple function returns the log (in natural Neper base) of the EggBox function for nd > 1.
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

    ! calculates the log Gamma function for a whole integer input. This is basically logGamma(positiveInteger+1)
    ! This function is mostly useful for large input integers
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

    ! calculates the Gamma function for a whole integer input. This is basically logGamma(intNum+1)
    pure function getFactorial(intNum)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getFactorial
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK), intent(in) :: intNum
        integer(IK)             :: i
        real(RK)                :: getFactorial
        getFactorial = 1._RK
        do i = 2,intNum
            getFactorial = getFactorial * i
        end do
    end function getFactorial

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! calculates the Gamma function for a half integer input.
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

    ! calculates the coefficient for the calculation of an nd-dimensional ellipsoid.
    ! see for example: http://math.stackexchange.com/questions/606184/volume-of-n-dimensional-ellipsoid
    ! see for example: https://en.wikipedia.org/wiki/Volume_of_an_n-ball
    ! see for example: https://en.wikipedia.org/wiki/Particular_values_of_the_Gamma_function
    ! getEllVolCoef = PI^(nd/2) / gamma(nd/2+1) where n is just a positive integer
    pure function getEllVolCoef(nd)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getEllVolCoef
#endif
        use Constants_mod, only: IK, RK, PI
        implicit none
        integer(IK), intent(in) :: nd
        integer(IK)             :: i,k
        real(RK)                :: getEllVolCoef
        if (mod(nd,2_IK)==0_IK) then  ! nd is even
            getEllVolCoef = PI
            do i = 2_IK, nd / 2_IK
                getEllVolCoef = getEllVolCoef * PI / i    ! nd = 2k ; getEllVolCoef = PI^(k) / Factorial(k)
            end do
        else    ! nd is an odd integer

            ! nd = 2k-1 ; gamma(nd/2 + 1) = gamma(k + 1/2) ; gamma(k+1/2) = sqrt(PI) * (2k)! / (4^k * k!)
            k = (nd + 1_IK) / 2_IK

            ! This is to avoid an extra unnecessary division of getEllVolCoef by PI
            getEllVolCoef = 4._RK / (k + 1_IK)

            do i = k+2_IK, 2_IK*k
                ! getEllVolCoef = PI^(k-1/2) / gamma(k+1/2) = PI^(k+1) * 4^k * k! / (2k)!
                getEllVolCoef = getEllVolCoef * PI * 4._RK / i
            end do

        end if
    end function getEllVolCoef

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns the log of an nd-dimensional unit ball.
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

    ! computes the volume of an nd-dimensional hyperellipsoid
    ! see for example: https://math.stackexchange.com/questions/2854930/volume-of-an-n-dimensional-ellipsoid
    pure function getLogVolEllipsoid(nd,logSqrtDetInvCovMat) result(logVolEllipsoid)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLogVolEllipsoid
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK) , intent(in)    :: nd
        real(RK)    , intent(in)    :: logSqrtDetInvCovMat
        real(RK)                    :: logVolEllipsoid
        logVolEllipsoid = getLogVolUnitBall(nd) + logSqrtDetInvCovMat
    end function getLogVolEllipsoid

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    ! This function has not been tested yet.
    ! There is really no reason to use this function for HPC solutions, as it can be likely done more efficiently.
    pure function getEllipsoidVolume(nEllipsoid,ndim,sqrtDeterminant)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getEllipsoidVolume
#endif
        use Constants_mod, only: IK, RK
        implicit none
        integer(IK)             :: i
        integer(IK), intent(in) :: nEllipsoid                   ! number of ellipsoids
        integer(IK), intent(in) :: ndim(nEllipsoid)             ! dimensions of the ellipsoids
        integer(IK), intent(in) :: sqrtDeterminant(nEllipsoid)  ! sqrtDeterminant of the covariance matrices of the ellipsoids
        real(RK)                :: getEllipsoidVolume(nEllipsoid)
        do i=1,nEllipsoid
            getEllipsoidVolume(i) = getEllVolCoef(ndim(i)) / sqrtDeterminant(i)
        end do 
    end function getEllipsoidVolume

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure function getLowerGamma(exponent,logGammaExponent,upperLim,tolerance) result(lowerGamma)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getLowerGamma
#endif
        ! returns the lower incomplete gamma function for the specified exponent and upper limit.
        ! tolerance represents the relative accuracy.
        ! NOTE: do not set tolerance to a number larger than 1, or else the code will crash spectacularly.
        ! logGammaExponent = log_gamma(exponent)
        ! NOTE: on input, both exponent and upperLim must be positive, otherwise, a -infinity will be returned to signal error.
        ! this algorithm borrows from the gammp implementation of the Numerical Recipes by Press et al 1992.
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

    pure function getUpperGamma(exponent,logGammaExponent,lowerLim,tolerance) result(upperGamma)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getUpperGamma
#endif
        ! returns the lower incomplete gamma function for the specified exponent and upper limit.
        ! tolerance represents the relative accuracy.
        ! NOTE: do not set tolerance to a number larger than 1, or else the code will crash spectacularly.
        ! NOTE: on input, both exponent and lowerLim must be positive, otherwise, a -infinity will be returned to signal error.
        ! this algorithm borrows from the gammq implementation of the Numerical Recipes by Press et al 1992.
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

    pure function getGammaSeries(exponent,logGammaExponent,upperLim,tolerance) result(gammaSeries)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaSeries
#endif
        ! Returns the lower incomplete gamma function P(exponent,upperLim) evaluated by its series representation as gamser.
        ! logGammaExponent is log(gamma(exponent)) on input.
        ! tolerance represents the relative accuracy.
        ! NOTE: do not set tolerance to a number larger than 1, or else the code will crash spectacularly.
        ! NOTE: if the algorithm fails to converge, gammaSeries will be set to negative infinity on output to signal error.
        ! this algorithm borrows from the gser implementation of the Numerical Recipes by Press et al 1992.
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

    pure function getGammaContFrac(exponent,logGammaExponent,lowerLim,tolerance) result(gammaContFrac)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getGammaContFrac
#endif
        ! Returns the incomplete gamma function Q(exponent,lowerLim) evaluated by its continued fraction representation as gammcf.
        ! logGammaExponent is log(gamma(exponent)) on input.
        ! ITMAX is the maximum allowed number of iterations.
        ! EPS represents the relative accuracy;
        ! FPMIN is a number near the smallest representable floating-point number.
        ! NOTE: if the algorithm fails to converge, gammaSeries will be set to negative infinity on output to signal error.
        ! NOTE: do not set tolerance to a number larger than 1, or else the code will crash spectacularly.
        ! this algorithm borrows from the gcf implementation of the Numerical Recipes by Press et al 1992.
        use Constants_mod, only: IK, RK, HUGE_RK
        implicit none
        real(RK), intent(in)            :: exponent, logGammaExponent, lowerLim
        real(RK), optional, intent(in)  :: tolerance
        real(RK)                        :: gammaContFrac
        integer(IK), parameter          :: ITMAX = 100
        real(RK), parameter             :: EPS_DEFAULT = epsilon(lowerLim), FPMIN_DEFAULT = tiny(lowerLim) / EPS_DEFAULT
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