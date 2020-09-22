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

module Integration_mod

    use Constants_mod, only: IK, RK

    implicit none

    real(RK), parameter :: ONE_THIRD = 1._RK / 3._RK

    character(len=117)  :: ErrorMessage(3) = &
    [ "@Integration_mod@doQuadRombClosed(): Too many steps in doQuadRombClosed().                                           " &
    , "@Integration_mod@doQuadRombClosed()@doPolInterp(): Input lowerLim, upperLim to doQuadRombClosed() are likely equal.  " &
    , "@Integration_mod@doQuadRombOpen()@doPolInterp(): Input lowerLim, upperLim to doQuadRombOpen() are likely equal.      " &
    ]

    abstract interface

        function integrand_proc(x) result(integrand)
            use Constants_mod, only: RK
            implicit none
            real(RK), intent(in)  :: x
            real(RK)              :: integrand
        end function integrand_proc

        subroutine integrator_proc(getFunc,lowerLim,upperLim,integral,refinementStage,numFuncEval)
            use Constants_mod, only: IK, RK
            import :: integrand_proc
            implicit none
            integer(IK) , intent(in)        :: refinementStage
            real(RK)    , intent(in)        :: lowerLim,upperLim
            integer(IK) , intent(out)       :: numFuncEval
            real(RK)    , intent(inout)     :: integral
            procedure(integrand_proc)       :: getFunc
        end subroutine integrator_proc

    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns the integral of function getFunc in the closed range [lowerLim,upperLim] using Adaptive Romberg extrapolation method.
    recursive subroutine doQuadRombClosed   ( getFunc &
                                            , lowerLim &
                                            , upperLim &
                                            , maxRelativeError &
                                            , nRefinement &
                                            , integral &
                                            , relativeError &
                                            , numFuncEval &
                                            , ierr &
                                            )
        !checked by Joshua Osborne on 5/28/2020 at 9:06pm
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: doQuadRombClosed
#endif
        use, intrinsic :: iso_fortran_env, only: DPI => int64
        implicit none
        integer(IK) , intent(in)    :: nRefinement
        real(RK)    , intent(in)    :: lowerLim,upperLim,maxRelativeError
        real(RK)    , intent(out)   :: integral, relativeError
        integer(IK) , intent(out)   :: ierr, numFuncEval
        integer(IK)                 :: refinementStage,km,numFuncEvalNew
        integer(IK), parameter      :: NSTEP = 31_IK
        real(RK)                    :: h(NSTEP+1),s(NSTEP+1)
        procedure(integrand_proc)   :: getFunc
        ierr = 0_IK
        km = nRefinement - 1_IK
        h(1) = 1._RK
        numFuncEval = 0_IK
        do refinementStage = 1, NSTEP
            call doQuadTrap(getFunc,lowerLim,upperLim,s(refinementStage),refinementStage,numFuncEvalNew)
            numFuncEval = numFuncEval + numFuncEvalNew
            if (refinementStage>=nRefinement) then
                call doPolInterp(h(refinementStage-km),s(refinementStage-km),nRefinement,0._RK,integral,relativeError,ierr)
                if ( abs(relativeError)<=maxRelativeError*abs(integral) .or. ierr/=0_IK ) return
            end if
            s(refinementStage+1)=s(refinementStage)
            h(refinementStage+1)=0.25_RK*h(refinementStage)
        end do
        ierr = 1_IK    ! too many steps in doQuadRombClosed()
  end subroutine doQuadRombClosed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! integrate getFunc() using Romberg adaptive extrapolation method on an Open interval (lowerLim,upperLim).
    recursive subroutine doQuadRombOpen ( getFunc &
                                        , integrate &
                                        , lowerLim &
                                        , upperLim &
                                        , maxRelativeError &
                                        , nRefinement &
                                        , integral &
                                        , relativeError &
                                        , numFuncEval &
                                        , ierr &
                                        )
        !checked by Joshua Osborne on 5/28/2020 at 9:03pm
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: doQuadRombOpen
#endif
        use, intrinsic :: iso_fortran_env, only: DPI => int64
        implicit none
        integer(IK) , intent(in)    :: nRefinement
        real(RK)    , intent(in)    :: lowerLim,upperLim,maxRelativeError
        real(RK)    , intent(out)   :: integral, relativeError
        integer(IK) , intent(out)   :: numFuncEval, ierr
        integer(IK)                 :: refinementStage,km,numFuncEvalNew
        integer(IK) , parameter     :: NSTEP = 20_IK
        real(RK)    , parameter     :: ONE_OVER_NINE = 1._RK / 9._RK
        real(RK)                    :: h(NSTEP+1),s(NSTEP+1)
        procedure(integrand_proc)   :: getFunc
        procedure(integrator_proc)  :: integrate
        ierr = 0_IK
        km = nRefinement-1
        h(1) = 1.0_RK
        numFuncEval = 0_IK
        do refinementStage = 1, NSTEP
            call integrate(getFunc,lowerLim,upperLim,s(refinementStage),refinementStage,numFuncEvalNew)
            numFuncEval = numFuncEval + numFuncEvalNew
            if (refinementStage>=nRefinement) then
                call doPolInterp(h(refinementStage-km),s(refinementStage-km),nRefinement,0.0_RK,integral,relativeError,ierr)
                if ( abs(relativeError)<=maxRelativeError*abs(integral) .or. ierr/=0_IK ) return
            end if
            s(refinementStage+1)=s(refinementStage)
            h(refinementStage+1)=h(refinementStage)*ONE_OVER_NINE
        end do
        ierr = 2_IK    ! too many steps in doQuadRombOpen()
  end subroutine doQuadRombOpen

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    recursive subroutine doPolInterp(xa,ya,ndata,x,y,dy,ierr)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: doPolInterp
#endif
        implicit none
        integer , intent(in)    :: ndata
        real(RK), intent(in)    :: x,xa(ndata),ya(ndata)
        real(RK), intent(out)   :: dy,y
        integer , intent(out)   :: ierr
        real(RK)                :: den,dif,dift,ho,hp,w,c(ndata),d(ndata)
        integer                 :: i,m,ns
        ierr = 0
        ns=1
        dif=abs(x-xa(1))
        do i=1,ndata
            dift=abs(x-xa(i))
            if (dift<dif) then
                ns=i
                dif=dift
            endif
            c(i)=ya(i)
            d(i)=ya(i)
        end do
        y=ya(ns)
        ns=ns-1
        do m=1,ndata-1
            do i=1,ndata-m
                ho=xa(i)-x
                hp=xa(i+m)-x
                w=c(i+1)-d(i)
                den=ho-hp
                if(den==0._RK) then
                    ierr = 3 ! failure in doPolInterp()
                    return
                end if
                den=w/den
                d(i)=hp*den
                c(i)=ho*den
            end do
            if (2*ns<ndata-m)then
                dy=c(ns+1)
            else
                dy=d(ns)
                ns=ns-1
            endif
            y=y+dy
        end do
    end subroutine doPolInterp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    recursive subroutine doQuadTrap(getFunc,lowerLim,upperLim,integral,refinementStage,numFuncEval)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: doQuadTrap
#endif
        implicit none
        integer(IK), intent(in)     :: refinementStage
        real(RK), intent(in)        :: lowerLim,upperLim
        real(RK), intent(inout)     :: integral
        integer(IK), intent(out)    :: numFuncEval
        integer(IK)                 :: iFuncEval
        real(RK)                    :: del,sum,tnm,x
        procedure(integrand_proc)   :: getFunc
        if (refinementStage==1) then
            numFuncEval = 2_IK
            integral = 0.5_RK*(upperLim-lowerLim)*(getFunc(lowerLim)+getFunc(upperLim))
        else
            numFuncEval = 2**(refinementStage-2)
            tnm = real(numFuncEval,kind=RK)
            del = (upperLim-lowerLim) / tnm
            x = lowerLim + 0.5_RK * del
            sum = 0._RK
            do iFuncEval = 1, numFuncEval
                sum=sum+getFunc(x)
                x=x+del
            end do
            integral=0.5_RK*(integral+(upperLim-lowerLim)*sum/tnm)
        endif
    end subroutine doQuadTrap

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    recursive subroutine midexp(getFunc,lowerLim,upperLim,integral,refinementStage,numFuncEval)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: midexp
#endif
        implicit none
        integer(IK) , intent(in)    :: refinementStage
        real(RK)    , intent(in)    :: lowerLim,upperLim
        integer(IK) , intent(out)   :: numFuncEval
        real(RK)    , intent(inout) :: integral
        procedure(integrand_proc)   :: getFunc
        real(RK)                    :: ddel,del,summ,x,lowerLimTrans,upperLimTrans
        real(RK)                    :: inverseThreeNumFuncEval
        integer(IK)                 :: iFuncEval
        upperLimTrans = exp(-lowerLim)
        lowerLimTrans = exp(-upperLim)
        if (refinementStage==1) then
            numFuncEval = 1_IK
            integral = (upperLimTrans-lowerLimTrans)*getTransFunc(0.5_RK*(lowerLimTrans+upperLimTrans))
        else
            numFuncEval = 3**(refinementStage-2)
            inverseThreeNumFuncEval = ONE_THIRD / numFuncEval
            del = (upperLimTrans-lowerLimTrans) * inverseThreeNumFuncEval
            ddel = del + del
            x = lowerLimTrans + 0.5_RK*del
            summ=0._RK
            do iFuncEval = 1,numFuncEval
                summ = summ + getTransFunc(x)
                x = x + ddel
                summ = summ + getTransFunc(x)
                x = x + del
            end do
            integral = ONE_THIRD * integral + (upperLimTrans-lowerLimTrans) * summ * inverseThreeNumFuncEval
            numFuncEval = 2_IK * numFuncEval
        end if
    contains
        function getTransFunc(x)
            implicit none
            real(RK), intent(in)    :: x
            real(RK)                :: getTransFunc
            getTransFunc = getFunc(-log(x)) / x
        end function getTransFunc
    end subroutine midexp

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! must be: lowerLim*upperLim > 0.0
    ! integrate a function on a semi-infinite interval
    subroutine midinf(getFunc,lowerLim,upperLim,integral,refinementStage,numFuncEval)
    !checked by Joshua Osborne on 5/28/2020 at 8:58pm
        implicit none
        real(RK)    , intent(in)    :: lowerLim,upperLim
        integer(IK) , intent(in)    :: refinementStage
        integer(IK) , intent(out)   :: numFuncEval
        real(RK)    , intent(inout) :: integral
        procedure(integrand_proc)   :: getFunc
        real(RK)                    :: lowerLimTrans, upperLimTrans, del, ddel, summ, x
        real(RK)                    :: inverseThreeNumFuncEval
        integer(IK)                 :: iFuncEval
        upperLimTrans = 1.0_RK / lowerLim
        lowerLimTrans = 1.0_RK / upperLim
        if (refinementStage == 1_IK) then
            numFuncEval = 1_IK
            integral = (upperLimTrans-lowerLimTrans) * getTransFunc(0.5_RK * (lowerLimTrans+upperLimTrans))
        else
            numFuncEval = 3**(refinementStage-2)
            inverseThreeNumFuncEval = ONE_THIRD / numFuncEval
            del = (upperLimTrans-lowerLimTrans) * inverseThreeNumFuncEval
            ddel = del + del
            x = lowerLimTrans + 0.5_RK * del
            summ = 0._RK
            do iFuncEval = 1, numFuncEval
                summ = summ + getTransFunc(x)
                x = x + ddel
                summ = summ + getTransFunc(x)
                x = x + del
            end do
            integral = ONE_THIRD * integral + (upperLimTrans-lowerLimTrans) * summ * inverseThreeNumFuncEval
            numFuncEval = 2_IK * numFuncEval
        end if
    contains
        function getTransFunc(x)
            implicit none
            real(RK), intent(in)    :: x
            real(RK)                :: getTransFunc
            getTransFunc = getFunc(1._RK/x) / x**2
        end function getTransFunc
    end subroutine midinf

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine midpnt(getFunc,lowerLim,upperLim,integral,refinementStage,numFuncEval)
    !checked by Joshua Osborne on 5/28/2020 at 8:55pm
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: midpnt
#endif
        implicit none
        integer(IK) , intent(in)    :: refinementStage
        real(RK)    , intent(in)    :: lowerLim, upperLim
        real(RK)    , intent(inout) :: integral
        integer(IK) , intent(out)   :: numFuncEval
        procedure(integrand_proc)   :: getFunc
        integer(IK)                 :: iFuncEval
        real(RK)                    :: ddel,del,summ,x
        real(RK)                    :: inverseThreeNumFuncEval
        if (refinementStage==1) then
            numFuncEval = 1_IK
            integral = (upperLim-lowerLim) * getFunc( 0.5_RK * (lowerLim+upperLim) )
        else
            numFuncEval = 3_IK**(refinementStage-2)
            inverseThreeNumFuncEval = ONE_THIRD / numFuncEval
            del = (upperLim-lowerLim) * inverseThreeNumFuncEval
            ddel = del+del
            x = lowerLim + 0.5_RK * del
            summ = 0._RK
            do iFuncEval = 1, numFuncEval
                summ = summ + getFunc(x)
                x = x + ddel
                summ = summ + getFunc(x)
                x = x + del
            end do
            integral = ONE_THIRD * integral + (upperLim-lowerLim) * summ * inverseThreeNumFuncEval
            numFuncEval = 2_IK * numFuncEval
        end if
    end subroutine midpnt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Integration_mod