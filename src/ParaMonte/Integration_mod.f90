!**********************************************************************************************************************************
!**********************************************************************************************************************************
!
!  ParaMonte: plain powerful parallel Monte Carlo library.
!
!  Copyright (C) 2012-present, The Computational Data Science Lab
!
!  This file is part of ParaMonte library. 
!
!  ParaMonte is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, version 3 of the License.
!
!  ParaMonte is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with ParaMonte.  If not, see <https://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!**********************************************************************************************************************************

module Integration_mod
  
    use Constants_mod, only: IK, RK
  
    implicit none

    character(len=117) :: ErrorMessage(3) = &
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
            real(RK)    , intent(out)       :: integral
            integer(IK) , intent(out)       :: numFuncEval
            procedure(integrand_proc)       :: getFunc
        end subroutine integrator_proc
  
    end interface

!***********************************************************************************************************************************
!***********************************************************************************************************************************
  
contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

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
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: doQuadRombClosed
#endif
        implicit none
        integer(IK), intent(in)     :: nRefinement
        real(RK), intent(in)        :: lowerLim,upperLim,maxRelativeError
        real(RK), intent(out)       :: integral, relativeError
        integer(IK), intent(out)    :: ierr, numFuncEval
        integer(IK)                 :: refinementStage,km,numFuncEvalNew
        integer(IK), parameter      :: NSTEP = 20_IK
        real(RK)                    :: h(NSTEP+1),s(NSTEP+1)
        procedure(integrand_proc)   :: getFunc
        ierr = 0_IK
        km = nRefinement - 1_IK
        h(1) = 1._RK
        numFuncEval = 0_IK
        do refinementStage=1,NSTEP
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

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    ! integrate getFunc() using Romberg's adaptive extrapolation method on an Open interval (lowerLim,upperLim).
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
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: doQuadRombOpen
#endif

        implicit none
        integer(IK), intent(in)     :: nRefinement
        real(RK), intent(in)        :: lowerLim,upperLim,maxRelativeError
        real(RK), intent(out)       :: integral, relativeError
        integer(IK), intent(out)    :: numFuncEval, ierr
        integer(IK)                 :: refinementStage,km,numFuncEvalNew
        integer, parameter          :: NSTEP = 20_IK
        real(RK), parameter         :: ONE_OVER_NINE = 1._RK / 9._RK
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

!***********************************************************************************************************************************
!***********************************************************************************************************************************

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

!***********************************************************************************************************************************
!***********************************************************************************************************************************
  
    recursive subroutine doQuadTrap(getFunc,lowerLim,upperLim,integral,refinementStage,numFuncEval)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: doQuadTrap
#endif
        implicit none
        integer(IK), intent(in)     :: refinementStage
        real(RK), intent(in)        :: lowerLim,upperLim
        real(RK), intent(out)       :: integral
        integer(IK), intent(out)    :: numFuncEval
        integer(IK)                 :: iFuncEval
        real(RK)                    :: del,sum,tnm,x
        procedure(integrand_proc)   :: getFunc
        if (refinementStage==1) then
            numFuncEval = 2_IK
            integral=0.5_RK*(upperLim-lowerLim)*(getFunc(lowerLim)+getFunc(upperLim))
        else
            numFuncEval=2**(refinementStage-2)
            tnm=real(numFuncEval,kind=RK)
            del=(upperLim-lowerLim)/tnm
            x=lowerLim+0.5*del
            sum=0._RK
            do iFuncEval = 1, numFuncEval
                sum=sum+getFunc(x)
                x=x+del
            end do
            integral=0.5_RK*(integral+(upperLim-lowerLim)*sum/tnm)
        endif
    end subroutine doQuadTrap

!***********************************************************************************************************************************
!***********************************************************************************************************************************
  
    recursive subroutine midexp(getFunc,lowerLim,upperLim,integral,refinementStage,numFuncEval)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: midexp
#endif
        implicit none
        integer(IK) , intent(in)    :: refinementStage
        real(RK)    , intent(in)    :: lowerLim,upperLim
        real(RK)    , intent(out)   :: integral
        integer(IK) , intent(out)   :: numFuncEval
        procedure(integrand_proc)   :: getFunc
        real(RK)                    :: ddel,del,summ,tnm,x,lowerLimTrans,upperLimTrans
        integer(IK)                 :: iFuncEval
        upperLimTrans = exp(-lowerLim)
        lowerLimTrans = exp(-upperLim)
        if (refinementStage==1) then
            numFuncEval = 1_IK
            integral = (upperLimTrans-lowerLimTrans)*getTransFunc(0.5_RK*(lowerLimTrans+upperLimTrans))
        else
            numFuncEval = 3**(refinementStage-2)
            tnm = numFuncEval
            del = (upperLimTrans-lowerLimTrans)/(3._RK*tnm)
            ddel = del + del
            x = lowerLimTrans + 0.5_RK*del
            summ=0._RK
            do iFuncEval = 1,numFuncEval
                summ = summ + getTransFunc(x)
                x = x + ddel
                summ = summ + getTransFunc(x)
                x = x + del
            end do
            integral = ( integral + (upperLimTrans-lowerLimTrans)*summ/tnm ) / 3._RK
            numFuncEval = 2_IK * numFuncEval
        end if
    contains
        function getTransFunc(x)
            real(RK), intent(in) :: x
            real(RK)             :: getTransFunc
            getTransFunc = getFunc(-log(x)) / x
        end function getTransFunc
    end subroutine midexp

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine midpnt(getFunc,lowerLimit,upperLimit,integral,refinementStage,numFuncEval)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: midpnt
#endif
        implicit none
        integer(IK) , intent(in)    :: refinementStage
        real(RK)    , intent(in)    :: lowerLimit,upperLimit
        real(RK)    , intent(out)   :: integral
        integer(IK) , intent(out)   :: numFuncEval
        procedure(integrand_proc)   :: getFunc
        integer                     :: iFuncEval
        real(RK) ddel,del,summ,tnm,x
        if (refinementStage==1) then
            numFuncEval = 1_IK
            integral = (upperLimit-lowerLimit) * getFunc( 0.5_RK * (lowerLimit+upperLimit) )
        else
            numFuncEval = 3**(refinementStage-2)
            tnm = numFuncEval
            del = (upperLimit-lowerLimit) / (3.0_RK*tnm)
            ddel = del+del
            x = lowerLimit + 0.5_RK*del
            summ = 0._RK
            do iFuncEval = 1,numFuncEval
                summ = summ + getFunc(x)
                x = x + ddel
                summ = summ + getFunc(x)
                x = x + del
            end do
            integral = ( integral + (upperLimit-lowerLimit)*summ/tnm ) / 3._RK
            numFuncEval = 2_IK * numFuncEval
        end if
    end subroutine midpnt

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module Integration_mod