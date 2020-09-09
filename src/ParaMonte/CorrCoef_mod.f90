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

module CorrCoef_mod

    use Constants_mod, only: RK, IK
    use Err_mod, only: Err_type
    implicit none

    character(*), parameter :: MODULE_NAME = "@CorrCoef_mod"

    type :: CorrCoefSpearman_type
        integer(IK)             :: ndata
        real(RK), allocatable   :: Data1(:), Data2(:)
        real(RK)                :: rho, rhoProb, dStarStar, dStarStarSignif, dStarStarProb
        type(Err_type)          :: Err
    contains
        procedure, nopass :: get => getCorrCoefSpearman
    end type CorrCoefSpearman_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Given two data arrays, data1(1:ndata) and data2(1:ndata), this routine returns:
    ! their sum-squared difference of ranks as dStarStar,
    ! the number of standard deviations by which dStarStar deviates from its null-hypothesis expected value as dStarStarSignif,
    ! the two-sided significance level of this deviation as dStarStarProb,
    ! the Spearman rank correlation as rho,
    ! the two-sided significance level of its deviation from zero as rhoProb.
    ! The external routines crank and sortAscending are used.
    ! A small value of either dStarStarProb or rhoProb indicates a significant correlation.
    subroutine getCorrCoefSpearman(ndata,Data1,Data2,rho,rhoProb,dStarStar,dStarStarSignif,dStarStarProb,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getCorrCoefSpearman
#endif
        use Constants_mod, only: IK, RK, SQRT2, SPR
        use Statistics_mod, only: getBetaCDF
        use Sort_mod, only: sortAscending

        implicit none

        integer(IK)     , intent(in)    :: ndata
        real(RK)        , intent(in)    :: Data1(ndata), Data2(ndata)
        real(RK)        , intent(out)   :: rho, rhoProb
        real(RK)        , intent(out)   :: dStarStar, dStarStarSignif, dStarStarProb
        type(Err_type)  , intent(out)   :: Err

        character(*)    , parameter     :: PROCEDURE_NAME = MODULE_NAME//"@getCorrCoefSpearman"
        real(RK)                        :: aved,df,en,en3n,fac,sf,sg,t,vard
        real(RK)                        :: WorkSpace1(ndata), WorkSpace2(ndata)

        WorkSpace1 = Data1
        WorkSpace2 = Data2
        call sortAscending(ndata,WorkSpace1,WorkSpace2,Err)
        if (Err%occurred) then
            Err%msg = PROCEDURE_NAME//Err%msg
            return
        end if
        call crank(ndata,WorkSpace1,sf)

        call sortAscending(ndata,WorkSpace2,WorkSpace1,Err)
        if (Err%occurred) then
            Err%msg = PROCEDURE_NAME//Err%msg
            return
        end if
        call crank(ndata,WorkSpace2,sg)

        WorkSpace1 = WorkSpace1 - WorkSpace2
        dStarStar = dot_product(WorkSpace1,WorkSpace1)
        en = ndata
        en3n = en**3-en
        aved = en3n/6.0_RK-(sf+sg)/12.0_RK
        fac = (1.0_RK-sf/en3n)*(1.0_RK-sg/en3n)
        vard = ((en-1.0_RK)*en**2*(en+1.0_RK)**2/36.0_RK)*fac
        dStarStarSignif = (dStarStar-aved)/sqrt(vard)
        dStarStarProb = erfc( real( abs(dStarStarSignif)/SQRT2 , kind=SPR ) )
        rho = (1.0_RK-(6.0_RK/en3n)*(dStarStar+(sf+sg)/12.0_RK))/sqrt(fac)
        fac = (1.0_RK+rho)*(1.0_RK-rho)
        if (fac > 0.0) then
            t = rho*sqrt((en-2.0_RK)/fac)
            df = en-2.0_RK
            rhoProb = getBetaCDF(0.5_RK*df,0.5_RK,df/(df+t**2))
        else
            rhoProb = 0.0
        end if

    contains

        subroutine crank(ndata,Array,s)
            use Constants_mod, only: IK, RK
            use Misc_mod, only : arth, copyArray
            implicit none
            integer(IK), intent(in) :: ndata
            real(RK), intent(out) :: s
            real(RK), intent(inout) :: Array(ndata)
            integer(IK) :: i,ndum,nties
            integer(IK) :: Tstart(ndata), Tend(ndata), Tie(ndata), Idx(ndata)
            Idx = arth(1,1,ndata)
            Tie = merge(1,0,Array == eoshift(Array,-1))
            Tie(1) = 0
            Array = Idx(:)
            if (all(Tie == 0)) then
                s = 0.0
                return
            end if
            call copyArray( pack(Idx,Tie<eoshift(Tie,1)), Tstart, nties, ndum )
            Tend(1:nties) = pack(Idx(:),Tie(:)>eoshift(Tie(:),1))
            do i = 1,nties
                Array(Tstart(i):Tend(i)) = (Tstart(i)+Tend(i))/2.0_RK
            end do
            Tend(1:nties) = Tend(1:nties)-Tstart(1:nties)+1
            s = sum(Tend(1:nties)**3-Tend(1:nties))
        end subroutine crank

    end subroutine getCorrCoefSpearman

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module CorrCoef_mod