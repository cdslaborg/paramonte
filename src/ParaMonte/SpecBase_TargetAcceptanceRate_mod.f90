!**********************************************************************************************************************************
!**********************************************************************************************************************************
!
!  ParaMonte: plain powerful parallel Monte Carlo library.
!
!  Copyright (C) 2012-present, The Computational Data Science Lab
!
!  This file is part of the ParaMonte library. 
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

module SpecBase_TargetAcceptanceRate_mod

    use Constants_mod, only: RK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecBase_TargetAcceptanceRate_mod"

    real(RK)                        :: TargetAcceptanceRate(2) ! namelist input

    type                            :: TargetAcceptanceRate_type
        logical                     :: scalingRequested
        real(RK)                    :: Val(2)
        real(RK)                    :: Def(2)
        real(RK)                    :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setTargetAcceptanceRate, checkForSanity, nullifyNameListVar
    end type TargetAcceptanceRate_type

    interface TargetAcceptanceRate_type
        module procedure            :: constructTargetAcceptanceRate
    end interface TargetAcceptanceRate_type

    private :: constructTargetAcceptanceRate, setTargetAcceptanceRate, nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructTargetAcceptanceRate(methodName) result(TargetAcceptanceRateObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructTargetAcceptanceRate
#endif
        use Constants_mod, only: NULL_RK
        use String_mod, only: num2str
        implicit none
        character(*), intent(in)            :: methodName
        type(TargetAcceptanceRate_type)     :: TargetAcceptanceRateObj
        TargetAcceptanceRateObj%scalingRequested = .true.
        TargetAcceptanceRateObj%Def = [ 0._RK, 1._RK ]
        TargetAcceptanceRateObj%null = NULL_RK
        TargetAcceptanceRateObj%desc = &
#include "SpecBase_TargetAcceptanceRate_desc.f90"
    end function constructTargetAcceptanceRate

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(TargetAcceptanceRateObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(TargetAcceptanceRate_type), intent(in)  :: TargetAcceptanceRateObj
        TargetAcceptanceRate = TargetAcceptanceRateObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setTargetAcceptanceRate(TargetAcceptanceRateObj,TargetAcceptanceRate)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setTargetAcceptanceRate
#endif
        use Constants_mod, only: RK
        implicit none
        class(TargetAcceptanceRate_type), intent(inout)     :: TargetAcceptanceRateObj
        real(RK), intent(in)                                :: TargetAcceptanceRate(2)
        logical                                             :: lowerLimitSet, upperLimitSet
        TargetAcceptanceRateObj%Val = TargetAcceptanceRate
        lowerLimitSet = TargetAcceptanceRateObj%Val(1)/=TargetAcceptanceRateObj%null
        upperLimitSet = TargetAcceptanceRateObj%Val(2)/=TargetAcceptanceRateObj%null
        if      ( lowerLimitSet .and. (.not. upperLimitSet) ) then
            TargetAcceptanceRateObj%Val(2) = TargetAcceptanceRateObj%Val(1)
        elseif  ( upperLimitSet .and. (.not. lowerLimitSet) ) then
            TargetAcceptanceRateObj%Val(1) = TargetAcceptanceRateObj%Val(2)
        elseif  ( .not. (lowerLimitSet .or.  upperLimitSet) ) then
            TargetAcceptanceRateObj%Val = TargetAcceptanceRateObj%Def
            TargetAcceptanceRateObj%scalingRequested = .false.
        elseif ( all(TargetAcceptanceRateObj%Val==TargetAcceptanceRateObj%Def) ) then
            TargetAcceptanceRateObj%scalingRequested = .false.
        end if
    end subroutine setTargetAcceptanceRate

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine checkForSanity(TargetAcceptanceRateObj,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: RK
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(TargetAcceptanceRate_type), intent(in) :: TargetAcceptanceRateObj
        type(Err_type), intent(inout)       :: Err
        character(*), parameter             :: PROCEDURE_NAME = "@checkForSanity()"
        if (.not. TargetAcceptanceRateObj%scalingRequested) return
        if ( any(TargetAcceptanceRateObj%val<0._RK) .or. any(TargetAcceptanceRateObj%val>1._RK) ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The target acceptance ratio limits targetAcceptanceRate [" // &
                        num2str(TargetAcceptanceRateObj%val) // "," // num2str(TargetAcceptanceRateObj%val) // &
                        "] cannot be less than 0 or larger than 1.\n\n"
        end if
        if ( all(TargetAcceptanceRateObj%val==0._RK) .or. all(TargetAcceptanceRateObj%val==1._RK) ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The target acceptance ratio limits targetAcceptanceRate [" // &
                        num2str(TargetAcceptanceRateObj%val) // "," // num2str(TargetAcceptanceRateObj%val) // &
                        "] cannot be both 0 or both 1.\n\n"
        end if
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecBase_TargetAcceptanceRate_mod