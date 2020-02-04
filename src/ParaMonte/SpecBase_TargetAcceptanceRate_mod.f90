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

module SpecBase_TargetAcceptanceRate_mod

    use Constants_mod, only: RK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecBase_TargetAcceptanceRate_mod"

    real(RK)                        :: targetAcceptanceRate ! namelist input

    type                            :: TargetAcceptanceRate_type
        logical                     :: scalingRequested
        real(RK)                    :: val
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
        TargetAcceptanceRateObj%null = NULL_RK
        TargetAcceptanceRateObj%desc = &
        "targetAcceptanceRate sets an optimal target for the ratio of the number of accepted objective function calls to the &
        &total number of function calls by " // methodName // ". By default, it is a real number between 0 and 1. &
        &If provided by the user, " // methodName // " will attempt (but not guarantee) to bring the average acceptance ratio &
        &of the sampler as close to the user-provided target ratio as possible. The success of " // methodName // &
        " in keeping the average acceptance ratio close to the requested target value depends heavily on:\n&
        &    1) the value of adaptiveUpdatePeriod; the larger, the easier.\n&
        &    2) the value of adaptiveUpdateCount; the larger, the easier.\n&
        &Note that the acceptance ratio adjustments will only occur every adaptiveUpdatePeriod sampling steps for a total &
        &number of adaptiveUpdateCount. &
        &There is no default value for targetAcceptanceRate, as the acceptance ratio is not directly adjusted during sampling."
    end function constructTargetAcceptanceRate

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(TargetAcceptanceRateObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(TargetAcceptanceRate_type), intent(in)  :: TargetAcceptanceRateObj
        targetAcceptanceRate = TargetAcceptanceRateObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setTargetAcceptanceRate(TargetAcceptanceRateObj,targetAcceptanceRate)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setTargetAcceptanceRate
#endif
        use Constants_mod, only: RK
        implicit none
        class(TargetAcceptanceRate_type), intent(inout)     :: TargetAcceptanceRateObj
        real(RK), intent(in)                                :: targetAcceptanceRate
        TargetAcceptanceRateObj%val = targetAcceptanceRate
        if (TargetAcceptanceRateObj%val==TargetAcceptanceRateObj%null) TargetAcceptanceRateObj%scalingRequested = .false.
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
        if ( TargetAcceptanceRateObj%val<=0._RK ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The target acceptance ratio targetAcceptanceRate (" // num2str(TargetAcceptanceRateObj%val) // &
                        ") cannot be less than or equal to 0.\n\n"
        end if
        if ( TargetAcceptanceRateObj%val>=1._RK ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The target acceptance ratio targetAcceptanceRate (" // num2str(TargetAcceptanceRateObj%val) // &
                        ") cannot be larger than or equal to 1.\n\n"
        end if
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecBase_TargetAcceptanceRate_mod