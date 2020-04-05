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

module SpecBase_ProgressReportPeriod_mod

    use Constants_mod, only: IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecBase_ProgressReportPeriod_mod"

    integer(IK)                     :: progressReportPeriod ! namelist input

    type                            :: ProgressReportPeriod_type
        integer(IK)                 :: val
        integer(IK)                 :: def
        integer(IK)                 :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setProgressReportPeriod, checkForSanity, nullifyNameListVar
    end type ProgressReportPeriod_type

    interface ProgressReportPeriod_type
        module procedure                :: constructProgressReportPeriod
    end interface ProgressReportPeriod_type

    private :: constructProgressReportPeriod, setProgressReportPeriod

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructProgressReportPeriod() result(ProgressReportPeriodObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructProgressReportPeriod
#endif
        use Constants_mod, only: IK, NULL_IK
        use String_mod, only: num2str
        implicit none
        type(ProgressReportPeriod_type) :: ProgressReportPeriodObj
        ProgressReportPeriodObj%def = 1000_IK
        ProgressReportPeriodObj%null = NULL_IK
        ProgressReportPeriodObj%desc = &
        "Every progressReportPeriod calls to the objective function, the sampling progress will be reported to the log file. &
        &Note that progressReportPeriod must be a positive integer. &
        &The default value is " // num2str(ProgressReportPeriodObj%def) // "."
    end function constructProgressReportPeriod

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(ProgressReportPeriodObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(ProgressReportPeriod_type), intent(inout)  :: ProgressReportPeriodObj
        progressReportPeriod = ProgressReportPeriodObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setProgressReportPeriod(ProgressReportPeriodObj,progressReportPeriod)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setProgressReportPeriod
#endif
        use Constants_mod, only: IK
        implicit none
        class(ProgressReportPeriod_type), intent(inout) :: ProgressReportPeriodObj
        integer(IK), intent(in)                         :: progressReportPeriod
        ProgressReportPeriodObj%val = progressReportPeriod
        if (ProgressReportPeriodObj%val==ProgressReportPeriodObj%null) then
            ProgressReportPeriodObj%val = ProgressReportPeriodObj%def
        end if
    end subroutine setProgressReportPeriod

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine checkForSanity(ProgressReportPeriod,Err,methodName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(ProgressReportPeriod_type), intent(in)    :: ProgressReportPeriod
        character(*), intent(in)                        :: methodName
        type(Err_type), intent(inout)                   :: Err
        character(*), parameter                         :: PROCEDURE_NAME = "@checkForSanity()"
        if ( ProgressReportPeriod%val<1 ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input value for variable progressReportPeriod must be a positive integer value. If you are not sure &
                        &about the appropriate value for this variable, simply drop it from the input. " // &
                        methodName // " will automatically assign an appropriate value to it.\n\n"
        end if
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecBase_ProgressReportPeriod_mod