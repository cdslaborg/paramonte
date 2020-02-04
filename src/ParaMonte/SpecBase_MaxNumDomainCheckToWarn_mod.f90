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

module SpecBase_MaxNumDomainCheckToWarn_mod

    use Constants_mod, only: IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecBase_MaxNumDomainCheckToWarn_mod"

    integer(IK)                     :: maxNumDomainCheckToWarn ! namelist input

    type                            :: MaxNumDomainCheckToWarn_type
        integer(IK)                 :: val
        integer(IK)                 :: def
        integer(IK)                 :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setMaxNumDomainCheckToWarn, checkForSanity, nullifyNameListVar
    end type MaxNumDomainCheckToWarn_type

    interface MaxNumDomainCheckToWarn_type
        module procedure                :: constructMaxNumDomainCheckToWarn
    end interface MaxNumDomainCheckToWarn_type

    private :: constructMaxNumDomainCheckToWarn, setMaxNumDomainCheckToWarn

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructMaxNumDomainCheckToWarn() result(MaxNumDomainCheckToWarnObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructMaxNumDomainCheckToWarn
#endif
        use Constants_mod, only: IK, NULL_IK
        use Decoration_mod, only: TAB
        use String_mod, only: num2str
        implicit none
        type(MaxNumDomainCheckToWarn_type) :: MaxNumDomainCheckToWarnObj
        MaxNumDomainCheckToWarnObj%def = 1000_IK
        MaxNumDomainCheckToWarnObj%null = NULL_IK
        MaxNumDomainCheckToWarnObj%desc = &
        "maxNumDomainCheckToWarn is an integer number beyond which the user will be warned about the newly-proposed points &
        &being excessively proposed outside the domain of the objective function. For every maxNumDomainCheckToWarn &
        &consecutively-proposed new points that fall outside the domain of the objective function, the user will be warned until &
        &maxNumDomainCheckToWarn = maxNumDomainCheckToStop, in which case the sampler returns a fatal error and the program stops &
        &globally. The counter for this warning message is reset after a proposal sample from within the domain of the &
        &objective function is obtained. The default value is " // num2str(MaxNumDomainCheckToWarnObj%def) // "."
    end function constructMaxNumDomainCheckToWarn

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(MaxNumDomainCheckToWarnObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(MaxNumDomainCheckToWarn_type), intent(in)     :: MaxNumDomainCheckToWarnObj
        maxNumDomainCheckToWarn = MaxNumDomainCheckToWarnObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setMaxNumDomainCheckToWarn(MaxNumDomainCheckToWarnObj,maxNumDomainCheckToWarn)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setMaxNumDomainCheckToWarn
#endif
        use Constants_mod, only: IK
        implicit none
        class(MaxNumDomainCheckToWarn_type), intent(inout)  :: MaxNumDomainCheckToWarnObj
        integer(IK), intent(in)                             :: maxNumDomainCheckToWarn
        MaxNumDomainCheckToWarnObj%val = maxNumDomainCheckToWarn
        if (MaxNumDomainCheckToWarnObj%val==MaxNumDomainCheckToWarnObj%null) then
            MaxNumDomainCheckToWarnObj%val = MaxNumDomainCheckToWarnObj%def
        end if
    end subroutine setMaxNumDomainCheckToWarn

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine checkForSanity(MaxNumDomainCheckToWarn,Err,methodName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(MaxNumDomainCheckToWarn_type), intent(in) :: MaxNumDomainCheckToWarn
        character(*), intent(in)                        :: methodName
        type(Err_type), intent(inout)                   :: Err
        character(*), parameter                         :: PROCEDURE_NAME = "@checkForSanity()"
        if ( MaxNumDomainCheckToWarn%val<1 ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input value for variable maxNumDomainCheckToWarn must be a positive integer. If you are not sure &
                        &about the appropriate value for this variable, simply drop it from the input. " // &
                        methodName // " will automatically assign an appropriate value to it.\n\n"
        end if
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecBase_MaxNumDomainCheckToWarn_mod