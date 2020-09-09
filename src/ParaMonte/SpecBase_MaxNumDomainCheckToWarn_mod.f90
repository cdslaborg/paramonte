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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(MaxNumDomainCheckToWarnObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(MaxNumDomainCheckToWarn_type), intent(in)     :: MaxNumDomainCheckToWarnObj
        maxNumDomainCheckToWarn = MaxNumDomainCheckToWarnObj%null
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecBase_MaxNumDomainCheckToWarn_mod