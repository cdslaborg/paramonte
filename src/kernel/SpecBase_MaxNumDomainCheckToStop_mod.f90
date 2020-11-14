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

module SpecBase_MaxNumDomainCheckToStop_mod

    use Constants_mod, only: IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecBase_MaxNumDomainCheckToStop_mod"

    integer(IK)                     :: maxNumDomainCheckToStop ! namelist input

    type                            :: MaxNumDomainCheckToStop_type
        integer(IK)                 :: val
        integer(IK)                 :: def
        integer(IK)                 :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setMaxNumDomainCheckToStop, checkForSanity, nullifyNameListVar
    end type MaxNumDomainCheckToStop_type

    interface MaxNumDomainCheckToStop_type
        module procedure                :: constructMaxNumDomainCheckToStop
    end interface MaxNumDomainCheckToStop_type

    private :: constructMaxNumDomainCheckToStop, setMaxNumDomainCheckToStop

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructMaxNumDomainCheckToStop() result(MaxNumDomainCheckToStopObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructMaxNumDomainCheckToStop
#endif
        use Constants_mod, only: IK, NULL_IK
        use String_mod, only: num2str
        implicit none
        type(MaxNumDomainCheckToStop_type) :: MaxNumDomainCheckToStopObj
        MaxNumDomainCheckToStopObj%def = 10000_IK
        MaxNumDomainCheckToStopObj%null = NULL_IK
        MaxNumDomainCheckToStopObj%desc = &
        "maxNumDomainCheckToStop is an integer number beyond which the program will stop globally with a fatal error message &
        &declaring that the maximum number of proposal-out-of-domain-bounds has reached. The counter for this global-stop request &
        &is reset after a proposal point is accepted as a sample from within the domain of the objective function. &
        &The default value is " // num2str(MaxNumDomainCheckToStopObj%def) // "."
    end function constructMaxNumDomainCheckToStop

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(MaxNumDomainCheckToStopObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(MaxNumDomainCheckToStop_type), intent(in)     :: MaxNumDomainCheckToStopObj
        maxNumDomainCheckToStop = MaxNumDomainCheckToStopObj%null
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setMaxNumDomainCheckToStop(MaxNumDomainCheckToStopObj,maxNumDomainCheckToStop)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setMaxNumDomainCheckToStop
#endif
        use Constants_mod, only: IK
        implicit none
        class(MaxNumDomainCheckToStop_type), intent(inout)  :: MaxNumDomainCheckToStopObj
        integer(IK), intent(in)                             :: maxNumDomainCheckToStop
        MaxNumDomainCheckToStopObj%val = maxNumDomainCheckToStop
        if (MaxNumDomainCheckToStopObj%val==MaxNumDomainCheckToStopObj%null) then
            MaxNumDomainCheckToStopObj%val = MaxNumDomainCheckToStopObj%def
        end if
    end subroutine setMaxNumDomainCheckToStop

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine checkForSanity(MaxNumDomainCheckToStop,Err,methodName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(MaxNumDomainCheckToStop_type), intent(in) :: MaxNumDomainCheckToStop
        character(*), intent(in)                        :: methodName
        type(Err_type), intent(inout)                   :: Err
        character(*), parameter                         :: PROCEDURE_NAME = "@checkForSanity()"
        if ( MaxNumDomainCheckToStop%val<1 ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input value for variable maxNumDomainCheckToStop must be a positive integer. &
                        &If you are not sure about the appropriate value for this variable, simply drop it from the input. " // &
                        methodName // " will automatically assign an appropriate value to it.\n\n"
        end if
    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecBase_MaxNumDomainCheckToStop_mod