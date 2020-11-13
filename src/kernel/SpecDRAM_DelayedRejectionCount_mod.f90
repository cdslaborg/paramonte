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

module SpecDRAM_DelayedRejectionCount_mod

    use Constants_mod, only: IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecDRAM_DelayedRejectionCount_mod"
    integer(IK) , parameter         :: MAX_DELAYED_REJECTION_COUNT = 1000_IK
    integer(IK) , parameter         :: MIN_DELAYED_REJECTION_COUNT = 0_IK

    integer(IK)                     :: delayedRejectionCount ! namelist input

    type                            :: DelayedRejectionCount_type
        integer(IK)                 :: val
        integer(IK)                 :: def
        integer(IK)                 :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setDelayedRejectionCount, checkForSanity, nullifyNameListVar
    end type DelayedRejectionCount_type

    interface DelayedRejectionCount_type
        module procedure            :: constructDelayedRejectionCount
    end interface DelayedRejectionCount_type

    private :: constructDelayedRejectionCount, setDelayedRejectionCount, checkForSanity, nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructDelayedRejectionCount(methodName) result(DelayedRejectionCountObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructDelayedRejectionCount
#endif
        use Constants_mod, only: IK, NULL_IK
        use String_mod, only: num2str
        use Decoration_mod, only: TAB
        implicit none
        character(*), intent(in)            :: methodName
        type(DelayedRejectionCount_type)    :: DelayedRejectionCountObj
        DelayedRejectionCountObj%def    = 0_IK
        DelayedRejectionCountObj%null   = NULL_IK
        DelayedRejectionCountObj%desc   = &
        num2str(MIN_DELAYED_REJECTION_COUNT) // " <= delayedRejectionCount <= " // num2str(MAX_DELAYED_REJECTION_COUNT) // &
        " is an integer that represents the total number of stages for which rejections of new proposals &
        &will be tolerated by "//methodName//" before going back to the previously accepted point (state). &
        &Possible values are:\n\n&
        &    delayedRejectionCount = 0\n\n&
        &            indicating no deployment of the delayed rejection algorithm.\n\n&
        &    delayedRejectionCount > 0\n\n&
        &            which implies a maximum delayedRejectionCount number of rejections will be tolerated.\n\n&
        &For example, delayedRejectionCount = 1, means that at any point during the sampling, if a proposal is rejected, "&
        //methodName//" will not go back to the last sampled state. Instead, it will continue to propose a new from the current &
        &rejected state. If the new state is again rejected based on the rules of "//methodName//", then the algorithm will not &
        &tolerate further rejections, because the maximum number of rejections to be tolerated has been set by the user to be &
        &delayedRejectionCount = 1. The algorithm then goes back to the original last-accepted state and will begin proposing &
        &new states from that location. &
        &The default value is delayedRejectionCount = "// num2str(DelayedRejectionCountObj%def) // "."
    end function constructDelayedRejectionCount

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(DelayedRejectionCountObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(DelayedRejectionCount_type), intent(in)  :: DelayedRejectionCountObj
        delayedRejectionCount = DelayedRejectionCountObj%null
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setDelayedRejectionCount(DelayedRejectionCountObj,delayedRejectionCount)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setDelayedRejectionCount
#endif
        use Constants_mod, only: IK
        implicit none
        class(DelayedRejectionCount_type), intent(inout)    :: DelayedRejectionCountObj
        integer(IK), intent(in)                             :: delayedRejectionCount
        DelayedRejectionCountObj%val = delayedRejectionCount
        if (DelayedRejectionCountObj%val==DelayedRejectionCountObj%null) DelayedRejectionCountObj%val = DelayedRejectionCountObj%def
    end subroutine setDelayedRejectionCount

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine checkForSanity(DelayedRejectionCountObj,Err,methodName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: IK
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(DelayedRejectionCount_type), intent(in)   :: DelayedRejectionCountObj
        character(*), intent(in)            :: methodName
        type(Err_type), intent(inout)       :: Err
        character(*), parameter             :: PROCEDURE_NAME = "@checkForSanity()"
        if ( DelayedRejectionCountObj%val<MIN_DELAYED_REJECTION_COUNT) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input requested value for delayedRejectionCount (" // num2str(DelayedRejectionCountObj%val) // ") &
                        &can not be negative. If you are not sure of the appropriate value for delayedRejectionCount, drop it &
                        &from the input list. " // methodName // " will automatically assign an appropriate value to it.\n\n"
        elseif ( DelayedRejectionCountObj%val>MAX_DELAYED_REJECTION_COUNT) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input requested value for delayedRejectionCount (" // num2str(DelayedRejectionCountObj%val) // ") &
                        &can not be > " // num2str(MAX_DELAYED_REJECTION_COUNT) // &
                        ". If you are not sure of the appropriate value for delayedRejectionCount, drop it &
                        &from the input list. " // methodName // " will automatically assign an appropriate value to it.\n\n"
        end if
    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecDRAM_DelayedRejectionCount_mod