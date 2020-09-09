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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        "targetAcceptanceRate sets an optimal target for the ratio of the number of accepted objective function calls to the &
        &total number of function calls by the " // methodName // " sampler. It is a real-valued array of length 2, whose elements &
        &determine the upper and lower bounds of the desired acceptance rate. When the acceptance rate of the sampler is outside the &
        &specified limits, the sampler's settings will be automatically adjusted to bring the overall acceptance rate to within the &
        &specified limits by the input variable targetAcceptanceRate. When assigned from within a dynamic-language programming &
        &environment, such as MATLAB or Python, or from within an input file, targetAcceptanceRate can also be a single real number &
        &between 0 and 1. In such case, the " // methodName // " sampler will constantly attempt (with no guarantee of success) &
        &to bring the average acceptance ratio of the sampler as close to the user-provided target ratio as possible. The success &
        &of " // methodName // " in keeping the average acceptance ratio close to the requested target value depends heavily on:\n\n&
        &    1) the value of adaptiveUpdatePeriod; the larger, the easier.\n&
        &    2) the value of adaptiveUpdateCount; the larger, the easier.\n\n&
        &Note that the acceptance ratio adjustments will only occur every adaptiveUpdatePeriod sampling steps for a total &
        &number of adaptiveUpdateCount. &
        &There is no default value for targetAcceptanceRate, as the acceptance ratio is not directly adjusted during sampling."
    end function constructTargetAcceptanceRate

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(TargetAcceptanceRateObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(TargetAcceptanceRate_type), intent(in)  :: TargetAcceptanceRateObj
        TargetAcceptanceRate = TargetAcceptanceRateObj%null
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecBase_TargetAcceptanceRate_mod