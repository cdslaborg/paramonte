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
!!!!       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!> \brief
!> This module contains the classes and procedures for setting up the `delayedRejectionScaleFactorVec` attribute of samplers of class [ParaDRAM_type](@ref paradram_mod::paradram_type).
!> For more information, see the description of this attribute in the body of the module.
!> \author Amir Shahmoradi

module SpecDRAM_DelayedRejectionScaleFactorVec_mod

    use Constants_mod, only: RK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecDRAM_DelayedRejectionScaleFactorVec_mod"

    real(RK), allocatable           :: DelayedRejectionScaleFactorVec(:) ! namelist input

    type                            :: DelayedRejectionScaleFactorVec_type
        real(RK), allocatable       :: Val(:)
        real(RK)                    :: def
        real(RK)                    :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set, checkForSanity, nullifyNameListVar
    end type DelayedRejectionScaleFactorVec_type

    interface DelayedRejectionScaleFactorVec_type
        module procedure            :: construct
    end interface DelayedRejectionScaleFactorVec_type

    private :: construct, set, checkForSanity, nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function construct(nd,methodName) result(self)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: construct
#endif
        use Decoration_mod, only: TAB
        use Constants_mod, only: IK, RK, NULL_RK
        implicit none
        integer(IK), intent(in)                     :: nd
        character(*), intent(in)                    :: methodName
        type(DelayedRejectionScaleFactorVec_type)   :: self
        self%def  = 0.5_RK**(1._RK/real(nd,kind=RK))  ! This gives a half volume to the covariance
        self%null = NULL_RK
        self%desc = &
        "delayedRejectionScaleFactorVec is a real-valued positive vector of length (1:delayedRejectionCount) by which &
        &the covariance matrix of the proposal distribution of "//methodName//" sampler is scaled when the Delayed Rejection (DR) &
        &scheme is activated (by setting delayedRejectionCount>0). At each ith stage of the DR process, &
        &the proposal distribution from the last stage is scaled by the factor delayedRejectionScaleFactorVec(i). &
        &Missing elements of the delayedRejectionScaleFactorVec in the input to "//methodName//" will be set to the default value. &
        &The default value at all stages is 0.5^(1/ndim) where ndim is the number of dimensions of the domain of the objective &
        &function. This default value effectively reduces the volume of the covariance matrix of the proposal distribution by &
        &half compared to the last DR stage."
    end function construct

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(self)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use SpecDRAM_DelayedRejectionCount_mod, only: MAX_DELAYED_REJECTION_COUNT
        use Constants_mod, only: IK
        implicit none
        class(DelayedRejectionScaleFactorVec_type), intent(in) :: self
        if (allocated(DelayedRejectionScaleFactorVec)) deallocate(DelayedRejectionScaleFactorVec)
        allocate(DelayedRejectionScaleFactorVec(MAX_DELAYED_REJECTION_COUNT), source = self%null)
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure subroutine set(self,delayedRejectionCount, DelayedRejectionScaleFactorVec)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: set
#endif
        use Constants_mod, only: IK, RK
        implicit none
        class(DelayedRejectionScaleFactorVec_type), intent(inout)   :: self
        integer(IK) , intent(in)                                    :: delayedRejectionCount
        real(RK)    , intent(in), optional                          :: DelayedRejectionScaleFactorVec(:) ! This must remain colon
        if (present(DelayedRejectionScaleFactorVec)) then
            self%Val = pack( DelayedRejectionScaleFactorVec, mask = DelayedRejectionScaleFactorVec /= self%null )
            if (size(self%Val) == 0) then
                deallocate(self%Val)
                allocate(self%Val(delayedRejectionCount), source = self%def)
            end if
        else
            if (allocated(self%Val)) then
                if (size(self%Val) == 0) then
                    deallocate(self%Val)
                    allocate(self%Val(delayedRejectionCount), source = self%def)
                end if
            end if
        end if
    end subroutine set

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine checkForSanity(self,Err,methodName,delayedRejectionCount)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: IK, RK
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(DelayedRejectionScaleFactorVec_type), intent(in)  :: self
        type(Err_type), intent(inout)                           :: Err
        character(*), intent(in)                                :: methodName
        integer(IK), intent(in)                                 :: delayedRejectionCount
        character(*), parameter                                 :: PROCEDURE_NAME = "@checkForSanity()"
        integer(IK)                                             :: i
        if ( size(self%Val)/=delayedRejectionCount ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME//PROCEDURE_NAME//": Error occurred. The length of the vector delayedRejectionScaleFactorVec (" // &
                        num2str(size(self%Val)) // ") is not equal to delayedRejectionCount = " // &
                        num2str(delayedRejectionCount) // ". If you are not sure how to set the values of &
                        &delayedRejectionScaleFactorVec, drop it from the input. " &
                        // methodName // " will automatically set the appropriate value for delayedRejectionScaleFactorVec.\n\n"
        end if
        do i = 1,size(self%Val)
            if ( self%Val(i)<=0._RK ) then
                Err%occurred = .true.
                Err%msg =   Err%msg // &
                            MODULE_NAME // PROCEDURE_NAME // ": Error occurred. The input value for the element " // &
                            num2str(i) // " of the variable delayedRejectionScaleFactorVec cannot be smaller than or equal to 0.\n\n"
            end if
        end do
        deallocate(DelayedRejectionScaleFactorVec)
    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecDRAM_DelayedRejectionScaleFactorVec_mod ! LCOV_EXCL_LINE