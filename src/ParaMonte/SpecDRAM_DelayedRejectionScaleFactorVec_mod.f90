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

module SpecDRAM_DelayedRejectionScaleFactorVec_mod

    use Constants_mod, only: RK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecDRAM_DelayedRejectionScaleFactorVec_mod"

    real(RK), allocatable           :: delayedRejectionScaleFactorVec(:) ! namelist input

    type                            :: DelayedRejectionScaleFactorVec_type
        real(RK), allocatable       :: Val(:)
        real(RK)                    :: def
        real(RK)                    :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setDelayedRejectionScaleFactorVec, checkForSanity, nullifyNameListVar
    end type DelayedRejectionScaleFactorVec_type

    interface DelayedRejectionScaleFactorVec_type
        module procedure            :: constructDelayedRejectionScaleFactorVec
    end interface DelayedRejectionScaleFactorVec_type

    private :: constructDelayedRejectionScaleFactorVec, setDelayedRejectionScaleFactorVec, checkForSanity, nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructDelayedRejectionScaleFactorVec(nd,methodName) result(DelayedRejectionScaleFactorVecObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructDelayedRejectionScaleFactorVec
#endif
        use Decoration_mod, only: TAB
        use Constants_mod, only: IK, RK, NULL_RK
        use String_mod, only: num2str
        implicit none
        integer(IK), intent(in)                     :: nd
        character(*), intent(in)                    :: methodName
        type(DelayedRejectionScaleFactorVec_type)   :: DelayedRejectionScaleFactorVecObj
        DelayedRejectionScaleFactorVecObj%def  = 0.5_RK**(1._RK/real(nd,kind=RK))  ! This gives a half volume to the covariance
        DelayedRejectionScaleFactorVecObj%null = NULL_RK
        DelayedRejectionScaleFactorVecObj%desc = &
        "delayedRejectionScaleFactorVec is a real-valued positive vector of length (1:delayedRejectionCount) by which &
        &the covariance matrix of the proposal distribution of "//methodName//" sampler is scaled when the Delayed Rejection (DR) &
        &scheme is activated (by setting delayedRejectionCount>0). At each ith stage of the DR process, &
        &the proposal distribution from the last stage is scaled by the factor delayedRejectionScaleFactorVec(i). &
        &Missing elements of the delayedRejectionScaleFactorVec in the input to "//methodName//" will be set to the default value. &
        &The default value at all stages is 0.5^(1/ndim) = "//num2str(DelayedRejectionScaleFactorVecObj%def)//", which reduces the &
        &volume of the covariance matrix of the proposal from the last DR stage by one half. The variable ndim represents the &
        &number of dimensions of the Domain of the objective function."
    end function constructDelayedRejectionScaleFactorVec

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(DelayedRejectionScaleFactorVecObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use SpecDRAM_DelayedRejectionCount_mod, only: MAX_DELAYED_REJECTION_COUNT
        use Constants_mod, only: IK
        implicit none
        class(DelayedRejectionScaleFactorVec_type), intent(in) :: DelayedRejectionScaleFactorVecObj
        if (allocated(delayedRejectionScaleFactorVec)) deallocate(delayedRejectionScaleFactorVec)
        allocate(delayedRejectionScaleFactorVec(MAX_DELAYED_REJECTION_COUNT))
        delayedRejectionScaleFactorVec = DelayedRejectionScaleFactorVecObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setDelayedRejectionScaleFactorVec(DelayedRejectionScaleFactorVecObj,delayedRejectionScaleFactorVec,delayedRejectionCount)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setDelayedRejectionScaleFactorVec
#endif
        use Constants_mod, only: IK, RK
        implicit none
        class(DelayedRejectionScaleFactorVec_type), intent(inout)   :: DelayedRejectionScaleFactorVecObj
        real(RK)    , intent(in)                                    :: delayedRejectionScaleFactorVec(:)
        integer(IK) , intent(in)                                    :: delayedRejectionCount
        integer(IK)                                                 :: i
        DelayedRejectionScaleFactorVecObj%Val = pack( delayedRejectionScaleFactorVec &
                                                    , mask = delayedRejectionScaleFactorVec /= DelayedRejectionScaleFactorVecObj%null )
        if ( size(DelayedRejectionScaleFactorVecObj%Val) == 0 .and. delayedRejectionCount > 0 ) then
            deallocate(DelayedRejectionScaleFactorVecObj%Val)
            allocate(DelayedRejectionScaleFactorVecObj%Val(delayedRejectionCount))
            do i = 1, delayedRejectionCount
                DelayedRejectionScaleFactorVecObj%Val(i) = DelayedRejectionScaleFactorVecObj%def
            end do
        end if
    end subroutine setDelayedRejectionScaleFactorVec

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine checkForSanity(DelayedRejectionScaleFactorVecObj,Err,methodName,delayedRejectionCount)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: IK, RK
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(DelayedRejectionScaleFactorVec_type), intent(in)  :: DelayedRejectionScaleFactorVecObj
        type(Err_type), intent(inout)                           :: Err
        character(*), intent(in)                                :: methodName
        integer(IK), intent(in)                                 :: delayedRejectionCount
        character(*), parameter                                 :: PROCEDURE_NAME = "@checkForSanity()"
        integer(IK)                                             :: i
        if ( size(DelayedRejectionScaleFactorVecObj%Val)/=delayedRejectionCount ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME//PROCEDURE_NAME//": Error occurred. The length of the vector delayedRejectionScaleFactorVec (" // &
                        num2str(size(DelayedRejectionScaleFactorVecObj%Val)) // ") is not equal to delayedRejectionCount = " // &
                        num2str(delayedRejectionCount) // ". If you are not sure how to set the values of &
                        &delayedRejectionScaleFactorVec, drop it from the input. " &
                        // methodName // " will automatically set the appropriate value for delayedRejectionScaleFactorVec.\n\n"
        end if
        do i = 1,size(DelayedRejectionScaleFactorVecObj%Val)
            if ( DelayedRejectionScaleFactorVecObj%Val(i)<=0._RK ) then
                Err%occurred = .true.
                Err%msg =   Err%msg // &
                            MODULE_NAME // PROCEDURE_NAME // ": Error occurred. The input value for the element " // &
                            num2str(i) // " of the variable delayedRejectionScaleFactorVec cannot be smaller than or equal to 0.\n\n"
            end if
        end do
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecDRAM_DelayedRejectionScaleFactorVec_mod