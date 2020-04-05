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

module SpecBase_DomainUpperLimitVec_mod

    use Constants_mod, only: RK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecBase_DomainUpperLimitVec_mod"

    real(RK), allocatable           :: domainUpperLimitVec(:) ! namelist input

    type                            :: DomainUpperLimitVec_type
        real(RK), allocatable       :: Val(:)
        real(RK)                    :: def
        real(RK)                    :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setDomainUpperLimitVec, checkForSanity, nullifyNameListVar
    end type DomainUpperLimitVec_type

    interface DomainUpperLimitVec_type
        module procedure            :: constructDomainUpperLimitVec
    end interface DomainUpperLimitVec_type

    private :: constructDomainUpperLimitVec, setDomainUpperLimitVec, checkForSanity, nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructDomainUpperLimitVec(methodName) result(DomainUpperLimitVecObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructDomainUpperLimitVec
#endif
        use Decoration_mod, only: TAB
        use Constants_mod, only: NULL_RK, POSINF_RK
        use String_mod, only: num2str
        implicit none
        type(DomainUpperLimitVec_type)  :: DomainUpperLimitVecObj
        character(*), intent(in)        :: methodName
        DomainUpperLimitVecObj%def  = POSINF_RK
        DomainUpperLimitVecObj%null = NULL_RK
        DomainUpperLimitVecObj%desc = &
        "domainUpperLimitVec represents the upper boundaries of the cubical domain of the objective function to be sampled. &
        &It is an ndim-dimensional vector of 64-bit real numbers, where ndim is the number of variables of the objective function. &
        &It is also possible to assign only select values of domainUpperLimitVec and leave the rest of the components to be assigned &
        &the default value. This is POSSIBLE ONLY when domainUpperLimitVec is defined inside the input file to "//methodName//". &
        &For example,\n\n&
        &    domainUpperLimitVec(3:5) = 100\n\n&
        &            will only set the upper limits of the third, fourth, and the fifth dimensions to 100, or,\n\n&
        &    domainUpperLimitVec(1) = 100, domainUpperLimitVec(2) = 1.e6 \n\n&
        &            will set the upper limit on the first dimension to 100, and 1.e6 on the second dimension, or,\n\n&
        &    domainUpperLimitVec = 3*2.5e100\n\n&
        &            will only set the upper limits on the first, second, and the third dimensions to 2.5*10^100, while the rest of &
                    &the upper limits for the missing dimensions will be automatically set to the default value.\n\n&
        &The default value for all elements of domainUpperLimitVec is: " // num2str(DomainUpperLimitVecObj%def) // "."
    end function constructDomainUpperLimitVec

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(DomainUpperLimitVecObj,nd)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use Constants_mod, only: IK
        implicit none
        class(DomainUpperLimitVec_type), intent(in) :: DomainUpperLimitVecObj
        integer(IK), intent(in)                     :: nd
        if (allocated(domainUpperLimitVec)) deallocate(domainUpperLimitVec)
        allocate(domainUpperLimitVec(nd))
        domainUpperLimitVec(1:nd) = DomainUpperLimitVecObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setDomainUpperLimitVec(DomainUpperLimitVecObj,domainUpperLimitVec)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setDomainUpperLimitVec
#endif
        use Constants_mod, only: IK, RK
        implicit none
        class(DomainUpperLimitVec_type), intent(inout)  :: DomainUpperLimitVecObj
        real(RK), intent(in)                            :: domainUpperLimitVec(:)
        DomainUpperLimitVecObj%Val = domainUpperLimitVec
        where ( DomainUpperLimitVecObj%Val == DomainUpperLimitVecObj%null )
            DomainUpperLimitVecObj%Val = DomainUpperLimitVecObj%def
        end where
    end subroutine setDomainUpperLimitVec

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine checkForSanity(DomainUpperLimitVecObj,Err,domainLowerLimitVec)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Err_mod, only: Err_type
        use Constants_mod, only: IK, POSINF_RK
        use String_mod, only: num2str
        implicit none
        class(DomainUpperLimitVec_type), intent(in) :: DomainUpperLimitVecObj
        real(RK), intent(in)                        :: domainLowerLimitVec(1:size(DomainUpperLimitVecObj%Val))
        type(Err_type), intent(inout)               :: Err
        character(*), parameter                     :: PROCEDURE_NAME = MODULE_NAME//"@checkForSanity()"
        integer(IK)                                 :: i
        do i = 1,size(DomainUpperLimitVecObj%Val(:))
            if ( DomainUpperLimitVecObj%Val(i) > POSINF_RK ) then
                Err%occurred = .true.
                Err%msg =   Err%msg // &
                            PROCEDURE_NAME // ": Error occurred. &
                            &The upper limit of the component " // &
                            num2str(i) // " of the variable domainUpperLimitVec (" // num2str(DomainUpperLimitVecObj%Val(i)) // ") &
                            &cannot be larger than the largest positive real number representable in the simulation (" // &
                            num2str(POSINF_RK) // ").\n\n"
            end if
            if ( DomainUpperLimitVecObj%Val(i)<=domainLowerLimitVec(i) ) then
                Err%occurred = .true.
                Err%msg =   Err%msg // &
                            PROCEDURE_NAME // ": Error occurred. &
                            &The input value for the upper limit of the component " // &
                            num2str(i) // " of the variable domainUpperLimitVec cannot be smaller than or equal to the input value &
                            &for the lower limit of the corresponding dimension as given by domainLowerLimitVec:\n" // &
                            "    domainLowerLimitVec(" // num2str(i) // ") = " // num2str(domainLowerLimitVec(i)) // "\n" // &
                            "    domainUpperLimitVec(" // num2str(i) // ") = " // num2str(DomainUpperLimitVecObj%Val(i)) // "\n\n"
            end if
        end do
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecBase_DomainUpperLimitVec_mod