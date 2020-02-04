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

module SpecBase_DomainLowerLimitVec_mod

    use Constants_mod, only: RK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecBase_DomainLowerLimitVec_mod"

    real(RK), allocatable           :: domainLowerLimitVec(:) ! namelist input

    type                            :: domainLowerLimitVec_type
        real(RK), allocatable       :: Val(:)
        real(RK)                    :: def
        real(RK)                    :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setDomainLowerLimitVec, checkForSanity, nullifyNameListVar
    end type DomainLowerLimitVec_type

    interface DomainLowerLimitVec_type
        module procedure            :: constructDomainLowerLimitVec
    end interface DomainLowerLimitVec_type

    private :: constructDomainLowerLimitVec, setDomainLowerLimitVec, checkForSanity, nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructDomainLowerLimitVec(methodName) result(DomainLowerLimitVecObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructDomainLowerLimitVec
#endif
        use Decoration_mod, only: TAB
        use Constants_mod, only: NULL_RK, NEGINF_RK
        use String_mod, only: num2str
        implicit none
        type(DomainLowerLimitVec_type)  :: DomainLowerLimitVecObj
        character(*), intent(in)        :: methodName
        DomainLowerLimitVecObj%def  = NEGINF_RK
        DomainLowerLimitVecObj%null = NULL_RK
        DomainLowerLimitVecObj%desc = &
        "domainLowerLimitVec represents the lower boundaries of the cubical domain of the objective function to be sampled. &
        &It is an ndim-dimensional vector of 64-bit real numbers, where ndim is the number of variables of the objective function. &
        &It is also possible to assign only select values of domainLowerLimitVec and leave the rest of the components to be assigned &
        &the default value. This is POSSIBLE ONLY when domainLowerLimitVec is defined inside the input file to "//methodName//". &
        &For example, having the following inside the input file, \n\n&
        &    domainLowerLimitVec(3:5) = -100\n\n&
        &            will only set the lower limits of the third, fourth, and the fifth dimensions to -100, or,\n\n&
        &    domainLowerLimitVec(1) = -100, domainLowerLimitVec(2) = -1.e6 \n\n&
        &            will set the lower limit on the first dimension to -100, and 1.e6 on the second dimension, or,\n\n&
        &    domainLowerLimitVec = 3*-2.5e100\n\n&
        &            will only set the lower limits on the first, second, and the third dimensions to -2.5*10^100, while the rest of &
                    &the lower limits for the missing dimensions will be automatically set to the default value.\n\n&
        &The default value for all elements of domainLowerLimitVec is: " // num2str(DomainLowerLimitVecObj%def) // "."
    end function constructDomainLowerLimitVec

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(DomainLowerLimitVecObj,nd)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use Constants_mod, only: IK
        implicit none
        class(DomainLowerLimitVec_type), intent(in)    :: DomainLowerLimitVecObj
        integer(IK), intent(in)                     :: nd
        if (allocated(domainLowerLimitVec)) deallocate(domainLowerLimitVec)
        allocate(domainLowerLimitVec(nd))
        domainLowerLimitVec(1:nd) = DomainLowerLimitVecObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setDomainLowerLimitVec(DomainLowerLimitVecObj,domainLowerLimitVec)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setDomainLowerLimitVec
#endif
        use Constants_mod, only: IK, RK
        implicit none
        class(DomainLowerLimitVec_type), intent(inout)  :: DomainLowerLimitVecObj
        real(RK), intent(in)                            :: domainLowerLimitVec(:)
        DomainLowerLimitVecObj%Val = domainLowerLimitVec
        where ( DomainLowerLimitVecObj%Val == DomainLowerLimitVecObj%null )
            DomainLowerLimitVecObj%Val = DomainLowerLimitVecObj%def
        end where
    end subroutine setDomainLowerLimitVec

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine checkForSanity(DomainLowerLimitVecObj,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Err_mod, only: Err_type
        use Constants_mod, only: IK, NEGINF_RK
        use String_mod, only: num2str
        implicit none
        class(DomainLowerLimitVec_type), intent(in)     :: DomainLowerLimitVecObj
        type(Err_type), intent(inout)                   :: Err
        character(*), parameter                         :: PROCEDURE_NAME = MODULE_NAME//"@checkForSanity()"
        integer(IK)                                     :: i
        do i = 1,size(DomainLowerLimitVecObj%Val(:))
            if ( DomainLowerLimitVecObj%Val(i) < NEGINF_RK ) then
                Err%occurred = .true.
                Err%msg =   Err%msg // &
                            PROCEDURE_NAME // ": Error occurred. The component " // num2str(i) // &
                            " of the variable domainLowerLimitVec (" // num2str(DomainLowerLimitVecObj%Val(i)) // ") &
                            &cannot be smaller than the smallest positive real number representable in the simulation (" // &
                            num2str(NEGINF_RK) // ").\n\n"
            end if
        end do
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecBase_DomainLowerLimitVec_mod