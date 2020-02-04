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

module SpecMCMC_RandomStartPointDomainLowerLimitVec_mod

    use Constants_mod, only: RK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecMCMC_RandomStartPointDomainLowerLimitVec_mod"

    real(RK), allocatable           :: RandomStartPointDomainLowerLimitVec(:) ! namelist input

    type                            :: RandomStartPointDomainLowerLimitVec_type
        real(RK), allocatable       :: Val(:)
        real(RK)                    :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setRandomStartPointDomainLowerLimitVec, checkForSanity, nullifyNameListVar
    end type RandomStartPointDomainLowerLimitVec_type

    interface RandomStartPointDomainLowerLimitVec_type
        module procedure            :: constructRandomStartPointDomainLowerLimitVec
    end interface RandomStartPointDomainLowerLimitVec_type

    private :: constructRandomStartPointDomainLowerLimitVec, setRandomStartPointDomainLowerLimitVec, checkForSanity, nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructRandomStartPointDomainLowerLimitVec(methodName) result(RandomStartPointDomainLowerLimitVecObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructRandomStartPointDomainLowerLimitVec
#endif
        use Constants_mod, only: NULL_RK
        use String_mod, only: num2str
        implicit none
        character(*), intent(in)                        :: methodName
        type(RandomStartPointDomainLowerLimitVec_type)  :: RandomStartPointDomainLowerLimitVecObj
        RandomStartPointDomainLowerLimitVecObj%null = NULL_RK
        RandomStartPointDomainLowerLimitVecObj%desc = &
        "RandomStartPointDomainLowerLimitVec represents the lower boundaries of the cubical domain from which the starting point(s) of &
        &the MCMC chain(s) will be initialized randomly (only if requested via the input variable randomStartPointRequested. &
        &This happens only when some or all of the elements of the input variable StartPoint are missing. &
        &In such cases, every missing value of input StartPoint will be set to the center point between RandomStartPointDomainLowerLimitVec &
        &and RandomStartPointDomainUpperLimit in the corresponding dimension. &
        &If RandomStartPointRequested=TRUE (or True, true, t, all case-insensitive), then the missing &
        &elements of StartPoint will be initialized to values drawn randomly from within the corresponding ranges specified by &
        &the input variable RandomStartPointDomainLowerLimitVec. &
        &As an input variable, RandomStartPointDomainLowerLimitVec is an ndim-dimensional vector of 64-bit real numbers, &
        &where ndim is the number of variables of the objective function. It is also possible to assign only select values of &
        &RandomStartPointDomainLowerLimitVec and leave the rest of the components to be assigned the default value. &
        &This is POSSIBLE ONLY when RandomStartPointDomainLowerLimitVec is defined inside the input file to "//methodName//". &
        &For example, having the following inside the input file, \n\n&
        &    RandomStartPointDomainLowerLimitVec(3:5) = -100\n\n&
        &            will only set the lower limits of the third, fourth, and the fifth dimensions to -100, or,\n\n&
        &    RandomStartPointDomainLowerLimitVec(1) = -100, RandomStartPointDomainLowerLimitVec(2) = -1.e6 \n\n&
        &            will set the lower limit on the first dimension to -100, and 1.e6 on the second dimension, or,\n\n&
        &    RandomStartPointDomainLowerLimitVec = 3*-2.5e100\n\n&
        &            will only set the lower limits on the first, second, and the third dimensions to -2.5*10^100, while the rest of &
                    &the lower limits for the missing dimensions will be automatically set to the default value.\n\n&
        &The default values for all elements of RandomStartPointDomainLowerLimitVec are taken from the corresponding values in the input &
        &variable domainLowerLimitVec."
    end function constructRandomStartPointDomainLowerLimitVec

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(RandomStartPointDomainLowerLimitVecObj,nd)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use Constants_mod, only: IK
        implicit none
        class(RandomStartPointDomainLowerLimitVec_type), intent(in)  :: RandomStartPointDomainLowerLimitVecObj
        integer(IK), intent(in)                             :: nd
        if (allocated(randomStartPointDomainLowerLimitVec)) deallocate(randomStartPointDomainLowerLimitVec)
        allocate(randomStartPointDomainLowerLimitVec(nd))
        randomStartPointDomainLowerLimitVec = RandomStartPointDomainLowerLimitVecObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setRandomStartPointDomainLowerLimitVec(RandomStartPointDomainLowerLimitVecObj,randomStartPointDomainLowerLimitVec,domainLowerLimitVec)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setRandomStartPointDomainLowerLimitVec
#endif
        use Constants_mod, only: IK, RK
        implicit none
        class(RandomStartPointDomainLowerLimitVec_type), intent(inout)  :: RandomStartPointDomainLowerLimitVecObj
        real(RK), intent(in)                                            :: randomStartPointDomainLowerLimitVec(:)
        real(RK), intent(in)                                            :: domainLowerLimitVec(:)
        RandomStartPointDomainLowerLimitVecObj%Val = randomStartPointDomainLowerLimitVec
        where (RandomStartPointDomainLowerLimitVecObj%Val==RandomStartPointDomainLowerLimitVecObj%null)
            RandomStartPointDomainLowerLimitVecObj%Val = domainLowerLimitVec
        end where
    end subroutine setRandomStartPointDomainLowerLimitVec

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine checkForSanity(RandomStartPointDomainLowerLimitVecObj,Err,methodName,domainLowerLimitVec)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: RK
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(RandomStartPointDomainLowerLimitVec_type), intent(in) :: RandomStartPointDomainLowerLimitVecObj
        real(RK), intent(in)                                        :: domainLowerLimitVec(:)
        character(*), intent(in)                                    :: methodName
        type(Err_type), intent(inout)                               :: Err
        character(*), parameter                                     :: PROCEDURE_NAME = "@checkForSanity()"
        integer                                                     :: i
        do i = 1,size(RandomStartPointDomainLowerLimitVecObj%Val(:))
            if ( RandomStartPointDomainLowerLimitVecObj%Val(i)<domainLowerLimitVec(i) ) then
                Err%occurred = .true.
                Err%msg =   Err%msg // &
                            MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                            &The component " // num2str(i) // " of the variable RandomStartPointDomainLowerLimitVec (" // &
                            num2str(RandomStartPointDomainLowerLimitVecObj%Val(i)) // &
                            ") cannot be smaller than the corresponding component of the variable &
                            &domainLowerLimitVec (" // num2str(domainLowerLimitVec(i)) // "). If you don't know &
                            &an appropriate value to set for RandomStartPointDomainLowerLimitVec, drop it from the input list. " // &
                            methodName // " will automatically assign an appropriate value to it.\n\n"
            end if
        end do
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecMCMC_RandomStartPointDomainLowerLimitVec_mod