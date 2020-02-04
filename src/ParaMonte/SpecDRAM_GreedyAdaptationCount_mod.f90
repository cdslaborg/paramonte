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

module SpecDRAM_GreedyAdaptationCount_mod

    use Constants_mod, only: IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecDRAM_GreedyAdaptationCount_mod"

    integer(IK)                     :: greedyAdaptationCount ! namelist input

    type                            :: GreedyAdaptationCount_type
        integer(IK)                 :: val
        integer(IK)                 :: def
        integer(IK)                 :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setGreedyAdaptationCount, checkForSanity, nullifyNameListVar
    end type GreedyAdaptationCount_type

    interface GreedyAdaptationCount_type
        module procedure            :: constructGreedyAdaptationCount
    end interface GreedyAdaptationCount_type

    private :: constructGreedyAdaptationCount, setGreedyAdaptationCount, checkForSanity, nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructGreedyAdaptationCount(methodName) result(GreedyAdaptationCountObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructGreedyAdaptationCount
#endif
        use Constants_mod, only: IK, NULL_IK
        use String_mod, only: num2str
        implicit none
        character(*), intent(in)            :: methodName
        type(GreedyAdaptationCount_type)    :: GreedyAdaptationCountObj
        GreedyAdaptationCountObj%def    = 0
        GreedyAdaptationCountObj%null   = NULL_IK
        GreedyAdaptationCountObj%desc   = &
        "If greedyAdaptationCount is set to a positive integer then the first greedyAdaptationCount number of &
        &the adaptive updates of the sampler will be made using only the 'unique' accepted points in the MCMC chain. &
        &This is useful for example, the function to be sampled by " // methodName // " is high dimensional, &
        &in which case, the adaptive updates to " // methodName // "'s sampler distribution will less likely lead to &
        &numerical instabilities, for example, a singular covariance matrix for the multivariate proposal sampler. &
        &The variable greedyAdaptationCount must be a non-negative integer, and not larger than the value of adaptiveUpdateCount. &
        &If it is larger, it will be automatically set to adaptiveUpdateCount for the simulation. &
        &The default value is " // num2str(GreedyAdaptationCountObj%def) // "."
    end function constructGreedyAdaptationCount

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(GreedyAdaptationCountObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(GreedyAdaptationCount_type), intent(in)  :: GreedyAdaptationCountObj
        greedyAdaptationCount = GreedyAdaptationCountObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setGreedyAdaptationCount(GreedyAdaptationCountObj,greedyAdaptationCount)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setGreedyAdaptationCount
#endif
        use Constants_mod, only: IK
        implicit none
        class(GreedyAdaptationCount_type), intent(inout)  :: GreedyAdaptationCountObj
        integer(IK), intent(in)                         :: greedyAdaptationCount
        GreedyAdaptationCountObj%val = greedyAdaptationCount
        if ( GreedyAdaptationCountObj%val==GreedyAdaptationCountObj%null ) then
            GreedyAdaptationCountObj%val = GreedyAdaptationCountObj%def
        end if
    end subroutine setGreedyAdaptationCount

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine checkForSanity(GreedyAdaptationCountObj,Err,methodName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: IK
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(GreedyAdaptationCount_type), intent(in)   :: GreedyAdaptationCountObj
        character(*), intent(in)            :: methodName
        type(Err_type), intent(inout)       :: Err
        character(*), parameter             :: PROCEDURE_NAME = "@checkForSanity()"
        if ( GreedyAdaptationCountObj%val<0_IK) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input requested value for greedyAdaptationCount (" // num2str(GreedyAdaptationCountObj%val) // ") &
                        &can not be negative. If you are not sure of the appropriate value for greedyAdaptationCount, drop it &
                        &from the input list. " // methodName // " will automatically assign an appropriate value to it.\n\n"
        end if
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecDRAM_GreedyAdaptationCount_mod