!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!   ParaMonte: plain powerful parallel Monte Carlo library.
!
!   Copyright (C) 2012-present, The Computational Data Science Lab
!
!   This file is part of the ParaMonte library.
!
!   ParaMonte is free software: you can redistribute it and/or modify it
!   under the terms of the GNU Lesser General Public License as published
!   by the Free Software Foundation, version 3 of the License.
!
!   ParaMonte is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the ParaMonte library. If not, see,
!
!       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
!
!   ACKNOWLEDGMENT
!
!   As per the ParaMonte library license agreement terms,
!   if you use any parts of this library for any purposes,
!   we ask you to acknowledge the use of the ParaMonte library
!   in your work (education/research/industry/development/...)
!   by citing the ParaMonte library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module SpecDRAM_AdaptiveUpdateCount_mod

    use Constants_mod, only: IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecDRAM_AdaptiveUpdateCount_mod"

    integer(IK)                     :: adaptiveUpdateCount ! namelist input

    type                            :: AdaptiveUpdateCount_type
        integer(IK)                 :: val
        integer(IK)                 :: def
        integer(IK)                 :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setAdaptiveUpdateCount, checkForSanity, nullifyNameListVar
    end type AdaptiveUpdateCount_type

    interface AdaptiveUpdateCount_type
        module procedure            :: constructAdaptiveUpdateCount
    end interface AdaptiveUpdateCount_type

    private :: constructAdaptiveUpdateCount, setAdaptiveUpdateCount, checkForSanity, nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   !function constructAdaptiveUpdateCount(methodName,chainSizeDef,adaptiveUpdatePeriodDef) result(AdaptiveUpdateCountObj)
    function constructAdaptiveUpdateCount(methodName) result(AdaptiveUpdateCountObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructAdaptiveUpdateCount
#endif
        use Constants_mod, only: IK, NULL_IK, POSINF_IK
        use String_mod, only: num2str
        implicit none
        character(*), intent(in)        :: methodName
        type(AdaptiveUpdateCount_type)  :: AdaptiveUpdateCountObj
       !integer(IK), intent(in)         :: chainSizeDef, adaptiveUpdatePeriodDef
       !AdaptiveUpdateCountObj%def  = chainSizeDef / ( 2 * adaptiveUpdatePeriodDef )  
        AdaptiveUpdateCountObj%def  = POSINF_IK
        AdaptiveUpdateCountObj%null = NULL_IK
        AdaptiveUpdateCountObj%desc = &
        "adaptiveUpdateCount represents the total number of adaptive updates that will be made &
        &to the parameters of the proposal distribution, to increase the efficiency of the sampler &
        &thus increasing the sampling efficiency of " // methodName // ". &
        &Every adaptiveUpdatePeriod number of calls to the objective function, the parameters of the proposal distribution &
        &will be updated until either the total number of adaptive updates reaches the value of adaptiveUpdateCount. &
        &This variable must be a non-negative integer. As a rule of thumb, it may be appropriate to &
        &set the input variable chainSize > 2 * adaptiveUpdatePeriod * adaptiveUpdateCount, &
        &to ensure ergodicity and stationarity of the MCMC sampler. &
        &If adaptiveUpdateCount=0, then the proposal distribution parameters will be fixed to &
        &the initial input values throughout the entire MCMC sampling. &
        &The default value is " // num2str(AdaptiveUpdateCountObj%def) // "."
    end function constructAdaptiveUpdateCount

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(AdaptiveUpdateCountObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(AdaptiveUpdateCount_type), intent(in)  :: AdaptiveUpdateCountObj
        adaptiveUpdateCount = AdaptiveUpdateCountObj%null
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setAdaptiveUpdateCount(AdaptiveUpdateCountObj,adaptiveUpdateCount)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setAdaptiveUpdateCount
#endif
        use Constants_mod, only: IK
        implicit none
        class(AdaptiveUpdateCount_type), intent(inout)  :: AdaptiveUpdateCountObj
        integer(IK), intent(in)                         :: adaptiveUpdateCount
        AdaptiveUpdateCountObj%val = adaptiveUpdateCount
        if ( AdaptiveUpdateCountObj%val==AdaptiveUpdateCountObj%null ) AdaptiveUpdateCountObj%val = AdaptiveUpdateCountObj%def
    end subroutine setAdaptiveUpdateCount

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine checkForSanity(AdaptiveUpdateCountObj,Err,methodName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Constants_mod, only: IK
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(AdaptiveUpdateCount_type), intent(in)   :: AdaptiveUpdateCountObj
        character(*), intent(in)            :: methodName
        type(Err_type), intent(inout)       :: Err
        character(*), parameter             :: PROCEDURE_NAME = "@checkForSanity()"
        if ( AdaptiveUpdateCountObj%val<0) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input requested value for adaptiveUpdateCount (" // num2str(AdaptiveUpdateCountObj%val) // ") &
                        &can not be negative. If you are not sure of the appropriate value for adaptiveUpdateCount, drop it &
                        &from the input list. " // methodName // " will automatically assign an appropriate value to it.\n\n"
        end if
    end subroutine checkForSanity

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecDRAM_AdaptiveUpdateCount_mod