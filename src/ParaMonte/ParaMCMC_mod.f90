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

module ParaMCMC_mod

    use Constants_mod, only: RK, IK
    use SpecMCMC_mod, only: SpecMCMC_type
    use ParaMonte_mod, only: ParaMonte_type, ParaMonteLogFuncMode_type, ParaMonteStatistics_type, Moment_type
    use ParaMCMCRefinedChain_mod, only: RefinedChain_type
    implicit none

    character(*), parameter                     :: MODULE_NAME = "@ParaMCMC_mod"

    type                                        :: ParaMCMC_Chain_type
        integer(IK)                             :: compact, verbose
    end type ParaMCMC_Chain_type

    type, extends(ParaMonteLogFuncMode_type)    :: ParaMCMC_LogFuncMode_type
        type(ParaMCMC_Chain_type)               :: Loc
    end type ParaMCMC_LogFuncMode_type

    type, extends(ParaMonteStatistics_type)     :: ParaMCMC_Statistics_type
        type(Moment_type)                       :: Chain
        type(ParaMCMC_Chain_type)               :: BurninLoc
        type(ParaMCMC_LogFuncMode_type)         :: LogFuncMode
    end type ParaMCMC_Statistics_type

    type, extends(ParaMonte_type)               :: ParaMCMC_type
       type(RefinedChain_type)                  :: RefinedChain
       type(SpecMCMC_type)                      :: SpecMCMC
    contains
        procedure, pass                         :: setupParaMCMC
    end type ParaMCMC_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setupParaMCMC(PMPM)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setupParaMCMC
#endif
        use SpecMCMC_mod, only: SpecMCMC_type
        implicit none
        class(ParaMCMC_type), intent(inout)    :: PMPM
        character(*), parameter :: PROCEDURE_NAME = "@setupParaMCMC()"
        PMPM%SpecMCMC = SpecMCMC_type(nd=PMPM%nd%val,methodName=PMPM%name)
    end subroutine setupParaMCMC

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module ParaMCMC_mod