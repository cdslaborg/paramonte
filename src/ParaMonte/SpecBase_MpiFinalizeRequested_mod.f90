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

module SpecBase_MpiFinalizeRequested_mod

    implicit none

    ! namelist input
    logical                         :: mpiFinalizeRequested

    type                            :: MpiFinalizeRequested_type
        logical                     :: val
        logical                     :: def
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setMpiFinalizeRequested, nullifyNameListVar
    end type MpiFinalizeRequested_type

    interface MpiFinalizeRequested_type
        module procedure            :: constructMpiFinalizeRequested
    end interface MpiFinalizeRequested_type

    private :: constructMpiFinalizeRequested, setMpiFinalizeRequested

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructMpiFinalizeRequested(methodName) result(MpiFinalizeRequestedObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructMpiFinalizeRequested
#endif
        use String_mod, only: log2str
        implicit none
        character(*), intent(in)        :: methodName
        type(MpiFinalizeRequested_type) :: MpiFinalizeRequestedObj
        MpiFinalizeRequestedObj%def = .true.
        MpiFinalizeRequestedObj%desc = &
        "In parallel " // methodName // " simulations via MPI communication libraries, &
        &if mpiFinalizeRequested = true (or T, both case-insensitive), then a call will be made to the MPI_Finalize() routine &
        &from inside " // methodName // " at the end of the simulation to finalize the MPI communications. Set this variable to &
        &false (or f, both case-insensitive) if you do not want " // methodName // " to finalize the MPI communications for you. &
        &This is a low-level simulation specification variable, relevant to simulations that directly involve MPI parallelism. &
        &If you do not have any MPI-routine calls in your main program, you can safely ignore this variable with its default value. &
        &Note that in non-MPI-enabled simulations, such as serial and Coarray-enabled simulations, the value of &
        &this variable is completely ignored. The default value is " // log2str(MpiFinalizeRequestedObj%def) // "."
    end function constructMpiFinalizeRequested

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(MpiFinalizeRequestedObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(MpiFinalizeRequested_type), intent(in) :: MpiFinalizeRequestedObj
        mpiFinalizeRequested = MpiFinalizeRequestedObj%def
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setMpiFinalizeRequested(MpiFinalizeRequestedObj,mpiFinalizeRequested)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setMpiFinalizeRequested
#endif
        implicit none
        class(MpiFinalizeRequested_type), intent(inout) :: MpiFinalizeRequestedObj
        logical, intent(in)                             :: mpiFinalizeRequested
        MpiFinalizeRequestedObj%val = mpiFinalizeRequested
    end subroutine setMpiFinalizeRequested

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecBase_MpiFinalizeRequested_mod