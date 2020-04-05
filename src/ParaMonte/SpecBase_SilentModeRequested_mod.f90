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

module SpecBase_SilentModeRequested_mod

    implicit none

    ! namelist input
    logical                         :: silentModeRequested

    type                            :: SilentModeRequested_type
        logical                     :: val
        logical                     :: def
        logical                     :: isFalse
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setSilentModeRequested, nullifyNameListVar
    end type SilentModeRequested_type

    interface SilentModeRequested_type
        module procedure            :: constructSilentModeRequested
    end interface SilentModeRequested_type

    private :: constructSilentModeRequested, setSilentModeRequested

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructSilentModeRequested(methodName) result(SilentModeRequestedObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructSilentModeRequested
#endif
        use String_mod, only: log2str
        implicit none
        character(*), intent(in)        :: methodName
        type(SilentModeRequested_type) :: SilentModeRequestedObj
        SilentModeRequestedObj%def = .false.
        SilentModeRequestedObj%isFalse = .true.
        SilentModeRequestedObj%desc = &
        "If silentModeRequested = true (or T, both case-insensitive), then the following contents will not be printed in the &
        &output report file of " // methodName // ":\n\n&
        &    - " // methodName // " interface, compiler, and platform specifications.\n&
        &    - " // methodName // " simulation specification-descriptions.\n\n&
        &The default value is " // log2str(SilentModeRequestedObj%def) // "."
    end function constructSilentModeRequested

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(SilentModeRequestedObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(SilentModeRequested_type), intent(in) :: SilentModeRequestedObj
        silentModeRequested = SilentModeRequestedObj%def
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setSilentModeRequested(SilentModeRequestedObj,silentModeRequested)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setSilentModeRequested
#endif
        implicit none
        class(SilentModeRequested_type), intent(inout)  :: SilentModeRequestedObj
        logical, intent(in)                             :: silentModeRequested
        SilentModeRequestedObj%val = silentModeRequested
        SilentModeRequestedObj%isFalse = .not. silentModeRequested
    end subroutine setSilentModeRequested

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecBase_SilentModeRequested_mod