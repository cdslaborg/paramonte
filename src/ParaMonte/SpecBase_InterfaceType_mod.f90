!***********************************************************************************************************************************
!***********************************************************************************************************************************
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
!***********************************************************************************************************************************
!***********************************************************************************************************************************

module SpecBase_InterfaceType_mod

    use Constants_mod, only: IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecBase_InterfaceType_mod"
    integer(IK), parameter          :: MAX_INTERFACETYPE_LEN = 511_IK

    character(:)    , allocatable   :: interfaceType ! namelist input

    type                            :: InterfaceType_type
        logical                     :: isFortran = .false.
        logical                     :: isMATLAB = .false.
        logical                     :: isPython = .false.
        logical                     :: isClang = .false.
        character(:), allocatable   :: val
        character(:), allocatable   :: def
        character(:), allocatable   :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setInterfaceType, nullifyNameListVar
    end type InterfaceType_type

    interface InterfaceType_type
        module procedure            :: constructInterfaceType
    end interface InterfaceType_type

    private :: constructInterfaceType, setInterfaceType, nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructInterfaceType() result(InterfaceTypeObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructInterfaceType
#endif
        use Constants_mod, only: NULL_SK
        use Decoration_mod, only: TAB
        use String_mod, only: num2str
        implicit none
        type(InterfaceType_type) :: InterfaceTypeObj
#if defined C_ENABLED
        InterfaceTypeObj%def = "The C Programming Language."
#elif defined FORTRAN_ENABLED
        InterfaceTypeObj%def = "The Fortran Programming Language."
#elif defined MATLAB_ENABLED
        InterfaceTypeObj%def = "The MATLAB Programming Language."
#elif defined PYTHON_ENABLED
        InterfaceTypeObj%def = "The Python Programming Language."
#else
#error "Undefined ParaMonte interface type in SpecBase_InterfaceType_mod.f90"
#endif
        if ( allocated(InterfaceTypeObj%null) ) deallocate(InterfaceTypeObj%null)
        allocate( character(len=MAX_INTERFACETYPE_LEN) :: InterfaceTypeObj%null )
        InterfaceTypeObj%null = repeat(NULL_SK, MAX_INTERFACETYPE_LEN)
        InterfaceTypeObj%desc = &
        "This is a ParaMonte internal variable used for providing information about other languages' interface with ParaMonte."
    end function constructInterfaceType

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(InterfaceTypeObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(InterfaceType_type), intent(in) :: InterfaceTypeObj
        interfaceType = InterfaceTypeObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setInterfaceType(InterfaceTypeObj,interfaceType)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterfaceType
#endif
        use String_mod, only: getLowerCase
        implicit none
        class(InterfaceType_type), intent(inout)    :: InterfaceTypeObj
        character(*), intent(in)                    :: interfaceType
        character(:), allocatable                   :: interfaceLowerCase
        if (allocated(InterfaceTypeObj%val)) deallocate(InterfaceTypeObj%val)
        InterfaceTypeObj%val = trim(adjustl(interfaceType))
        if (InterfaceTypeObj%val==trim(adjustl(InterfaceTypeObj%null))) then
            InterfaceTypeObj%val=InterfaceTypeObj%def
        end if
        interfaceLowerCase = getLowerCase(InterfaceTypeObj%val)
        if (index(interfaceLowerCase,"fortran")/=0) then
            InterfaceTypeObj%isFortran = .true.
        elseif (index(interfaceLowerCase,"matlab")/=0) then
            InterfaceTypeObj%isMATLAB = .true.
        elseif (index(interfaceLowerCase,"python")/=0) then
            InterfaceTypeObj%isPython = .true.
        else
            InterfaceTypeObj%isClang = .true.
        end if
    end subroutine setInterfaceType

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecBase_InterfaceType_mod