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

module SpecBase_Description_mod

    use Constants_mod, only: IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecBase_Description_mod"
    integer(IK), parameter          :: MAX_DESCRIPTION_LEN = 4096

    character(:)    , allocatable   :: description ! namelist input

    type                            :: Description_type
        character(:), allocatable   :: val
        character(:), allocatable   :: def
        character(:), allocatable   :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setDescription, nullifyNameListVar
    end type Description_type

    interface Description_type
        module procedure            :: constructDescription
    end interface Description_type

    private :: constructDescription, setDescription, nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructDescription(methodName) result(DescriptionObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructDescription
#endif
        use Constants_mod, only: NULL_SK
        use Decoration_mod, only: TAB
        use String_mod, only: num2str
        implicit none
        type(Description_type)  :: DescriptionObj
        character(*)            :: methodName
        DescriptionObj%def = "Nothing provided by the user."
        if ( allocated(DescriptionObj%null) ) deallocate(DescriptionObj%null)
        allocate( character(len=MAX_DESCRIPTION_LEN) :: DescriptionObj%null )
        DescriptionObj%null = repeat(NULL_SK, MAX_DESCRIPTION_LEN)
        DescriptionObj%desc = &
        "The variable 'description' contains general information about the specific " // methodName // " simulation that &
        &is going to be performed. It has no effects on the simulation and serves only as a general description of the simulation &
        &for future reference. The " // methodName // " parser automatically recognizes the C-style '\\n' escape sequence as the &
        &new-line character, and '\\\\' as the backslash character '\\' if they used in the description. For example, '\\\\n' &
        &will be converted to '\\n' on the output, while '\\n' translates to the new-line character. Other C escape sequences &
        &are neither supported nor needed. The default value for description is '" // DescriptionObj%def // "'."
    end function constructDescription

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(DescriptionObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(Description_type), intent(in) :: DescriptionObj
        description = DescriptionObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setDescription(DescriptionObj,description)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setDescription
#endif
        implicit none
        class(Description_type), intent(inout)  :: DescriptionObj
        character(*), intent(in)                :: description
        if (allocated(DescriptionObj%val)) deallocate(DescriptionObj%val)
        DescriptionObj%val = trim(adjustl(description))
        if (DescriptionObj%val==trim(adjustl(DescriptionObj%null))) DescriptionObj%val=trim(adjustl(DescriptionObj%def))
    end subroutine setDescription

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecBase_Description_mod