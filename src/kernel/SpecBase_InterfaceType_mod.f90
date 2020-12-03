!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
!!!!   MIT License
!!!!
!!!!   ParaMonte: plain powerful parallel Monte Carlo library.
!!!!
!!!!   Copyright (C) 2012-present, The Computational Data Science Lab
!!!!
!!!!   This file is part of the ParaMonte library.
!!!!
!!!!   Permission is hereby granted, free of charge, to any person obtaining a 
!!!!   copy of this software and associated documentation files (the "Software"), 
!!!!   to deal in the Software without restriction, including without limitation 
!!!!   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
!!!!   and/or sell copies of the Software, and to permit persons to whom the 
!!!!   Software is furnished to do so, subject to the following conditions:
!!!!
!!!!   The above copyright notice and this permission notice shall be 
!!!!   included in all copies or substantial portions of the Software.
!!!!
!!!!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
!!!!   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
!!!!   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
!!!!   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
!!!!   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
!!!!   OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
!!!!   OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!!!!
!!!!   ACKNOWLEDGMENT
!!!!
!!!!   ParaMonte is an honor-ware and its currency is acknowledgment and citations.
!!!!   As per the ParaMonte library license agreement terms, if you use any parts of 
!!!!   this library for any purposes, kindly acknowledge the use of ParaMonte in your 
!!!!   work (education/research/industry/development/...) by citing the ParaMonte 
!!!!   library as described on this page:
!!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module SpecBase_InterfaceType_mod

    use Constants_mod, only: IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecBase_InterfaceType_mod"
    integer(IK), parameter          :: MAX_INTERFACETYPE_LEN = 511_IK

    character(:)    , allocatable   :: interfaceType ! namelist input

    type                            :: InterfaceType_type
        logical                     :: isC = .false.
        logical                     :: isCPP = .false.
        logical                     :: isFortran = .false.
        logical                     :: isMATLAB = .false.
        logical                     :: isPython = .false.
        logical                     :: isR = .false.
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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructInterfaceType() result(InterfaceTypeObj)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructInterfaceType
#endif
        use Constants_mod, only: NULL_SK
        use Decoration_mod, only: TAB
        use String_mod, only: num2str
        implicit none
        type(InterfaceType_type) :: InterfaceTypeObj
#if defined C_ENABLED
        InterfaceTypeObj%def = "The C Programming Language."
        InterfaceTypeObj%isC = .true.
#elif defined CPP_ENABLED
        InterfaceTypeObj%def = "The C++ Programming Language."
        InterfaceTypeObj%isCPP = .true.
#elif defined FORTRAN_ENABLED
        InterfaceTypeObj%def = "The Fortran Programming Language."
        InterfaceTypeObj%isFortran = .true. ! index(interfaceLowerCase,"fortran") /= 0
#elif defined MATLAB_ENABLED
        InterfaceTypeObj%def = "The MATLAB Programming Language."
        InterfaceTypeObj%isMATLAB = .true.
#elif defined PYTHON_ENABLED
        InterfaceTypeObj%def = "The Python Programming Language."
        InterfaceTypeObj%isPython = .true.
#elif defined R_ENABLED
        InterfaceTypeObj%def = "The R Programming Language."
        InterfaceTypeObj%isR = .true.
#else
#error "Undefined ParaMonte interface type in SpecBase_InterfaceType_mod.f90"
#endif
        if ( allocated(InterfaceTypeObj%null) ) deallocate(InterfaceTypeObj%null)
        allocate( character(len=MAX_INTERFACETYPE_LEN) :: InterfaceTypeObj%null )
        InterfaceTypeObj%null = repeat(NULL_SK, MAX_INTERFACETYPE_LEN)
        InterfaceTypeObj%desc = &
        "This is a ParaMonte internal variable used for providing information about other languages' interface with ParaMonte."
    end function constructInterfaceType

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(InterfaceTypeObj)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(InterfaceType_type), intent(in) :: InterfaceTypeObj
        interfaceType = InterfaceTypeObj%null
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setInterfaceType(InterfaceTypeObj,interfaceType)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setInterfaceType
#endif
        use String_mod, only: getLowerCase
        implicit none
        class(InterfaceType_type), intent(inout)    :: InterfaceTypeObj
        character(*), intent(in)                    :: interfaceType
        if (allocated(InterfaceTypeObj%val)) deallocate(InterfaceTypeObj%val)
        InterfaceTypeObj%val = trim(adjustl(interfaceType))
        if (InterfaceTypeObj%val==trim(adjustl(InterfaceTypeObj%null))) then
            InterfaceTypeObj%val=InterfaceTypeObj%def
        end if
    end subroutine setInterfaceType

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecBase_InterfaceType_mod ! LCOV_EXCL_LINE