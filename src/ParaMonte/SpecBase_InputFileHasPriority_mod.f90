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

module SpecBase_InputFileHasPriority_mod

    implicit none

    ! namelist input
    logical                         :: inputFileHasPriority

    type                            :: InputFileHasPriority_type
        logical                     :: val
        logical                     :: def
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setInputFileHasPriority, nullifyNameListVar
    end type InputFileHasPriority_type

    interface InputFileHasPriority_type
        module procedure            :: constructInputFileHasPriority
    end interface InputFileHasPriority_type

    private :: constructInputFileHasPriority, setInputFileHasPriority

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructInputFileHasPriority(methodName) result(InputFileHasPriorityObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructInputFileHasPriority
#endif
        use String_mod, only: log2str
        implicit none
        character(*), intent(in)        :: methodName
        type(InputFileHasPriority_type) :: InputFileHasPriorityObj
        InputFileHasPriorityObj%def = .false.
        InputFileHasPriorityObj%desc = &
        "If inputFileHasPriority = true (or T, both case-insensitive), then all " // methodName // " variables will be read from &
        &the input file provided by the user, and the parameter specifications from within the programming language &
        &environment (if any are made) will be completely ignored. &
        &If inputFileHasPriority = false (or F, both case-insensitive), then all of " // methodName // " variable values that are &
        &taken from the user-specified input file will be overwritten by their corresponding input values that are set from within &
        &the user's programming environment (if any is provided). Note that this feature is useful when, for example, some " // &
        methodName // " variables have to computed and specified at runtime and therefore, cannot be specified prior to the &
        &program execution. Currently, this functionality (i.e., prioritizing the input file values to input-procedure-argument &
        &values) is available only in the Fortran-interface to the " // methodName // ". &
        &The default value is " // log2str(InputFileHasPriorityObj%def) // "."
    end function constructInputFileHasPriority

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(InputFileHasPriorityObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(InputFileHasPriority_type), intent(in)    :: InputFileHasPriorityObj
        inputFileHasPriority = InputFileHasPriorityObj%def
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setInputFileHasPriority(InputFileHasPriorityObj,inputFileHasPriority)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setInputFileHasPriority
#endif
        implicit none
        class(InputFileHasPriority_type), intent(inout) :: InputFileHasPriorityObj
        logical, intent(in)                             :: inputFileHasPriority
        InputFileHasPriorityObj%val = inputFileHasPriority
    end subroutine setInputFileHasPriority

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecBase_InputFileHasPriority_mod