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

module SpecBase_OutputColumnWidth_mod

    use Constants_mod, only: IK
    implicit none

    integer(IK), parameter          :: OUTPUT_COL_WIDTH = 63

    integer(IK)                     :: outputColumnWidth ! namelist input

    character(*), parameter         :: MODULE_NAME = "@SpecBase_OutputColumnWidth_mod"

    type                            :: OutputColumnWidth_type
        integer(IK)                 :: val
        integer(IK)                 :: def
        integer(IK)                 :: null
        character(:), allocatable   :: str
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setOutputColumnWidth, checkForSanity, nullifyNameListVar
    end type OutputColumnWidth_type

    interface OutputColumnWidth_type
        module procedure                :: constructOutputColumnWidth
    end interface OutputColumnWidth_type

    private :: constructOutputColumnWidth, setOutputColumnWidth

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructOutputColumnWidth(methodName) result(OutputColumnWidthObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructOutputColumnWidth
#endif
        use Constants_mod, only: IK, NULL_IK
        use String_mod, only: num2str
        implicit none
        character(*), intent(in)  :: methodName
        type(OutputColumnWidth_type) :: OutputColumnWidthObj
        OutputColumnWidthObj%def = 0_IK
        OutputColumnWidthObj%null = NULL_IK
        OutputColumnWidthObj%desc = &
        "The variable outputColumnWidth is a non-negative integer number that determines the width of &
        &the data columns in " // methodName // " formatted output files that have tabular structure. &
        &If it is set to zero, " // methodName // " will ensure to set the width of each output element &
        &to the minimum possible width without losing the requested output precision. In other words, &
        &setting outputColumnWidth = 0 will result in the smallest-size for the formatted output files that are in ASCII format. &
        &The default value is " // num2str(OutputColumnWidthObj%def) // "."
    end function constructOutputColumnWidth

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(OutputColumnWidthObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(OutputColumnWidth_type), intent(in) :: OutputColumnWidthObj
        outputColumnWidth = OutputColumnWidthObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setOutputColumnWidth(OutputColumnWidthObj,outputColumnWidth)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setOutputColumnWidth
#endif
        use String_mod, only: num2str
        use Constants_mod, only: IK
        implicit none
        class(OutputColumnWidth_type), intent(inout)  :: OutputColumnWidthObj
        integer(IK), intent(in)                    :: outputColumnWidth
        OutputColumnWidthObj%val = outputColumnWidth
        if (OutputColumnWidthObj%val==OutputColumnWidthObj%null) then
            OutputColumnWidthObj%val = OutputColumnWidthObj%def
        end if
        OutputColumnWidthObj%str = num2str(OutputColumnWidthObj%val)
    end subroutine setOutputColumnWidth

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine checkForSanity(OutputColumnWidthObj,Err,methodName,outputRealPrecision)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Err_mod, only: Err_type
        use Constants_mod, only: IK
        use String_mod, only: num2str
        implicit none
        class(OutputColumnWidth_type), intent(in) :: OutputColumnWidthObj
        character(*), intent(in)               :: methodName
        integer(IK) , intent(in)               :: outputRealPrecision
        type(Err_type), intent(inout)          :: Err
        character(*), parameter                :: PROCEDURE_NAME = "@checkForSanity()"
        if ( OutputColumnWidthObj%val<0_IK ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input value for variable outputColumnWidth must be a non-negative integer. &
                        &If you are not sure about the appropriate value for this variable, simply drop it from the input. " // &
                        methodName // " will automatically assign an appropriate value to it.\n\n"
        elseif ( OutputColumnWidthObj%val>0_IK .and. OutputColumnWidthObj%val<outputRealPrecision+7_IK ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. The input value for variable outputColumnWidth &
                        &must be equal to or greater than the input value for outputRealPrecision + 7. &
                        &If you are not sure about the appropriate value for this variable, either set it to zero on input, &
                        &or simply drop it from the input. " // methodName // &
                        " will automatically assign an appropriate value to it.\n\n"
        end if
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecBase_OutputColumnWidth_mod