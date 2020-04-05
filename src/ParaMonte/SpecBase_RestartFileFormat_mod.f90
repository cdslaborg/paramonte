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

module SpecBase_RestartFileFormat_mod

    use Constants_mod, only: IK
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecBase_RestartFileFormat_mod"

    integer(IK), parameter          :: MAX_LEN_RESTART_FILE_FORMAT = 63

    character(MAX_LEN_RESTART_FILE_FORMAT) :: restartFileFormat

    type                            :: RestartFileFormat_type
        logical                     :: isBinary
        logical                     :: isAscii
        character(6)                :: binary
        character(5)                :: ascii
        character(:), allocatable   :: def
        character(:), allocatable   :: val
        character(:), allocatable   :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setRestartFileFormat, checkForSanity, nullifyNameListVar
    end type RestartFileFormat_type

    interface RestartFileFormat_type
        module procedure            :: constructRestartFileFormat
    end interface RestartFileFormat_type

    private :: constructRestartFileFormat, setRestartFileFormat, checkForSanity, nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructRestartFileFormat(methodName) result(RestartFileFormatObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructRestartFileFormat
#endif
        use Constants_mod, only: NULL_SK, FILE_EXT, FILE_TYPE
        use String_mod, only: num2str
        implicit none
        character(*), intent(in)    :: methodName
        type(RestartFileFormat_type)  :: RestartFileFormatObj

        RestartFileFormatObj%isBinary = .false.
        RestartFileFormatObj%isAscii = .false.
        RestartFileFormatObj%binary = FILE_TYPE%binary
        RestartFileFormatObj%ascii = FILE_TYPE%ascii
        RestartFileFormatObj%def = RestartFileFormatObj%binary

        RestartFileFormatObj%null = repeat(NULL_SK, MAX_LEN_RESTART_FILE_FORMAT)

        RestartFileFormatObj%desc = &
        "restartFileFormat is a string variable that represents the format of the output restart file(s) which are used to restart &
        &an interrupted "// methodName //" simulation. The string value must be enclosed by either single or double quotation &
        &marks when provided as input. Two values are possible:\n\n&
        &    restartFileFormat = '" // RestartFileFormatObj%binary // "'\n\n&
        &            This is the binary file format which is not human-readable, but preserves the exact values of the &
                    &specification variables required for the simulation restart. This full accuracy representation is required &
                    &to exactly reproduce an interrupted simulation. The binary format is also normally the fastest mode of restart file &
                    &generation. Binary restart files will have the " // FILE_EXT%binary // " file extensions.\n\n&
        &    restartFileFormat = '" // RestartFileFormatObj%ascii // "'\n\n&
        &            This is the ASCII (text) file format which is human-readable but does not preserve the full accuracy of &
                    &the specification variables required for the simulation restart. It is also a significantly slower mode of &
                    &restart file generation, compared to the binary format. Therefore, its usage should be limited to situations where &
                    &the user wants to track the dynamics of simulation specifications throughout the simulation time. &
                    &ASCII restart file(s) will have the " // FILE_EXT%ascii //" file extensions.\n\n&
        &The default value is restartFileFormat = '" // RestartFileFormatObj%def // "'. Note that the input values are case-insensitive."
    end function constructRestartFileFormat

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(RestartFileFormatObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        implicit none
        class(RestartFileFormat_type), intent(in) :: RestartFileFormatObj
        restartFileFormat = RestartFileFormatObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setRestartFileFormat(RestartFileFormatObj,restartFileFormat)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setRestartFileFormat
#endif
        use String_mod, only: getLowerCase
        implicit none
        class(RestartFileFormat_type), intent(inout)    :: RestartFileFormatObj
        character(*), intent(in)                        :: restartFileFormat
        RestartFileFormatObj%val = trim(adjustl(restartFileFormat))
        if ( RestartFileFormatObj%val==trim(adjustl(RestartFileFormatObj%null)) ) then
            RestartFileFormatObj%val = trim(adjustl(RestartFileFormatObj%def))
        end if
        if (getLowerCase(RestartFileFormatObj%val)==getLowerCase(RestartFileFormatObj%binary)) RestartFileFormatObj%isBinary = .true.
        if (getLowerCase(RestartFileFormatObj%val)==getLowerCase(RestartFileFormatObj%ascii)) RestartFileFormatObj%isAscii = .true.
    end subroutine setRestartFileFormat

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine checkForSanity(RestartFileFormat,Err,methodName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: checkForSanity
#endif
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        class(RestartFileFormat_type), intent(in)   :: RestartFileFormat
        character(*), intent(in)                    :: methodName
        type(Err_type), intent(inout)               :: Err
        character(*), parameter                     :: PROCEDURE_NAME = "@checkForSanity()"
        if ( .not.(RestartFileFormat%isBinary .or. RestartFileFormat%isAscii) ) then
            Err%occurred = .true.
            Err%msg =   Err%msg // &
                        MODULE_NAME // PROCEDURE_NAME // ": Error occurred. &
                        &The input requested restart file format ('" // RestartFileFormat%val // &
                        "') represented by the variable restartFileFormat cannot be anything other than '" // &
                        RestartFileFormat%binary // "' or '" // RestartFileFormat%ascii // "'. If you don't know an appropriate &
                        &value for RestartFileFormat, drop it from the input list. " // methodName // &
                        " will automatically assign an appropriate value to it.\n\n"
        end if
    end subroutine checkForSanity

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecBase_RestartFileFormat_mod