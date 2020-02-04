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

module SpecBase_OutputFileName_mod

    use Path_mod, only: Path_type
    implicit none

    character(*), parameter         :: MODULE_NAME = "@SpecBase_OutputFileName_mod"

    character(:), allocatable       :: outputFileName ! namelist input

    type, extends(Path_type)        :: OutputFileName_type
        character(:), allocatable   :: def
        character(:), allocatable   :: namePrefix
        character(:), allocatable   :: pathPrefix
        character(:), allocatable   :: null
        character(:), allocatable   :: desc
    contains
        procedure, pass             :: set => setOutputFileName, nullifyNameListVar
    end type OutputFileName_type

    interface OutputFileName_type
        module procedure            :: constructOutputFileName
    end interface OutputFileName_type

    private :: constructOutputFileName, setOutputFileName, nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    function constructOutputFileName(methodName) result(OutputFileNameObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructOutputFileName
#endif

        use Path_mod, only: MAX_FILE_PATH_LEN
        use Constants_mod, only: NULL_SK
        use Decoration_mod, only: TAB

        implicit none

        character(*), intent(in)    :: methodName
        type(OutputFileName_type)   :: OutputFileNameObj
        character(8)                :: date
        character(10)               :: time

        call date_and_time(date,time)
        OutputFileNameObj%def = methodName // "_run_" // date // "_" // time(1:6) // "_" // time(8:10)

        OutputFileNameObj%null = repeat(NULL_SK, MAX_FILE_PATH_LEN)
        OutputFileNameObj%desc = &
        "outputFileName contains the path and the base of the filename for " // methodName // " output files. &
        &If not provided by the user, the default outputFileName is constructed from the current date and time:\n\n" &
        // TAB // methodName // "_run_yyyymmdd_hhmmss_mmm\n\n&
        &where yyyy, mm, dd, hh, mm, ss, mmm stand respectively for the current year, month, day, hour, minute, second, &
        &and millisecond. In such a case, the default directory for the output files will be the current working directory of " &
        // methodName // ". If outputFileName is provided, but ends with a separator character '/' or '\' (as in Linux or Windows OS), &
        &then its value will be used as the directory to which " // methodName // " output files will be written. In this case, &
        &the output file naming convention described above will be used. Also, the given directory will be automatically created &
        &if it does not exist already."

    end function constructOutputFileName

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine nullifyNameListVar(OutputFileNameObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use Path_mod, only: Path_type, MAX_FILE_PATH_LEN
        implicit none
        class(OutputFileName_type), intent(in) :: OutputFileNameObj
        outputFileName = OutputFileNameObj%null
    end subroutine nullifyNameListVar

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setOutputFileName(OutputFileNameObj,outputFileName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setOutputFileName
#endif
        implicit none
        class(OutputFileName_type), intent(inout)   :: OutputFileNameObj
        character(*)                                :: outputFileName
        OutputFileNameObj%original = trim(adjustl(outputFileName))
        if ( trim(adjustl(OutputFileNameObj%original))==trim(adjustl(OutputFileNameObj%null)) ) then
            OutputFileNameObj%original = OutputFileNameObj%def
        end if

        ! set the outputFileName the same on all images. This becomes relevant in two scenarios:
        ! 1. when outputFileName is missing as input.
        ! 2. when outputFileName is present, but is not the same on all images (for example when called from Python)

#if defined CAF_ENABLED
        block
            character(63), save :: co_defaultOutputFileName[*]
            if (this_image()==1) then
                co_defaultOutputFileName = OutputFileNameObj%def
                sync images(*)
            else
                sync images(1)
                OutputFileNameObj%def = trim(adjustl(co_defaultOutputFileName[1]))
            end if
        end block
#elif defined MPI_ENABLED
        block
            use mpi
            integer :: ierrMPI
            character(63) :: co_defaultOutputFileName
            co_defaultOutputFileName = OutputFileNameObj%def
            ! bcast co_defaultOutputFileName from image one to all others
            call mpi_bcast  ( co_defaultOutputFileName  & ! buffer
                            , 63                        & ! count
                            , mpi_character             & ! datatype
                            , 0                         & ! root
                            , mpi_comm_world            & ! comm
                            , ierrMPI                   & ! ierr
                            )
            OutputFileNameObj%def = trim(adjustl(co_defaultOutputFileName))
        end block
#endif

    end subroutine setOutputFileName

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module SpecBase_OutputFileName_mod