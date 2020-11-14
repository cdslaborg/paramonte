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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine nullifyNameListVar(OutputFileNameObj)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: nullifyNameListVar
#endif
        use Path_mod, only: Path_type, MAX_FILE_PATH_LEN
        implicit none
        class(OutputFileName_type), intent(in) :: OutputFileNameObj
        outputFileName = OutputFileNameObj%null
    end subroutine nullifyNameListVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setOutputFileName(OutputFileNameObj,outputFileName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setOutputFileName
#endif
        implicit none
        class(OutputFileName_type), intent(inout)   :: OutputFileNameObj
        character(*)                                :: outputFileName
        OutputFileNameObj%original = trim(adjustl(outputFileName))
        if ( OutputFileNameObj%original==trim(adjustl(OutputFileNameObj%null)) ) then
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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module SpecBase_OutputFileName_mod