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
!!!!       https://github.com/cdslaborg/paramonte/blob/main/ACKNOWLEDGMENT.md
!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief This module contains classes and procedures contents of external files.
!>  \author Amir Shahmoradi

module FileContents_mod

    use JaggedArray_mod, only: CharVec_type
    use Err_mod, only: Err_type
    implicit none

    character(*), parameter                 :: MODULE_NAME = "@FileContents_mod"

    !> The FileContents_type class.
    type :: FileContents_type
        integer                             :: numRecord    !< number f records (lines) in file.
        type(CharVec_type), allocatable     :: Line(:)      !< The list of lines in the file.
        type(Err_type)                      :: Err          !< The error object.
    contains
        procedure, nopass                   :: getNumRecordInFile
        procedure, nopass                   :: getFileContents
    end type FileContents_type

    interface FileContents_type
        module procedure :: constructFileContents
    end interface FileContents_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The constructor of the [FileContents_type](@ref filecontents_type) class.
    !>
    !> @param[in] filePath      :   The path to the file.
    !> @param[in] delEnabled    :   A logical flag indicating whether the file should be deleted upon return (**optional**, default = `.false.`).
    !>
    !> \return
    !> `FileContents` : An object of [FileContents_type](@ref filecontents_type) class.
    function constructFileContents(filePath,delEnabled) result(FileContents)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: constructFileContents
#endif
        implicit none
        character(*), intent(in)            :: filePath
        logical     , intent(in), optional  :: delEnabled
        type(FileContents_type)             :: FileContents
        call getFileContents(filePath,FileContents%Line,FileContents%numRecord,FileContents%Err,delEnabled)
        if (FileContents%Err%occurred) FileContents%Err%msg = "@constructFileContents()" // FileContents%Err%msg
    end function constructFileContents

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Returns the entire content of a file as an array of strings.
    !>
    !> @param[in]   path        :   The path to the file.
    !> @param[out]  Contents    :   A list of lines in the file. Each array element corresponds to one line (record) in the file.
    !> @param[out]  numRecord   :   The number of lines in the file.
    !> @param[out]  Err         :   An object of [Err_type](@ref err_mod::err_type) indicating whether error has occurred during the file IO.
    !> @param[out]  delEnabled  :   An optional logical value indicating whether the file should be deleted upon successful reading of it.
    subroutine getFileContents(path, Contents, numRecord, Err, delEnabled)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getFileContents
#endif
        use JaggedArray_mod, only: CharVec_type
        use String_mod, only: num2str
        use Err_mod, only: Err_type
        implicit none
        character(len=*), intent(in)                    :: path
        type(CharVec_type), allocatable, intent(out)    :: Contents(:)
        logical, intent(in), optional                   :: delEnabled
        type(Err_type), intent(out)                     :: Err
        integer, intent(out)                            :: numRecord

        character(99999)                                :: record
        character(:), allocatable                       :: closeStatus
        integer                                         :: fileUnit, irecord, i
        logical                                         :: delEnabledDefault, isOpen, fileExists

        character(*), parameter                         :: PROCEDURE_NAME = MODULE_NAME//"@getFileContents()"

        Err%occurred = .false.
        Err%msg = ""

        delEnabledDefault = .false.
        if (present(delEnabled)) then
            delEnabledDefault = delEnabled
        end if

        call getNumRecordInFile(path,numRecord,Err)
        if (Err%occurred) then
        ! LCOV_EXCL_START
            Err%msg = PROCEDURE_NAME // Err%msg
            return
        end if
        ! LCOV_EXCL_STOP

        allocate(Contents(numRecord))

        open( newunit = fileUnit & ! LCOV_EXCL_LINE
            , status = "old" & ! LCOV_EXCL_LINE
            , file = path & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        do irecord = 1,numRecord
            read(fileUnit,"(A)", iostat=Err%stat) record
            if (Err%stat==0) then
                Contents(irecord)%record = trim(adjustl(record))
            ! LCOV_EXCL_START
            elseif (is_iostat_end(Err%stat)) then
                Err%occurred = .true.
                Err%msg =   PROCEDURE_NAME // ": End-of-file error occurred while expecting " // &
                            num2str(numRecord-irecord+1) // " records in file='" // path // "'."
                return
            elseif (is_iostat_eor(Err%stat)) then
                Err%occurred = .true.
                Err%msg =   PROCEDURE_NAME // ": End-of-record error occurred while reading line number " // &
                            num2str(irecord) // " from file='" // path // "'."
                return
            else
                Err%occurred = .true.
                Err%msg =   PROCEDURE_NAME // ": Unknown error occurred while reading line number " // &
                            num2str(irecord) // " from file='" // path // "'."
                return
            end if
            ! LCOV_EXCL_STOP
        end do

        ! attempt to delete the file repeatedly. This is important on windows systems as the file often remains locked.

        blockDeleteFile: if (delEnabledDefault) then

            ! Attempt to delete the file repeatedly

            loopDeleteFile: do i = 1, 100
                inquire(file=path, opened=isOpen, exist=fileExists, iostat = Err%stat)
                if (Err%stat==0) then
                    if (fileExists) then
                        if (.not. isOpen) then
                            open( newunit = fileUnit    & ! LCOV_EXCL_LINE
                                , status = "replace"    & ! LCOV_EXCL_LINE
                                , iostat = Err%stat     & ! LCOV_EXCL_LINE
                                , file = path           & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
                                , SHARED                & ! LCOV_EXCL_LINE
#endif
                                )
                        end if
                        close(fileUnit, status="delete", iostat = Err%stat) ! parallel processes cannot delete the same file. Err%stat is required to handle exceptions.
                    else
                        return
                    end if
                end if
            end do loopDeleteFile

        end if blockDeleteFile

        ! If file deletion fails or if it should not be deleted, then simply close the file and return.

        close(fileUnit, status = "keep", iostat = Err%stat) ! parallel processes cannot delete the same file. Err%stat is required to handle exceptions.

        !if (Err%stat>0) then
        !! LCOV_EXCL_START
        !    Err%occurred = .true.
        !    Err%msg = PROCEDURE_NAME // "Error occurred while attempting to close or delete the open file='" // path // "'."
        !    return
        !end if
        !! LCOV_EXCL_STOP

    end subroutine getFileContents

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Returns the number of lines in a file.
    !>
    !> @param[in]   filePath    :   The path to the file.
    !> @param[out]  numRecord   :   The number of lines in the file.
    !> @param[out]  Err         :   An object of [Err_type](@ref err_mod::err_type) indicating whether error has occurred during the file IO.
    !> @param[in]   exclude     :   A string. If any line matches `exclude`, it will NOT be counted.
    subroutine getNumRecordInFile(filePath,numRecord,Err,exclude)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getNumRecordInFile
#endif
        use Constants_mod, only: IK
        implicit none
        character(len=*), intent(in)                :: filePath
        integer(IK)     , intent(out)               :: numRecord
        type(Err_type)  , intent(out)               :: Err
        character(*)    , intent(in)    , optional  :: exclude
        character(len=1)                            :: record
        integer                                     :: fileUnit
        logical                                     :: fileExists, fileIsOpen, excludeIsPresent
        integer                                     :: iostat

        character(*), parameter                     :: PROCEDURE_NAME = "@getNumRecordInFile()"

        Err%occurred = .false.
        Err%msg = ""
        excludeIsPresent = present(exclude)

        ! Check if file exists
        ! GFortran 7.3 bug: If file is not open, compiler assumes internal file if `number` is specified, causing runtime error.

        inquire( file=filePath, exist=fileExists, opened=fileIsOpen, number=fileUnit, iostat=Err%stat )
        if (Err%stat/=0) then
        ! LCOV_EXCL_START
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file='" // filePath // "'."
            return
        end if
        ! LCOV_EXCL_STOP

        if (.not.fileExists) then
        ! LCOV_EXCL_START
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": The input file='" // filePath // "' does not exist."
            return
        end if
        ! LCOV_EXCL_STOP

        if (fileIsOpen) close(unit=fileUnit,iostat=Err%stat)
        if (Err%stat>0) then
        ! LCOV_EXCL_START
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred while attempting to close the open input file='" // filePath // "'."
            return
        end if
        ! LCOV_EXCL_STOP

        open( newunit = fileUnit & ! LCOV_EXCL_LINE
            , iostat = Err%stat & ! LCOV_EXCL_LINE
            , file = filePath & ! LCOV_EXCL_LINE
            , status = "old" & ! LCOV_EXCL_LINE
#if defined INTEL_COMPILER_ENABLED && defined OS_IS_WINDOWS
            , SHARED & ! LCOV_EXCL_LINE
#endif
            )
        if (Err%stat>0) then
        ! LCOV_EXCL_START
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred while opening input file='" // filePath // "'."
            return
        end if
        ! LCOV_EXCL_STOP

        numRecord = 0_IK
        do
            read(fileUnit,'(A)',iostat=iostat) record
            if(iostat==0) then
                if (excludeIsPresent) then
                    if (trim(adjustl(record))/=exclude) numRecord = numRecord + 1_IK
                else
                    numRecord = numRecord + 1_IK
                end if
                cycle
            elseif(is_iostat_end(iostat) .or. is_iostat_eor(iostat)) then
                exit
            ! LCOV_EXCL_START
            else
                Err%occurred = .true.
                Err%stat = iostat
                Err%msg = PROCEDURE_NAME // ": Error occurred while reading input file='" // filePath // "'."
                return
            end if
            ! LCOV_EXCL_STOP
        end do

        close(fileUnit,iostat=Err%stat)
        if (Err%stat>0) then
        ! LCOV_EXCL_START
            Err%occurred = .true.
            Err%msg =   PROCEDURE_NAME // ": Error occurred while attempting to close the open input file='" // &
                        filePath // "' after counting the number of records in file."
            return
        end if
        ! LCOV_EXCL_STOP

    end subroutine getNumRecordInFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module FileContents_mod ! LCOV_EXCL_LINE