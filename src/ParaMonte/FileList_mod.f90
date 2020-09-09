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

module FileList_mod

    use String_mod, only: CharVec_type
    use Constants_mod, only: IK, RK
    use Err_mod, only: Err_type
    implicit none

    public

    character(*), parameter :: MODULE_NAME = "@FileList_mod"

    type :: FileList_type
        character(:), allocatable       :: searchStr
        character(:), allocatable       :: orderStr
        character(:), allocatable       :: excludeStr
        integer(IK)                     :: count
        type(CharVec_type), allocatable :: File(:)
        type(Err_type)                  :: Err
    contains
        procedure, nopass               :: get => getFileList
    end type FileList_type

    interface FileList_type
        module procedure :: construct
    end interface FileList_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function construct(searchStr,orderStr,excludeStr,OS) result(FileList)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: construct
#endif
        use System_mod, only: OS_type
        implicit none
        character(*), intent(in), optional  :: searchStr
        character(*), intent(in), optional  :: excludeStr
        character(*), intent(in), optional  :: orderStr
        type(OS_type), intent(in), optional :: OS
        type(FileList_type)                 :: FileList
        if (present(searchStr)) then
            FileList%searchStr = searchStr
        else
            FileList%searchStr = ""
        end if
        if (present(orderStr)) then
            FileList%orderStr = orderStr
        else
            FileList%orderStr = ""
        end if
        if (present(excludeStr)) then
            FileList%excludeStr = excludeStr
        else
            FileList%excludeStr = ""
        end if
        call getFileList(FileList%File,FileList%Err,FileList%count,FileList%searchStr,FileList%orderStr,FileList%excludeStr,OS)
    end function construct

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns a list of files that match searchStr
    ! searchStr can contain the path of the folder of interest to be searched.
    subroutine getFileList(FileList,Err,count,searchStr,orderStr,excludeStr,OS)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getFileList
#endif

        use, intrinsic :: iso_fortran_env, only: output_unit
        use System_mod, only: OS_type, executeCmd, sleep, removeFile
        use Constants_mod, only: IK, RK, MAX_REC_LEN
        use String_mod, only: getLowerCase, num2str
        !use Path_mod, only: winifyPath, linifyPath
        use JaggedArray_mod, only: CharVec_type
        use DateTime_mod, only: DateTime_type
        use Err_mod, only: Err_type

        implicit none

        character(*), intent(in), optional              :: searchStr, orderStr, excludeStr
        type(CharVec_type), allocatable, intent(out)    :: FileList(:)
        type(Err_type), intent(out)                     :: Err
        integer(IK), intent(out), optional              :: count
        type(OS_type), intent(in), optional             :: OS

        character(:), allocatable                       :: search,order,exclude !,searchModified
        character(:), allocatable                       :: command,filename,stdErr,recordTrimmed
        character(MAX_REC_LEN)                          :: record
        integer(IK)                                     :: fileUnit
        integer(IK)                                     :: counter,nrecord,nskip,fileCounter
        logical                                         :: fileIsOpen, OSisWindows

        character(*), parameter                         ::  PROCEDURE_NAME = "@getFileList()"

        Err%occurred = .false.
        Err%msg = ""

        if (present(searchStr)) then
            search = trim(adjustl(searchStr))
        else
            search = ""
        end if

        if (present(orderStr)) then
            if ( trim(adjustl(orderStr))/="" ) then
                order = getLowerCase( trim(adjustl(orderStr) ))
            else
                order = "name"
            end if
        else
            order = "name"
        end if

        ! check if the input order request is supported

        if (order/="name" .and. order/="date") then
            Err%occurred = .true.
            Err%msg =   PROCEDURE_NAME // ": Error Occurred. The requested search order orderStr='" // &
                        order // "' is not supported."
            return
        end if

        if (present(excludeStr)) then
            exclude = trim(adjustl(excludeStr))
        else
            exclude = ""
        end if

        if (present(OS)) then
            OSisWindows = OS%isWindows
        else
            block
                type(OS_type) :: OS
                call OS%query()
                if (OS%Err%occurred) then
                    Err = OS%Err
                    Err%msg =   PROCEDURE_NAME // &
                                ": Error occurred while attempting to query OS type in search of files containing '" // &
                                search // "'.\n" // Err%msg
                    return
                end if
                OSisWindows = OS%isWindows
            end block
        end if

        if (OSisWindows) then

            !call winifyPath(search,searchModified,Err)
            !if (Err%occurred) then
            !    Err%msg =   PROCEDURE_NAME // ": Error occurred while attempting to modify searchStr='" // search // &
            !                "' according to the OS type.\n" // Err%msg
            !    return
            !end if

            ! see: https://www.computerhope.com/dirhlp.htm
            ! /b Uses bare format (no heading information or summary).
            ! /A[:Attributes] List only files with the specified file attributes. Attributes is a series of letters:
            !       D : Directories.
            !       R : Read-only files.
            !       H : Hidden files.
            !       A : Files ready for archiving.
            !       S : System files.
            !       - : Prefix meaning "not".
            ! /O[:SortOrder]	List files in sorted order, indicated by SortOrder:
            !       N : By name (alphabetic).
            !       S : By size (smallest first).
            !       E : By extension (alphabetic).
            !       D : By date and time (earliest first).
            !       G : Group directories first.
            !       - : Prefix to reverse order.
            !       A : By Last Access Date (earliest first).

            if (order=="name") then  ! ascending in name
                command = "dir /b /a-d " // search
            elseif (order=="date") then   ! newest will be first
                command = "dir /b /a-d /o:-d " // search
            end if
            if ( len(exclude)>0 ) command = command // " | findstr /v /i " // exclude

        else    ! It is not windows: either Mac or Linux

            ! Assume Bash environment:
            ! see: http://pubs.opengroup.org/onlinepubs/9699919799/utilities/ls.html
            ! 1 causes each file to be printed on one line
            ! p causes directories to have "/" at the end
            ! t sorts files by modification date, most recent first.
            ! r reverses the sort order (not present here)

            if (order=="name") then  ! ascending in name
                !command = "ls -1 " // searchModified
                command = "ls -1p " // search
            elseif (order=="date") then   ! newest will be first
                command = "ls -1pt " // search
            end if
            if ( len(exclude)>0 ) command = command // " --ignore=" // trim(adjustl(exclude))

        end if

        ! generate a brand new, non-existing filename
        block
            use System_mod, only: RandomFileName_type
            type(RandomFileName_type)   :: RFN
            RFN = RandomFileName_type(key="getFileList")
            if (RFN%Err%occurred) then
                RFN%Err%msg = PROCEDURE_NAME // RFN%Err%msg
                return
            end if
            filename = RFN%path
        end block
        stdErr = filename // ".stderr"

        call executeCmd( command = command//" > "//filename//" 2> "//stdErr, Err=Err )
        if (Err%occurred) then
            Err%msg =   PROCEDURE_NAME // &
                        ": Error occurred while attempting to write the search results to external file.\n" // Err%msg
            return
        end if

        ! now count the number of records in file:

        inquire(file=filename,opened=fileIsOpen,number=fileUnit,iostat=Err%stat)    ! check if the file already exists
        if (Err%stat/=0) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the open status of file='" // filename // "'."
            return
        end if

        if (fileIsOpen) close(fileUnit,iostat=Err%stat)
        if (Err%stat/=0) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred while attempting to close the open file='" // filename // "'."
            return
        end if

        call sleep(seconds=0.1_RK,Err=Err)
        if (Err%occurred) then
            Err%msg = PROCEDURE_NAME // Err%msg
            return
        end if

        open(newunit=fileUnit,file=filename,status="old",iostat=Err%stat)
        if (Err%stat>0) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Unknown error occurred while opening file='" // filename // "'."
            return
        end if

        nskip = 0   ! check filename is not among records
        nrecord = 0 ! number of filenames in the file
        do
            read(fileUnit,"(A)",iostat=Err%stat) record
            if ( is_iostat_eor(Err%stat) ) then
                Err%occurred = .true.
                Err%msg  = PROCEDURE_NAME // ": End-Of-Record error condition occurred while attempting to read &
                         & from file='" // filename // "'."
                return
            elseif ( is_iostat_end(Err%stat) ) then
                exit
            elseif ( Err%stat>0 ) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Unknown error condition occurred while attempting to read &
                        & from file='" // filename // "'."
                return
            else
                recordTrimmed = trim(adjustl(record))
                if(filename==recordTrimmed) nskip = nskip + 1
                nrecord = nrecord + 1
                cycle
            end if
        end do
        close(fileUnit,iostat=Err%stat)
        if (Err%stat/=0) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred while attempting to close the open file='" // filename // "'."
            return
        end if

        allocate(FileList(nrecord-nskip))

        call sleep(seconds=0.1_RK,Err=Err)
        if (Err%occurred) then
            Err%msg = PROCEDURE_NAME // Err%msg
            return
        end if

        open(newunit=fileUnit,file=filename,status="old",iostat=Err%stat)
        if (Err%stat>0) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Unknown error occurred while opening file='" // filename // "'."
            return
        end if

        fileCounter = 0
        do counter = 1,nrecord
            read(fileUnit,"(A)",iostat=Err%stat) record
            if ( is_iostat_eor(Err%stat) ) then
                Err%occurred = .true.
                Err%msg  = PROCEDURE_NAME // ": End-Of-Record error condition occurred while attempting to read &
                         & from file='" // filename // "'."
                return
            elseif ( is_iostat_end(Err%stat) ) then
                exit
            elseif ( Err%stat>0 ) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Unknown error condition occurred while attempting to read &
                        & from file='" // filename // "'."
                return
            end if
            recordTrimmed = trim(adjustl(record))
            if(filename/=recordTrimmed) then
                fileCounter = fileCounter + 1
                FileList(fileCounter)%record = trim(adjustl(record))
            end if
        end do

        if (present(count)) count = fileCounter

        close(fileUnit,iostat=Err%stat)
        if (Err%stat/=0) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred while attempting to close the open file='" // filename // "'."
            return
        end if

        ! remove the files
        call removeFile(filename,OSisWindows,Err)
        if (Err%occurred) then
            Err%msg = PROCEDURE_NAME // Err%msg
            return
        end if
        call removeFile(stdErr,OSisWindows,Err)
        if (Err%occurred) then
            Err%msg = PROCEDURE_NAME // Err%msg
            return
        end if

        !call sleep(seconds=0.1_RK,Err=Err)
        !if (Err%occurred) then
        !    Err%msg = PROCEDURE_NAME // Err%msg
        !    return
        !end if
        !if (OSisWindows) then  ! it is Windows cmd
        !    command = "del " // filename // "; del " // stdErr
        !else
        !    command = "rm " // filename // "; rm " // stdErr
        !end if        
        !call executeCmd( command = command//" > "//filename//" 2> "//stdErr, Err=Err )
        !if (Err%occurred) then
        !    Err%msg =   PROCEDURE_NAME // &
        !                ": Error occurred while attempting to deleting the external file='" // filename // "'.\n" // Err%msg
        !    return
        !end if

    end subroutine getFileList

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end module FileList_mod