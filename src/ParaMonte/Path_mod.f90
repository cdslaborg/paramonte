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

module Path_mod

    use Constants_mod, only: IK
    use Err_mod, only: Err_type
    implicit none

    character(*), parameter :: MODULE_NAME = "@Path_mod"

    integer(IK), parameter  :: MAX_FILE_PATH_LEN = 2047

    ! Windows reserved characters (not allowed in filenames):
    character(*), parameter :: WINDOWS_RESERVED_CHAR = "<>:" // '"' // "|?*" ! /\

#if defined IFORT_ENABLED

    character(*), parameter :: SHELL_ESCAPE_CHAR =  &
                                                    " " // & ! — space character
                                                    "!" // & ! — history expansion.
                                                    '"' // & ! — shell syntax.
                                                    "#" // & ! — comment start when preceded by whitespace; zsh wildcards.
                                                    "$" // & ! — shell syntax.
                                                    "&" // & ! — shell syntax.
                                                    "'" // & ! — shell syntax.
                                                    "(" // & ! — even in the middle of a word: ksh extended globs (also available in bash and zsh); zsh wildcards.
                                                    ")" // & ! - even in the middle of a word: ksh extended globs (also available in bash and zsh); zsh wildcards.
                                                    "*" // & ! — sh wildcard.
                                                    "," // & ! — only inside brace expansion.
                                                    ";" // & ! — shell syntax.
                                                    "<" // & ! — shell syntax.
                                                    "=" // & ! — in zsh, when it is at the beginning of a file name (filename expansion with PATH lookup).
                                                    ">" // & ! — shell syntax.
                                                    "?" // & ! — sh wildcard.
                                                    "[" // & ! — sh wildcard.
                                                    "\" // & ! — shell syntax.
                                                    "]" // & ! — you may get away with leaving it unquoted.
                                                    "^" // & ! — history expansion; zsh wildcard.
                                                    "`" // & ! — shell syntax.
                                                    "{" // & ! — brace expansion.
                                                    "|" // & ! — shell syntax.
                                                    "}" // & ! — needs to be escaped in zsh, other shells are more lenient when there is no matching opening brace.
                                                    "~"      ! — home directory expansion when at the beginning of a filename; zsh wildcard; safe when it is the last character.

#else

    ! stupid gfortran gives error on the above syntax
    character(*), parameter :: SHELL_ESCAPE_CHAR = " !"//'"#$&'//"'()*,;<=>?[\]^`{|}~"

#endif

    type :: Path_type
        character(:), allocatable       :: original, modified, dir, name, base, ext
        character(1)                    :: slashOS    ! the type of slash with which the original path was "modified"
        type(Err_type)                  :: Err
    contains
        procedure, pass                 :: query => queryPath
        procedure, nopass               :: modify => modifyPath, getSlashOS
        procedure, nopass               :: getDirNameExt, getDirFullName, getNameExt
        procedure, nopass               :: winify => winifyPath, linify => linifyPath
        procedure, nopass               :: mkdir
    end type Path_type

    interface Path_type
        module procedure :: constructPath
    end interface Path_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructPath(inputPath,OS) result(Path)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructPath
#endif
        use System_mod, only: OS_type
        implicit none
        type(Path_type)                     :: Path
        character(*), intent(in)            :: inputPath
        type(OS_type), intent(in), optional :: OS
        call Path%query(inputPath,OS)
    end function constructPath

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Construct a Path object as output. On output, check Err%occurred before using the output Path.
    subroutine queryPath(Path,inputPath,OS)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: queryPath
#endif
        use Err_mod, only: Err_type
        use Constants_mod, only: IK
        use System_mod, only: OS_type
        use String_mod, only: replaceStr
        implicit none
        class(Path_type), intent(inout)     :: Path
        character(*), intent(in), optional  :: inputPath
        type(OS_type), intent(in), optional :: OS
        logical                             :: OSisWindows

        character(*), parameter             :: PROCEDURE_NAME = "@queryPath()"

        Path%Err%occurred = .false.
        Path%Err%msg = ""

        if (present(inputPath)) then
            Path%original = trim(adjustl(inputPath))
        elseif (.not.allocated(Path%original)) then
            Path%Err%occurred = .true.
            Path%Err%msg = PROCEDURE_NAME // ": Error occurred. Neither inputPath argument is given as input,&
                                             & nor Path%original is allocated to construct the Path object."
            return
        else
            if ( len(trim(adjustl(Path%original)))==0 ) then
                Path%Err%occurred = .true.
                Path%Err%msg = PROCEDURE_NAME // ": Error occurred. Neither inputPath argument is given as input,&
                                                 & nor Path%original has a non-blank length > 0 to construct the Path object."
                return
            end if
        end if

        if (present(OS)) then
            Path%slashOS = OS%slash
            OSisWindows = OS%isWindows
        else
            block
                type(OS_type) :: OS
                call OS%query()
                if (OS%Err%occurred) then
                    Path%Err%stat = OS%Err%stat
                    Path%Err%occurred = OS%Err%occurred
                    Path%Err%msg = PROCEDURE_NAME // ": Error occurred while querying OS type.\n" // Path%Err%msg
                end if
                Path%slashOS = OS%slash
                OSisWindows = OS%isWindows
            end block
            if (Path%Err%occurred) return
        end if

        if (OSisWindows) then
            call winifyPath(Path%original,Path%modified,Path%Err)
            if (Path%Err%occurred) then
                Path%Err%msg = PROCEDURE_NAME//": Error occurred while making path='"//Path%original//"' compatible with Windows OS.\n"//Path%Err%msg
                return
            end if
        else
            ! if the path contains both / and \, then assume that it is already in linux style
            if (index(Path%original,"/")==0) then ! path is given in Windows style
                call linifyPath(Path%original,Path%modified)
            else
                Path%modified = Path%original
            end if
        end if

        call Path%getDirNameExt( Path%modified, Path%slashOS, Path%dir, Path%name, Path%ext )
        Path%base = Path%dir // Path%name

    end subroutine queryPath

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! This code assumes that the input path is a Linux path. Windows paths like .\(paramonte)\paramonte.nml will be horribly
    ! treated by this routine as \( also represents a Linux escape character. The result will be .(paramonte)\paramonte.nml
    ! this routine strictly assumes that there is no dangling \ in the input Linux path, and if there is, then either it is used
    ! to escape the special shell characters, or otherwise, the path is a Windows path.
    pure subroutine winifyPath(inputPath,outputPath,Err)!,ignoreWindowsReservedChars)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: winifyPath
#endif
        use Err_mod, only: Err_type
        use Constants_mod, only: IK
        use String_mod, only: replaceStr
        implicit none
        character(len=*), intent(in)            :: inputPath
        character(:), allocatable, intent(out)  :: outputPath
        type(Err_type), intent(out)             :: Err
       !logical, intent(in), optional           :: ignoreWindowsReservedChars
       !logical                                 :: reservedCharInspectionNeeded
        character(:), allocatable               :: outputPathNew
        integer(IK)                             :: i, j

        character(*), parameter                 :: PROCEDURE_NAME = "@winifyPath()"

        Err%occurred = .false.
        Err%msg = ""

        ! check if any character in the input path is Windows Reserved Character:

        !reservedCharInspectionNeeded = .true.
        !if (present(ignoreWindowsReservedChars)) reservedCharInspectionNeeded = .not. ignoreReservedChars
        !if (reservedCharInspectionNeeded) then
        !    do i = 1, len(WINDOWS_RESERVED_CHAR)
        !       if ( index(inputPath,WINDOWS_RESERVED_CHAR(i:i)) /= 0 ) then
        !           Err%occurred = .true.
        !           Err%msg =   PROCEDURE_NAME // ": Error occurred. Invalid Windows character '" // &
        !                       WINDOWS_RESERVED_CHAR(i:i) // "' detected in the input file path='" // inputPath // "'."
        !           return
        !       end if
        !    end do
        !end if

        !if ( index(inputPath,"\\") /= 0 ) then
        !    Err%occurred = .true.
        !    Err%msg = PROCEDURE_NAME // ": Error occurred. Invalid Windows character '\' corresponding to '\\' detected &
        !            & in the input file path='" // inputPath // "'."
        !    return
        !end if

        ! note that multiple \ character in sequence is meaningless in Linux (basically \\ reduces to \), 
        ! and in Windows means the same as a single \. Therefore, reduce all sequential \ characters to a single \.
        outputPath = trim(adjustl(inputPath))
        loopRemoveMultipleSlash: do
            outputPathNew = replaceStr(outputPath,"\\","\")
            if (outputPathNew==outputPath) exit loopRemoveMultipleSlash
            outputPath = outputPathNew
        end do loopRemoveMultipleSlash

        ! Now check for the presence of any Linux Shell Escape Character in the input path without a preceding \.
        ! If there is any, this would imply that the input path is a Windows path,
        ! otherwise a escape character without preceding \ would be invalid in Linux.
        do i = 1, len(SHELL_ESCAPE_CHAR)
            if (SHELL_ESCAPE_CHAR(i:i)/="\") then
                do j = 1, len(outputPath)
                    if (outputPath(j:j)==SHELL_ESCAPE_CHAR(i:i)) then
                        if (j==1) then
                            !write(*,*) "secLoc = 1"
                            return  ! it is a windows path, so no need for further winifying
                        elseif (outputPath(j-1:j-1)/="\") then
                            !write(*,*) "secLoc= ", j-1
                            return  ! it is a windows path, so no need for further winifying
                        end if
                    end if
                end do
            end if
        end do

        ! By now, there is no way but to assume that the path is indeed a Linux path
        ! Thus, correct for any Linux Shell Escape Character in the input path:
        do i = 1, len(SHELL_ESCAPE_CHAR)
            outputPath = replaceStr(outputPath,"\"//SHELL_ESCAPE_CHAR(i:i),SHELL_ESCAPE_CHAR(i:i))
        end do

        ! Now remove any remaining backslash in the input path:
        ! commented out: it is assumed that there are no dangling \ in the Linux path
        !outputPath = replaceStr(outputPath,"\","")

        ! check if the file name contains white space. if so, put the entire name in quotations
        if ( index(outputPath," ") /= 0 ) then
            outputPath = '"' // outputPath  // '"'
        end if
        outputPath = replaceStr(outputPath,"/","\")

    end subroutine winifyPath

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    pure subroutine linifyPath(inputPath,outputPath)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: linifyPath
#endif
        use Constants_mod, only: IK
        use String_mod, only: replaceStr
        implicit none
        character(*), intent(in)                :: inputPath
        character(:), allocatable, intent(out)  :: outputPath
        character(:), allocatable               :: outputPathNew
        integer(IK)                             :: i

        character(*), parameter                 :: PROCEDURE_NAME = "@linifyPath()"

        ! check if the path is sandwiched between quotation marks. If so, remove them:
        outputPath = trim(adjustl(inputPath))
        i = len(outputPath)
        if (i==0) return
        if ( i>1 ) then
            if ( (outputPath(1:1)=='"' .and. outputPath(i:i)=='"') .or. (outputPath(1:1)=="'" .and. outputPath(i:i)=="'") )then
                outputPathNew = outputPath(2:i-1)
            else
                outputPathNew = outputPath
            end if
        end if

        ! First change all backslashes to forward slash:
        outputPath = replaceStr(outputPathNew,"\","/")

        ! Now correct for any Linux Shell Escape Character in the input path:
        do i = 1, len(SHELL_ESCAPE_CHAR)
            if (SHELL_ESCAPE_CHAR(i:i)/="\") then
                outputPathNew = replaceStr(outputPath,SHELL_ESCAPE_CHAR(i:i),"\"//SHELL_ESCAPE_CHAR(i:i))
                outputPath = outputPathNew
            end if
        end do

        !! Now correct for any white spaces in outputPath:
        !outputPath = replaceStr(outputPath," ","\ ")

    end subroutine linifyPath

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Returns the slash type of the OS (backslash in Windows, forward-slash in NonWindows).
    subroutine getSlashOS(slashOS,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getSlashOS
#endif
        use Err_mod, only: Err_type
        use System_mod, only: OS_type
        use String_mod, only: replaceStr
        implicit none
        character(1)  , intent(out) :: slashOS
        type(Err_type), intent(out) :: Err

        type(OS_type)               :: OS
        character(*), parameter     :: PROCEDURE_NAME = "@getSlashOS()"

        Err%occurred = .false.
        Err%msg = ""

        call OS%query()

        if (OS%Err%occurred) then
            Err = OS%Err
            Err%msg = PROCEDURE_NAME // ": Error occurred while fetching the OS slash character.\n" // Err%msg
            return
        end if

        if (OS%isWindows) then
            slashOS = "\"
        else
            slashOS = "/"
        end if

    end subroutine getSlashOS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine modifyPath(inputPath,outputPath,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: modifyPath
#endif
        use Err_mod, only: Err_type
        use Constants_mod, only: IK
        use System_mod, only: OS_type
        use String_mod, only: replaceStr
        implicit none
        character(len=*), intent(in) :: inputPath
        character(:), allocatable, intent(out) :: outputPath
        type(Err_type), intent(out) :: Err

        type(OS_type) :: OS

        character(*), parameter :: PROCEDURE_NAME = "@modifyPath()"

        outputPath = trim(adjustl(inputPath))

        Err%occurred = .false.
        Err%msg = ""

        call OS%query()

        if (OS%Err%occurred) then
            Err = OS%Err
            Err%msg = PROCEDURE_NAME // ": Error occurred while modifying inputPath='" // outputPath // "'.\n" // Err%msg
            return
        end if

        if (OS%isWindows) then
            call winifyPath(inputPath,outputPath,Err)
            if (Err%occurred) then
                Err%msg =  PROCEDURE_NAME // ": Error occurred while making path='" // &
                           inputPath // "' compatible with Windows OS.\n" // Err%msg
                return
            end if
        else
            call linifyPath(inputPath,outputPath)
        end if

    end subroutine modifyPath

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine getDirNameExt(path,slash,dir,name,ext)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirNameExt
#endif
        implicit none
        character(*)             , intent(in)   :: path
        character(1)             , intent(in)   :: slash
        character(:), allocatable, intent(out)  :: dir,name,ext
        character(:), allocatable               :: fullName
        call getDirFullName(path,slash,dir,fullName)
        call getNameExt(fullName,name,ext)
    end subroutine getDirNameExt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns the directory and full filename of the input path.
    subroutine getDirFullName(path,slash,dir,fullName)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirFullName
#endif
        use Constants_mod, only: IK
        implicit none
        character(*)             , intent(in)   :: path
        character(1)             , intent(in)   :: slash
        character(:), allocatable, intent(out)  :: dir,fullName

        integer(IK)                             :: pathLen, slashPos

        pathLen = len(path)

        if ( pathLen==0 ) then
            dir=""; fullName=""
            return
        end if

        slashPos = index(path,slash,back=.true.)

        if (slashPos==0) then   ! it is all filename
            dir = ""
            fullName = path
        elseif (slashPos==pathLen) then   ! it is all directory
            dir = path
            fullName = ""
            return
        else
            dir = path(1:slashPos)
            fullName = path(slashPos+1:pathLen)
        end if

    end subroutine getDirFullName

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! returns the name and file extension of the input full file name.
    subroutine getNameExt(fullName,name,ext)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getNameExt
#endif
        use Constants_mod, only: IK
        implicit none
        character(*), intent(in)                :: fullName
        character(:), allocatable, intent(out)  :: name,ext
        integer(IK)                             :: dotPos, lenFilename
        lenFilename = len(fullName)
        if (lenFilename==0) then
            name = ""; ext = ""
            return
        else
            dotPos = index(fullName,".",back=.true.)
            if ( dotPos==0 .or. dotPos==lenFilename ) then     ! there is no extension
                name = fullName
                ext = ""
            elseif ( dotPos==1 ) then   ! it is all extension, no file name.
                name = ""
                ext  = fullName(1:)
            else
                name = fullName(1:dotPos-1)
                ext  = fullName(dotPos:)
            end if
        end if
    end subroutine getNameExt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! this routine does not check for OS type
    function mkdir(dirPath,isWindows,wait) result(Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: mkdir
#endif
        use System_mod, only: SysCmd_type
        use String_mod, only: num2str
        use Err_mod, only: Err_type
        implicit none
        character(*), parameter         :: PROCEDURE_NAME = "@mkdir()"
        character(*), intent(in)        :: dirPath
        logical, intent(in), optional   :: isWindows, wait
        type(Err_type)                  :: Err
        type(SysCmd_type)               :: SysCmd
        Err%occurred= .false.
        if (present(isWindows)) then
            if (isWindows) then
                SysCmd = SysCmd_type('mkdir "'//dirPath//'" >nul 2>&1',wait) ! path has to be enclosed with "" to allow nested mkdir
            else
                SysCmd = SysCmd_type("mkdir -p "//dirPath//" > /dev/null 2>&1",wait) ! -p enables nested mkdir
            end if
        else
            SysCmd = SysCmd_type("mkdir "//dirPath,wait)
        end if
        if (SysCmd%Err%occurred) then
            Err%occurred = .true.
            Err%stat = SysCmd%Err%stat
            Err%msg = PROCEDURE_NAME // SysCmd%Err%msg // "\nexecute_command_line() exitstat: " // num2str(SysCmd%exitstat)
        end if
    end function mkdir

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Path_mod