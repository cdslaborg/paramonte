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

!> This module contains classes and procedures for manipulating system file/folder paths.
!>  \author Amir Shahmoradi

module Path_mod

    use Constants_mod, only: IK ! LCOV_EXCL_LINE
    use Err_mod, only: Err_type
    implicit none

    character(*), parameter :: MODULE_NAME = "@Path_mod"

    integer(IK), parameter  :: MAX_FILE_PATH_LEN = 2047

    !> Windows reserved characters (not allowed in filenames):
    character(*), parameter :: WINDOWS_RESERVED_CHAR = "<>:" // '"' // "|?*" ! /\

#if defined INTEL_COMPILER_ENABLED

    character(*), parameter :: SHELL_ESCAPE_CHAR =  &
                                                    " " // & ! space character
                                                    "!" // & ! history expansion.
                                                    '"' // & ! shell syntax.
                                                    "#" // & ! comment start when preceded by whitespace; zsh wildcards.
                                                    "$" // & ! shell syntax.
                                                    "&" // & ! shell syntax.
                                                    "'" // & ! shell syntax.
                                                    "(" // & ! even in the middle of a word: ksh extended globs (also available in bash and zsh); zsh wildcards.
                                                    ")" // & ! even in the middle of a word: ksh extended globs (also available in bash and zsh); zsh wildcards.
                                                    "*" // & ! sh wildcard.
                                                    "," // & ! only inside brace expansion.
                                                    ";" // & ! shell syntax.
                                                    "<" // & ! shell syntax.
                                                    "=" // & ! in zsh, when it is at the beginning of a file name (filename expansion with PATH lookup).
                                                    ">" // & ! shell syntax.
                                                    "?" // & ! sh wildcard.
                                                    "[" // & ! sh wildcard.
                                                    "\" // & ! shell syntax.
                                                    "]" // & ! you may get away with leaving it unquoted.
                                                    "^" // & ! history expansion; zsh wildcard.
                                                    "`" // & ! shell syntax.
                                                    "{" // & ! brace expansion.
                                                    "|" // & ! shell syntax.
                                                    "}" // & ! needs to be escaped in zsh, other shells are more lenient when there is no matching opening brace.
                                                    "~"      ! home directory expansion when at the beginning of a filename; zsh wildcard; safe when it is the last character.

#else

    ! stupid gfortran (possibly version 8.3) gives error on the above syntax
    character(*), parameter :: SHELL_ESCAPE_CHAR = " !"//'"#$&'//"'()*,;<=>?[\]^`{|}~"

#endif

    ! The `Path_type` class.
    type :: Path_type
        character(:), allocatable       :: original     !< The original path.
        character(:), allocatable       :: modified     !< The modified path based on the OS/platform type.
        character(:), allocatable       :: dir          !< The directory segment of the path.
        character(:), allocatable       :: name         !< The name of the file, if any exists in the path.
        character(:), allocatable       :: base         !< The base of the file name, if any exists in the path.
        character(:), allocatable       :: ext          !< The file extension, if any exists in the path (including the dot separator).
        character(1)                    :: shellSlash   !< The type of the separator (forward/backward slash) with which the original path is *modified*.
        type(Err_type)                  :: Err          !< An object of class [Err_type](@ref err_mod::err_type) containing error handling tools.
    contains
        procedure, pass                 :: query
        procedure, nopass               :: modify
        procedure, nopass               :: getDirNameExt, getDirFullName, getNameExt
        procedure, nopass               :: winify, linify
        procedure, nopass               :: mkdir
        procedure, nopass               :: isdir
    end type Path_type

    interface Path_type
        module procedure :: constructPath
    end interface Path_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This is the constructor of the class [Path_type](@ref path_type).\n
    !> Return an object of class [Path_type](@ref path_type) given the input specifications.
    !>
    !> \param[in]   inputPath   :   The input path.
    !> \param[in]   OS          :   An object of class [OS_type](@ref system_mod::os_type) containing information about the operating system (**optional**).
    !>
    !> \return
    !> `Path` : An object of class [Path_type](@ref path_type) containing the path properties and methods.
    function constructPath(inputPath,OS) result(Path)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
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

    !> \brief
    !> This procedure is a method of the class [Path_type](@ref path_type).\n
    !> Construct an object of class [Path_type](@ref path_type) as output.
    !>
    !> \param[inout]    Path        :   An object of class [Path_type](@ref path_type) containing the path properties and methods.
    !> \param[in]       inputPath   :   The input path (**optional**). If provided, it will overwrite `Path%original`.
    !> \param[in]       OS          :   An object of class [OS_type](@ref system_mod::os_type) containing information about the operating system (**optional**).
    !>
    !> \warning
    !> On output, do not forget to check the value `Path%%Err%%occurred` before using the output `Path`.
    subroutine query(Path,inputPath,OS)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: query
#endif
        use Err_mod, only: Err_type
        use Constants_mod, only: IK
        use System_mod, only: OS_type
        use String_mod, only: replaceStr
        implicit none
        class(Path_type), intent(inout)     :: Path
        character(*), intent(in), optional  :: inputPath
        type(OS_type), intent(in), optional :: OS
        logical                             :: isUnixShell

        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME//"@query()"

        Path%Err%occurred = .false.
        Path%Err%msg = ""

        if (present(inputPath)) then
            Path%original = trim(adjustl(inputPath))
        elseif (.not.allocated(Path%original)) then
            Path%Err%occurred = .true.
            Path%Err%msg = PROCEDURE_NAME//": Error occurred. Neither inputPath argument is given as input, nor Path%original is allocated to construct the Path object."
            return
        else
            if ( len(trim(adjustl(Path%original)))==0 ) then
                Path%Err%occurred = .true.
                Path%Err%msg = PROCEDURE_NAME//": Error occurred. Neither inputPath argument is given as input, nor Path%original has a non-blank length > 0 to construct the Path object."
                return
            end if
        end if

        if (present(OS)) then
            Path%shellSlash = OS%Shell%slash
            isUnixShell = OS%Shell%isUnix
        else
            block
                type(OS_type) :: OS
                call OS%query()
                if (OS%Err%occurred) then
                ! LCOV_EXCL_START
                    Path%Err%stat = OS%Err%stat
                    Path%Err%occurred = OS%Err%occurred
                    Path%Err%msg = PROCEDURE_NAME // ": Error occurred while querying OS type.\n" // Path%Err%msg
                end if
                ! LCOV_EXCL_STOP
                Path%shellSlash = OS%Shell%slash
                isUnixShell = OS%Shell%isUnix
            end block
            if (Path%Err%occurred) return
        end if

        if (isUnixShell) then
            ! if the path contains both / and \, then assume that it is already in linux style
            if (index(Path%original,"/")==0) then ! path is given in Windows style
                Path%modified = linify(Path%original)
            else
                Path%modified = Path%original
            end if
#if defined OS_IS_WINDOWS
        else
            Path%modified = winify(Path%original)
#endif
        end if

        call Path%getDirNameExt( Path%modified, Path%shellSlash, Path%dir, Path%name, Path%ext )
        Path%base = Path%dir // Path%name

    end subroutine query

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a static method of the class [Path_type](@ref path_type).\n
    !> Convert the the input path to the modified path according to the rules of the Windows operating system.
    !>
    !> \param[in]       inputPath   :   The input path. If provided, it will overwrite `Path%original`.
    !>
    !> \return
    !> `outputPath` : The output modified path which conforms to the rules of the Windows OS.
    !>
    !> \warning
    !> This code assumes that the input path is a Linux path. Windows paths like `.\(paramonte)\paramonte.nml` will be horribly
    !> treated by this routine as `\(` also represents a Linux escape character. The result will be `.(paramonte)\paramonte.nml`.
    !>
    !> \warning
    !> This routine strictly assumes that there is no dangling `\` in the input Linux path, and if there is,
    !> then either it is used to escape the special shell characters, or otherwise, the path is a Windows path.
    pure function winify(inputPath) result(outputPath) !,Err)!,ignoreWindowsReservedChars)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: winify
#endif
       !use Err_mod, only: Err_type
        use Constants_mod, only: IK
        use String_mod, only: replaceStr
        implicit none
        character(len=*), intent(in)            :: inputPath
        character(:), allocatable               :: outputPath
       !type(Err_type), intent(out)             :: Err
       !logical, intent(in), optional           :: ignoreWindowsReservedChars
       !logical                                 :: reservedCharInspectionNeeded
       !character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME//"@winify()"
        character(:), allocatable               :: outputPathDummy
        integer(IK)                             :: i, j, outputPathLen

        !Err%occurred = .false.
        !Err%msg = ""

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
            outputPathDummy = replaceStr(outputPath,"\\","\")
            if (outputPathDummy==outputPath) exit loopRemoveMultipleSlash
            outputPath = outputPathDummy
        end do loopRemoveMultipleSlash
        outputPathLen = len(outputPath)

        ! Now check for the presence of any Linux Shell Escape Character in the input path without a preceding \.
        ! If there is any, this would imply that the input path is a Windows path,
        ! otherwise a escape character without preceding \ would be invalid in Linux.

        if (outputPathLen==1_IK) then
            if (outputPath=="/") outputPath = "\"
            return
        else
            do i = 1, len(SHELL_ESCAPE_CHAR)
                if (SHELL_ESCAPE_CHAR(i:i)/="\") then
                    do j = 2, outputPathLen
                        if (outputPath(j:j)==SHELL_ESCAPE_CHAR(i:i)) then
                            if (outputPath(j-1:j-1)/="\") return ! no escaping has occurred. Therefore, it is a windows path, there is no need for further winifying.
                        end if
                    end do
                end if
            end do
        end if

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

    end function winify

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This `pure` procedure is a static method of the class [Path_type](@ref path_type).\n
    !> Convert the the input path to the modified path according to the rules of the Unix operating systems.
    !>
    !> \param[in]   inputPath   :   The input path. If provided, it will overwrite `Path%original`.
    !>
    !> \return
    !> `outputPath` : The output modified path which conforms to the rules of the Unix OS.
    pure function linify(inputPath) result(outputPath)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: linify
#endif
        use Constants_mod, only: IK
        use String_mod, only: replaceStr
        implicit none
        character(*), intent(in)    :: inputPath
        character(:), allocatable   :: outputPath
        character(:), allocatable   :: outputPathDummy
        integer(IK)                 :: i

        ! check if the path is sandwiched between quotation marks. If so, remove them:
        outputPath = trim(adjustl(inputPath))
        i = len(outputPath)
        if (i==0) return
        if ( i>1 ) then
            if ( (outputPath(1:1)=='"' .and. outputPath(i:i)=='"') .or. (outputPath(1:1)=="'" .and. outputPath(i:i)=="'") )then
                outputPathDummy = outputPath(2:i-1)
            else
                outputPathDummy = outputPath
            end if
        end if

        ! First change all backslashes to forward slash:
        outputPath = replaceStr(outputPathDummy,"\","/")

        ! Now correct for any Linux Shell Escape Character in the input path:
        do i = 1, len(SHELL_ESCAPE_CHAR)
            if (SHELL_ESCAPE_CHAR(i:i)/="\") then
                outputPathDummy = replaceStr(outputPath,SHELL_ESCAPE_CHAR(i:i),"\"//SHELL_ESCAPE_CHAR(i:i))
                outputPath = outputPathDummy
            end if
        end do

        !! Now correct for any white spaces in outputPath:
        !outputPath = replaceStr(outputPath," ","\ ")

    end function linify

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a static method of the class [Path_type](@ref path_type).\n
    !> Modify the input path to conform to the rules of the current inferred operating system.
    !>
    !> \param[in]       inputPath   :   The input path. If provided, it will overwrite `Path%original`.
    !> \param[out]      outputPath  :   The output modified path which conforms to the rules of the current OS.
    !> \param[out]      Err         :   An object of class [Err_type](@ref err_mod::err_type) containing error handling tools.
    subroutine modify(inputPath,outputPath,Err)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: modify
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

        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME//"@modify()"

        outputPath = trim(adjustl(inputPath))

        Err%occurred = .false.
        Err%msg = ""

        call OS%query()

        if (OS%Err%occurred) then
        ! LCOV_EXCL_START
            Err = OS%Err
            Err%msg = PROCEDURE_NAME // ": Error occurred while modifying inputPath='" // outputPath // "'.\n" // Err%msg
            return
        end if
        ! LCOV_EXCL_STOP

        if (OS%Shell%isUnix) then
            outputPath = linify(inputPath)
#if defined OS_IS_WINDOWS
        else
            outputPath = winify(inputPath)
#endif
        end if

    end subroutine modify

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a static method of the class [Path_type](@ref path_type).\n
    !> Split the input path to directory, base file name, and the file extension, based on the input OS slash.
    !>
    !> \param[in]       path    :   The input path.
    !> \param[in]       slash   :   The separator used by the operating system to delimit segments of a path.
    !> \param[out]      dir     :   The directory segment of the path.
    !> \param[out]      name    :   The base file name segment of the path.
    !> \param[out]      ext     :   The file extension segment of the path.
    subroutine getDirNameExt(path,slash,dir,name,ext)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirNameExt
#endif
        implicit none
        character(*)             , intent(in)   :: path
        character(1)             , intent(in)   :: slash
        character(:), allocatable, intent(out)  :: dir, name, ext
        character(:), allocatable               :: fullName
        call getDirFullName(path,slash,dir,fullName)
        call getNameExt(fullName,name,ext)
    end subroutine getDirNameExt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a static method of the class [Path_type](@ref path_type).\n
    !> Return the directory and full filename (including the file extension) of the input path.
    !>
    !> \param[in]       path        :   The input path.
    !> \param[in]       slash       :   The separator used by the operating system to delimit segments of a path.
    !> \param[out]      dir         :   The directory segment of the path.
    !> \param[out]      fullName    :   The full file name and extension segment of the path.
    subroutine getDirFullName(path,slash,dir,fullName)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: getDirFullName
#endif
        use Constants_mod, only: IK
        implicit none
        character(*)             , intent(in)   :: path
        character(1)             , intent(in)   :: slash
        character(:), allocatable, intent(out)  :: dir, fullName

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

    !> \brief
    !> This procedure is a static method of the class [Path_type](@ref path_type).\n
    !> Return the name and file extension of the input full file name.
    !>
    !> \param[in]       fullName    :   The full file name and extension of the path.
    !> \param[out]      name        :   The name segment of the file.
    !> \param[out]      ext         :   The extension segment of the file (including the dot separator).
    subroutine getNameExt(fullName,name,ext)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
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
            else
                name = fullName(1:dotPos-1)
                ext  = fullName(dotPos:)
            end if
        end if
    end subroutine getNameExt

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a static method of the class [Path_type](@ref path_type).\n
    !> Make the requested (nested) directory (recursively, if needed).
    !>
    !> \param[in]   dirPath     :   The full directory path.
    !> \param[in]   isUnixShell :   The logical flag indicating whether the OS is Windows (**optional**).
    !>                              If not present, the runtime shell type will be inferred by the procedure.
    !> \param[in]   wait        :   The logical flag indicating whether the procedure should wait for the system
    !>                              operation to complete and return (**optional**, default = `.true.`).
    !>
    !> \return
    !> `Err` : An object of class [Err_type](@ref err_mod::err_type), indicating whether an error has occurred while creating the directory.
    !>
    !> \author
    !> Last updated by Amir Shahmoradi, Tuesday 3:09 AM, Dec 8, 2020, Dallas, TX
    function mkdir(dirPath,isUnixShell,wait) result(Err)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: mkdir
#endif
        use Constants_mod, only: IK, NLC
        use System_mod, only: SysCmd_type, OS_type
        use String_mod, only: num2str
        use Err_mod, only: Err_type
        implicit none
        character(*), parameter         :: PROCEDURE_NAME = MODULE_NAME//"@mkdir()"
        character(*), intent(in)        :: dirPath
        logical, intent(in), optional   :: isUnixShell, wait
        type(Err_type)                  :: Err
        type(SysCmd_type)               :: SysCmd
        type(OS_type)                   :: OS
        logical                         :: isUnixShellDefault
        character(:), allocatable       :: command
        integer(IK)                     :: itry



        Err%occurred = .false.

        if (present(isUnixShell)) then
            isUnixShellDefault = isUnixShell
        else
            OS%Err%occurred = .false.
            call OS%query()
            if (OS%Err%occurred) then
                ! LCOV_EXCL_START
                command = "mkdir "//dirPath
                ! LCOV_EXCL_STOP
            else
                isUnixShellDefault = OS%Shell%isUnix
            end if
        end if

        if (.not. allocated(command)) then
            if (isUnixShellDefault) then
                command = "mkdir -p "//dirPath//" > /dev/null 2>&1" ! -p enables nested mkdir
#if defined OS_IS_WINDOWS
            else
                command = 'mkdir "'//dirPath//'" >nul 2>&1' ! WARNING: path has to be enclosed with "" to allow nested mkdir
#endif
            end if
        end if

        ! Try to create the folder for 10 times, and fail if all attempts fail.

        loopTry: do itry = 1, 10
            SysCmd = SysCmd_type(command, wait)
            if (SysCmd%Err%occurred .and. .not. isdir(dirPath)) cycle loopTry
            deallocate(command)
            return
        end do loopTry

        ! LCOV_EXCL_START
        Err%occurred = .true.
        Err%stat = SysCmd%Err%stat
        Err%msg = PROCEDURE_NAME // SysCmd%Err%msg //NLC//"execute_command_line() exitstat: " // num2str(SysCmd%exitstat)
        ! LCOV_EXCL_STOP

    end function mkdir

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a static method of the class [Path_type](@ref path_type).\n
    !> Return `.true.` if the input path is a directory, otherwise, return `.false.`.
    !>
    !> \param[in]   path    :   The full directory path.
    !>
    !> \return
    !> `pathIsDir` : A logical output variable indicating whether the input path is a directory.
    !>
    !> \author
    !> Amir Shahmoradi, Tuesday 3:09 AM, Dec 8, 2020, Dallas, TX
    function isdir(path) result(pathIsDir)
#if INTEL_COMPILER_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN)
        !DEC$ ATTRIBUTES DLLEXPORT :: isdir
#endif
        implicit none
        character(*), intent(in)        :: path
        logical                         :: pathIsDir
#if defined INTEL_COMPILER_ENABLED
        inquire(directory = path, exist = pathIsDir)
#elif defined GNU_COMPILER_ENABLED
        inquire(file = path, exist = pathIsDir)
#else
#error "This procedure does not currently support compilers other than Intel ifort and GNU gfortran."
#endif
    end function isdir

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module Path_mod ! LCOV_EXCL_LINE