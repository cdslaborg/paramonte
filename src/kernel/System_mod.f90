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

!>  \brief This module contains classes and procedures relevant to the system operations.
!>  @author Amir Shahmoradi

module System_mod

    use JaggedArray_mod, only: CharVec_type
    use Constants_mod, only: IK, NLC
    use Err_mod, only: Err_type
    implicit none

    character(*), parameter             :: MODULE_NAME = "@System_mod"
    integer(IK) , parameter             :: MAX_OS_NAME_LEN = 63_IK

    !> The `RandomFileName_type` class.
    type :: RandomFileName_type
        character(:), allocatable       :: path                     !< The full path to the randomly-generated unique file name.
        character(:), allocatable       :: dir                      !< The directory within which is the unique new file is supposed to be generated.
        character(:), allocatable       :: key                      !< The optionally user-specified file prefix for the unique file name.
        character(:), allocatable       :: ext                      !< The optionally user-specified file extension.
        type(Err_type)                  :: Err                      !< An object of class [Err_type](@ref err_mod::err_type).
    end type RandomFileName_type

    !> The `RandomFileName_type` constructor.
    interface RandomFileName_type
        module procedure                :: getRandomFileName
    end interface RandomFileName_type

    !> The `SystemInfo_type` class.
    type :: SystemInfo_type
        integer(IK)                     :: nRecord                  !< The number of elements of the vector `List`.
        type(CharVec_type), allocatable :: List(:)                  !< An array of length `nRecord` of strings, each element of which represents
                                                                    !! one line in the output system information.
        type(Err_type)                  :: Err                      !< An object of class [Err_type](@ref err_mod::err_type) indicating whether
                                                                    !! any error has occurred during information collection.
    contains
        procedure, nopass               :: get => getSystemInfo
    end type SystemInfo_type

    !> The `SystemInfo_type` constructor.
    interface SystemInfo_type
        module procedure                :: constructSystemInfo
    end interface SystemInfo_type

    !> The Shell name type.
    type, private :: ShellName_type
        character(:), allocatable       :: current                  !< The name of the current runtime shell.
        character(:), allocatable       :: default                  !< The name of the default runtime shell.
    end type ShellName_type

    !> The `Shell_type` class.
    type :: Shell_type
        logical                         :: isSh         = .false.   !< The logical value indicating whether the shell is Unix sh.
        logical                         :: isCMD        = .false.   !< The logical value indicating whether the shell is Windows CMD.
        logical                         :: isZsh        = .false.   !< The logical value indicating whether the shell is Unix zsh.
        logical                         :: isCsh        = .false.   !< The logical value indicating whether the shell is Unix csh.
        logical                         :: isBash       = .false.   !< The logical value indicating whether the shell is Unix Bash.
        logical                         :: isPowerShell = .false.   !< The logical value indicating whether the shell is Windows PowerShell.
        logical                         :: isUnix       = .false.   !< The logical value indicating whether the shell is Unix-like.
        character(:), allocatable       :: name                     !< The name of or path to the current shell.
        type(Err_type)                  :: Err                      !< An object of class [Err_type](@ref err_mod::err_type) indicating
                                                                    !! whether error has occurred during the query.
    contains
        procedure, pass                 :: query => queryRuntimeShell
    end type Shell_type

    !> The `OS_type` class.
    type :: OS_type
        character(:), allocatable       :: name                     !< The name of the operating system.
        character(:), allocatable       :: slash                    !< The file/folder name separator used by the OS.
        logical                         :: isWindows = .false.      !< Logical variable indicating whether the OS is Windows.
        logical                         :: isDarwin = .false.       !< Logical variable indicating whether the OS is Darwin (macOS).
        logical                         :: isLinux = .false.        !< Logical variable indicating whether the OS is Linux.
        type(Shell_type)                :: Shell                    !< An object of class [Shell_type](@ref shell_type) containing
                                                                    !! information about the runtime shell name and type.
        type(Err_type)                  :: Err                      !< An object of class [Err_type](@ref err_mod::err_type) indicating whether
                                                                    !! error has occurred during the object initialization.
    contains
        procedure, pass                 :: query => queryOS
    end type OS_type

    !> The `EnvVar_type` class.
    type :: EnvVar_type
        character(:), allocatable       :: name
        character(:), allocatable       :: value
        integer                         :: length
        type(Err_type)                  :: Err
    contains
        procedure, nopass               :: get => getEnvVar
    end type EnvVar_type

    !> The `CmdArg_type` class.
    type :: CmdArg_type
        character(:), allocatable       :: cmd      !< A string containing the full command line obtained via `get_command()` Fortran intrinsic subroutine.
        type(CharVec_type), allocatable :: Arg(:)   !< A list of `(0:CmdArg_type%count)` elements, each of which represents one command line argument,
                                                    !! including the main command as the zeroth element.
        integer                         :: count    !< The number of command line arguments, excluding the main (zeroth) command.
        type(Err_type)                  :: Err      !< An object of class [Err_type](@ref err_mod::err_type) indicating
                                                    !! whether error has occurred during the object initialization.
    contains
        procedure, pass                 :: query => queryCmdArg
    end type CmdArg_type

    !> The `SysCmd_type` class.
    type :: SysCmd_type
        character(:), allocatable       :: cmd      !< The command to be executed by the program in the terminal.
        logical                         :: wait     !< Indicated if the program should wait for the terminal to return the control to it.
        integer                         :: exitstat !< The exit status from the terminal.
        type(Err_type)                  :: Err      !< An object of class [Err_type](@ref err_mod::err_type) indicating
                                                    !! whether error has occurred during the object initialization.
    contains
        procedure, pass                 :: run => runSysCmd
    end type SysCmd_type

    !> The `SysCmd_type` constructor.
    interface SysCmd_type
        module procedure :: constructSysCmd
    end interface SysCmd_type

    ! cache the OS query result to speed up code

    logical                             :: mv_osCacheEnabled = .false.
    logical                             :: mv_shCacheEnabled = .false.
    type(OS_type)                       :: mv_OS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The constructor of the class [SystemInfo_type](@ref systeminfo_type).
    !> Return a comprehensive report of the system information.
    !>
    !> \param[in]   OS      :   An object of class [OS_type](@ref os_type) (optional).
    !>
    !> \return
    !> `SystemInfo` : An object of class [SystemInfo_type](@ref systeminfo_type) containing the system information.
    function constructSystemInfo(OS) result(SystemInfo)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructSystemInfo
#endif
        implicit none
        type(SystemInfo_type)               :: SystemInfo
        type(OS_type), intent(in), optional :: OS
        call getSystemInfo( SystemInfo%List, SystemInfo%Err, OS, SystemInfo%nRecord )
    end function constructSystemInfo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Query all attributes of the [OS_type](@ref os_type) class: `name`, `slash`, `isWindows`, `Err`.
    !>
    !> \param[out]  OS                  :   An object of class [OS_type](@ref os_type).
    !> \param[in]   shellQueryEnabled   :   A logical variable indicating if the type and name of the current
    !>                                      runtime shell should be queried or not (optional, default = `.false.`).
    subroutine queryOS(OS, shellQueryEnabled)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: queryOS
#endif

        use String_mod, only: num2str, getLowerCase
        use Constants_mod, only: IK, RK
        use Err_mod, only: Err_type
        implicit none
        class(OS_type)  , intent(out)           :: OS
        logical         , intent(in), optional  :: shellQueryEnabled
        character(*)    , parameter             :: PROCEDURE_NAME = MODULE_NAME // "@queryOS()"
        logical                                 :: shellQueryEnabledDefault
#if !defined OS_IS_WINDOWS && !defined OS_IS_DARWIN && !defined OS_IS_LINUX
        character(:)    , allocatable           :: osname
#endif

        shellQueryEnabledDefault = .false.
        if (present(shellQueryEnabled)) shellQueryEnabledDefault = shellQueryEnabled
        OS%Err%occurred = .false.
        OS%Err%msg = ""

        if (mv_osCacheEnabled) then

            OS%name         = mv_OS%name
            OS%slash        = mv_OS%slash
            OS%isWindows    = mv_OS%isWindows
            OS%isDarwin     = mv_OS%isDarwin
            OS%isLinux      = mv_OS%isLinux

            if (mv_shCacheEnabled) then
                OS%Shell    = mv_OS%Shell
            else
                mv_shCacheEnabled = .true.
                call OS%Shell%query(OS%isWindows)
                if (OS%Shell%Err%occurred) then
                    OS%Err = OS%Shell%Err
                    return
                end if
                mv_OS%Shell = OS%Shell
            end if

            return

        end if

#if defined OS_IS_WINDOWS

        OS%isWindows = .true.
        OS%name = "Windows"
        OS%slash = "\"

#elif defined OS_IS_DARWIN

        OS%isDarwin = .true.
        OS%name = "Darwin"
        OS%slash = "/"

#elif defined OS_IS_LINUX

        OS%isLinux = .true.
        OS%name = "Linux"
        OS%slash = "/"

#else

        if (allocated(OS%name)) deallocate(OS%name); allocate( character(MAX_OS_NAME_LEN) :: OS%name )
        call getEnvVar( name="OS", value=OS%name, Err=OS%Err )
        if (OS%Err%occurred) then
            OS%Err%msg = PROCEDURE_NAME // ": Error occurred while querying OS type." // NLC // OS%Err%msg
            OS%name = ""
            return
        end if

        OS%name = trim(adjustl(OS%name))

        blockOS: if (len(OS%name)>=7_IK) then

            if (getLowerCase(OS%name(1:7))=="windows") then

                OS%isWindows = .true.
                OS%isDarwin = .false.
                OS%isLinux = .false.
                OS%slash = "\"

            end if

        else blockOS ! it is either Linux- or Darwin- based OS

            if (allocated(OS%name)) deallocate( OS%name )
            allocate( character(MAX_OS_NAME_LEN) :: OS%name )
            OS%isWindows = .false.
            OS%slash = "/"

            if (allocated(OS%name)) deallocate(OS%name); allocate( character(MAX_OS_NAME_LEN) :: OS%name )
            call getEnvVar( name="OSTYPE", value=OS%name, Err=OS%Err )
            if (OS%Err%occurred) then
                OS%Err%msg = PROCEDURE_NAME // ": Error occurred while querying OS type." // NLC // OS%Err%msg
                OS%name = ""
                return
            end if

            OS%name = trim(adjustl(OS%name))
            osname = getLowerCase(OS%name)

            blockNonWindowsOS: if (index(osname,"darwin")/=0) then

                OS%isDarwin = .true.
                OS%isLinux = .false.
                return

            elseif (index(osname,"linux")/=0) then blockNonWindowsOS

                OS%isDarwin = .false.
                OS%isLinux = .true.
                return

            else blockNonWindowsOS

                if (allocated(OS%name)) deallocate(OS%name); allocate( character(MAX_OS_NAME_LEN) :: OS%name )

                blockUnknownOS: block

                    integer                     :: fileUnit
                    type(RandomFileName_type)   :: RFN
                    RFN = RandomFileName_type(key="queryOS")
                    if (RFN%Err%occurred) then
                        OS%Err = RFN%Err
                        OS%Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring OS type." // NLC // OS%Err%msg
                        OS%name = ""
                        return
                    end if

                    call executeCmd( command="uname > "//RFN%path, Err=OS%Err )
                    if (OS%Err%occurred) then
                        OS%Err%msg = PROCEDURE_NAME // ": Error occurred while executing command 'uname > "// RFN%path // "'." // NLC // OS%Err%msg
                        OS%name = ""
                        return
                    end if

                    open(newunit=fileUnit,file=RFN%path,status="old",iostat=OS%Err%stat)
                    if (OS%Err%stat>0) then
                        OS%Err%occurred = .true.
                        OS%Err%msg =    PROCEDURE_NAME // ": Unknown error occurred while opening file = '" // RFN%path // "'."
                        OS%name = ""
                        return
                    end if

                    read(fileUnit,*,iostat=OS%Err%stat) OS%name

                    if ( is_iostat_eor(OS%Err%stat) ) then
                        OS%Err%occurred = .true.
                        OS%Err%msg =    PROCEDURE_NAME // ": End-Of-Record error condition occurred while attempting to read &
                                        &the Operating System's name from file = '" // RFN%path // "'."
                        OS%name = ""
                        return
                    elseif ( is_iostat_end(OS%Err%stat) ) then
                        OS%Err%occurred = .true.
                        OS%Err%msg =    PROCEDURE_NAME // ": End-Of-File error condition occurred while attempting to read &
                                        &the Operating System's name from file = '" // RFN%path // "'."
                        OS%name = ""
                        return
                    elseif ( OS%Err%stat>0 ) then
                        OS%Err%occurred = .true.
                        OS%Err%msg =    PROCEDURE_NAME // ": Unknown error condition occurred while attempting to read &
                                        &the Operating System's name from file = '" // RFN%path // "'."
                        OS%name = ""
                        return
                    end if

                    close(fileUnit, status = "delete")

                    OS%name = trim(adjustl(OS%name))
                    osname = getLowerCase(OS%name)
                    if (index(osname,"darwin")/=0) then
                        OS%isDarwin = .true.
                        OS%isLinux = .false.
                    elseif (index(osname,"linux")/=0) then
                        OS%isLinux = .true.
                        OS%isDarwin = .false.
                    else
                        OS%isLinux = .false.
                        OS%isDarwin = .false.
                    end if

                end block blockUnknownOS

            end if blockNonWindowsOS

        end if blockOS

#endif


        mv_osCacheEnabled = .true.
        mv_OS%name      = OS%name
        mv_OS%slash     = OS%slash
        mv_OS%isWindows = OS%isWindows
        mv_OS%isDarwin  = OS%isDarwin
        mv_OS%isLinux   = OS%isLinux

        if (shellQueryEnabledDefault) then

            if (mv_shCacheEnabled) then
                OS%Shell    = mv_OS%Shell
            else
                mv_shCacheEnabled = .true.
                call OS%Shell%query(OS%isWindows)
                if (OS%Shell%Err%occurred) then
                    OS%Err = OS%Shell%Err
                    return
                end if
                mv_OS%Shell = OS%Shell
            end if

        end if

    end subroutine queryOS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine queryRuntimeShell(Shell, isWindowsOS)

        use FileContents_mod, only: FileContents_type

        implicit none

        class(Shell_type), intent(inout)    :: Shell
        logical, intent(in)                 :: isWindowsOS

        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@queryRuntimeShell()"

        type(RandomFileName_type)           :: RFN
        type(FileContents_type)             :: FileContents
        character(:), allocatable           :: command
        logical                             :: fileExists

        ! create a random output file name

        RFN = RandomFileName_type(key="queryShell")
        if (RFN%Err%occurred) then
            Shell%Err = RFN%Err
            Shell%Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring OS type." // NLC // Shell%Err%msg
            Shell%name = ""
            return
        end if

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! define the shell command. First try the bash command,
        ! as it does not lead to oddities on Windows terminal.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !command = "echo $0 >" // RFN%path // " 2>&1 && echo $SHELL >" // RFN%path // " 2>&1"
        command = "echo $0 >" // RFN%path // " 2>&1"
        call executeCmd( command = command, Err = Shell%Err )
        inquire(file = RFN%path, exist = fileExists)
        if (Shell%Err%occurred .or. .not. fileExists) then
            Shell%Err%msg = PROCEDURE_NAME // ": Error occurred while executing the Unix command "// command // NLC // Shell%Err%msg
            Shell%name = ""
            return
        end if

        ! read the command output

        FileContents = FileContents_type(RFN%path, delEnabled = .true.)
        if (FileContents%Err%occurred) then
            Shell%Err%occurred = .true.
            Shell%Err%msg = PROCEDURE_NAME // FileContents%Err%msg
            Shell%name = ""
            return
        end if

        if (FileContents%numRecord>0_IK) then
            Shell%name      = trim(adjustl(FileContents%Line(1)%record))
            Shell%isZsh     = index(Shell%name,"zsh") > 0
            Shell%isCsh     = index(Shell%name,"csh") > 0
            Shell%isBash    = index(Shell%name,"bash") > 0
            Shell%isSh      = .false.; if (.not. (Shell%isBash .or. Shell%isZsh .or. Shell%isCsh)) Shell%isSh = index(Shell%name,"sh") > 0
            Shell%isUnix    = (.not. isWindowsOS) .or. Shell%isBash .or. Shell%isZsh .or. Shell%isCsh .or. Shell%isSh
        end if

        if (Shell%isUnix) return

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! define the shell command, this time for Windows Batch.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (isWindowsOS) then

            command = "(dir 2>&1 *`|echo CMD >temp.txt);&<# rem #>echo PowerShell >temp.txt 2>&1"

            call executeCmd( command = command, Err = Shell%Err )
            if (Shell%Err%occurred .or. .not. fileExists) then
                Shell%Err%msg = PROCEDURE_NAME // ": Error occurred while executing the Windows command "// command // NLC // Shell%Err%msg
                Shell%name = ""
                return
            end if

            ! read the command output

            FileContents = FileContents_type(RFN%path, delEnabled = .true.)
            if (FileContents%Err%occurred) then
                Shell%Err%occurred = .true.
                Shell%Err%msg = PROCEDURE_NAME // FileContents%Err%msg
                Shell%name = ""
                return
            end if

            if (FileContents%numRecord>0_IK) then
                Shell%isCMD = index(Shell%name,"CMD") > 0
                Shell%isPowerShell = index(Shell%name,"PowerShell") > 0
            end if

        end if

    end subroutine queryRuntimeShell

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Generate a unique file path in the requested directory for temporary usage.
    !>
    !> \param[in]   dir : The requested directory within which the unique new file is supposed to be generated (optional).
    !> \param[in]   key : The requested input file name prefix (optional, default = "RandomFileName").
    !> \param[in]   ext : The requested input file extension (optional, default = ".rfn", standing for random file name).
    !>
    !> \return
    !> `RFN` : An object of class [RandomFileName_type](@ref randomfilename_type) containing the attributes of the random file name.
    function getRandomFileName(dir,key,ext) result(RFN)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRandomFileName
#endif
        use Constants_mod, only: IK, RK
        use DateTime_mod, only: DateTime_type
        use String_mod, only: num2str
        implicit none
        character(*), intent(in), optional  :: dir, key, ext
        type(RandomFileName_type)           :: RFN

        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@getRandomFileName()"

        integer(IK)                         :: counter
        logical                             :: fileExists
        type(DateTime_type)                 :: DT

        if (present(dir)) then
            RFN%dir = dir
        else
            RFN%dir = ""
        end if
        if (present(key)) then
            RFN%key = key
        else
            RFN%key = "RandomFileName"
        end if
        if (present(ext)) then
            RFN%ext = ext
        else
            RFN%ext = ".rfn"
        end if

        counter = 0
        do

            counter = counter + 1
            call DT%query()
#if defined CAF_ENABLED
            RFN%path = RFN%dir // RFN%key // '_' // DT%date // '_' // DT%time // '_process_' // num2str(this_image())   // '_' // num2str(counter) // RFN%ext
#elif defined MPI_ENABLED
            block
            use mpi
            integer :: imageID, ierrMPI
            call mpi_comm_rank(mpi_comm_world, imageID, ierrMPI)
            RFN%path = RFN%dir // RFN%key // '_' // DT%date // '_' // DT%time // '_process_' // num2str(imageID+1)      // '_' // num2str(counter) // RFN%ext
            end block
#else
            RFN%path = RFN%dir // RFN%key // '_' // DT%date // '_' // DT%time // '_process_' // num2str(1_IK)           // '_' // num2str(counter) // RFN%ext
#endif
            inquire(file=RFN%path,exist=fileExists,iostat=RFN%Err%stat)    ! check if the file already exists
            if (RFN%Err%stat/=0) then
                RFN%Err%occurred = .true.
                RFN%Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file = '" // RFN%path // "'."
                RFN%path = ""
                return
            end if
            if (counter>1000) then
                RFN%Err%occurred = .true.
                RFN%Err%msg = PROCEDURE_NAME//": Unbelievable! "//num2str(counter)//" filenames were tested and all seem to exist."
                RFN%path = ""
                return
            end if
            if (fileExists) cycle
            exit

        end do

    end function getRandomFileName

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Return the value of the requested input environmental variable.
    !>
    !> \param[in]   name    :   The requested environmental variable name.
    !> \param[out]  value   :   The value of the requested environmental variable name.
    !> \param[out]  length  :   The length of the value of the requested environmental variable name.
    !> \param[out]  Err     :   An object of class [Err_type](@ref err_mod::err_type)
    !!                          indicating whether any error has occurred during information collection.
    subroutine getEnvVar(name,value,length,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getEnvVar
#endif
        use Constants_mod, only: IK, MAX_REC_LEN
        use Err_mod, only: Err_type
        implicit none
        character(*), intent(in)                    :: name
        character(:), allocatable, intent(out)      :: value
        integer(IK) , intent(out), optional         :: length
        type(Err_type), intent(out), optional       :: Err

        character(*), parameter                     :: PROCEDURE_NAME = MODULE_NAME // "@getEnvVar()"
        allocate( character(MAX_REC_LEN) :: value )

        Err%occurred = .false.

        if (present(Err)) then
            if (len_trim(adjustl(name))==0) then
                Err%occurred = .true.
                Err%msg =   PROCEDURE_NAME // ": The input environment variable must have a non-zero length."
                return
            end if
            call get_environment_variable(name=name,value=value,length=length,status=Err%stat)
            if (Err%stat==2) then
                Err%occurred = .true.
                Err%msg =   PROCEDURE_NAME // ": Error occurred while fetching the value of the environment variable " // &
                            name // ". The processor does not support environment variables."
                return
            elseif (Err%stat>2) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME//": Unknown error occurred while fetching the value of the environment variable "//name//"."
                return
            end if
        else
            call get_environment_variable(name=name,value=value,length=length)
        end if

        value = trim(adjustl(value))

    end subroutine getEnvVar

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> The [SysCmd_type](@ref syscmd_type) class constructor.
    !> Execute the input system command `cmd` and return.
    !>
    !> \param[in]   cmd : The requested input system command to be executed.
    !> \param[in]   wait : A logical value indicating whether the program should wait for the control to be returned to it by the terminal.
    !>
    !> \return
    !> `SysCmd` : An object of class [SysCmd_type](@ref syscmd_type) containing the attributes and the statistics of the system command execution.
    function constructSysCmd(cmd,wait) result(SysCmd)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructSysCmd
#endif
        implicit none
        character(*), intent(in)        :: cmd
        logical, intent(in), optional   :: wait
        type(SysCmd_type)               :: SysCmd
        SysCmd%cmd = cmd
        SysCmd%exitstat = -huge(0)
        if (present(wait)) then
            SysCmd%wait = wait
        else
            SysCmd%wait = .true.
        end if
        call SysCmd%run()
    end function constructSysCmd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> A method of the [SysCmd_type](@ref syscmd_type) class.
    !> Execute the requested system command and return.
    !>
    !> \param[inout] SysCmd : An object of class [SysCmd_type](@ref syscmd_type) containing the attributes and
    !!                        the statistics of the system command execution.
    subroutine runSysCmd(SysCmd)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: runSysCmd
#endif
        use Constants_mod, only: MAX_REC_LEN
        implicit none
        class(SysCmd_type), intent(inout)   :: SysCmd
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@runSysCmd()"
        if (allocated(SysCmd%Err%msg)) deallocate(SysCmd%Err%msg)
        allocate( character(MAX_REC_LEN) :: SysCmd%Err%msg )
        call execute_command_line   ( SysCmd%cmd &
                                    , wait=SysCmd%wait &
                                    , exitstat=SysCmd%exitstat &
                                    , cmdstat=SysCmd%Err%stat &
                                    , cmdmsg=SysCmd%Err%msg )
        if (SysCmd%Err%stat==0) then
            SysCmd%Err%occurred = .false.
            return
        elseif (SysCmd%Err%stat==-1) then
            SysCmd%Err%occurred = .true.
            SysCmd%Err%msg =    PROCEDURE_NAME // &
                                ": Error occurred. The processor does not support command execution of the command: " // SysCmd%cmd
            return
        elseif (SysCmd%Err%stat==-2 .and. SysCmd%wait) then
            SysCmd%Err%occurred = .true.
            SysCmd%Err%msg =    PROCEDURE_NAME // &
                                ": Error occurred. The processor had to wait for the execution of the command: " // &
                                SysCmd%cmd // ", but the processor does not support asynchronous command execution."
            return
        elseif (SysCmd%Err%stat>0 .and. SysCmd%wait) then
            SysCmd%Err%occurred = .true.
            SysCmd%Err%msg =    PROCEDURE_NAME // &
                                ": Unknown error occurred while attempting to execute the command: " // SysCmd%cmd // &
                                ". The compiler/processor's explanatory message: " // trim(adjustl(SysCmd%Err%msg))
            return
        end if
    end subroutine runSysCmd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Execute the input system command `cmd` and return.
    !>
    !> \param[in]       command     :   The command to executed in the terminal.
    !> \param[in]       wait        :   A logical argument indicating whether the program should wait until the control is
    !!                                  returned to it or should not wait (optional, default = `.true.`).
    !> \param[inout]    exitstat    :   An integer indicating the exit status flag upon exiting the terminal.
    !> \param[out]      Err         :   An object of class [Err_type](@ref err_mod::err_type)
    !!                                  indicating whether any error has occurred during information collection.
    !>
    !> \remark
    !> This is the procedural implementation of the object-oriented [runSysCmd](@ref runsyscmd) method,
    !! kept here only for legacy usage.
    subroutine executeCmd(command,wait,exitstat,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: executeCmd
#endif
        use Constants_mod, only: MAX_REC_LEN
        use Err_mod, only: Err_type
        implicit none
        character(*), intent(in)                :: command
        logical     , intent(in)    , optional  :: wait
        integer     , intent(inout) , optional  :: exitstat
        type(Err_type), intent(out) , optional  :: Err

        logical                                 :: waitDefault
        integer                                 :: exitstatDefault

        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME // "@executeCmd()"

        if (present(wait)) then
            waitDefault = wait
        else
            waitDefault = .true.
        end if

        if (present(exitstat)) then
            exitstatDefault = exitstat
        else
            exitstatDefault = -huge(0_IK)
        end if

        if (present(Err)) then

            Err%occurred = .false.
            allocate( character(MAX_REC_LEN) :: Err%msg )

            call execute_command_line( command, wait=waitDefault, exitstat=exitstatDefault, cmdstat=Err%stat, cmdmsg=Err%msg )
            if (Err%stat==0_IK) then
                return
            elseif (Err%stat==-1_IK) then
                Err%occurred = .true.
                Err%msg =   PROCEDURE_NAME // &
                            ": Error occurred. The processor does not support command execution of the command: " // command
                return
            elseif (Err%stat==-2_IK .and. waitDefault) then
                Err%occurred = .true.
                Err%msg =   PROCEDURE_NAME // ": Error occurred. The processor had to wait for the execution of the command: " // &
                            command // ", but the processor does not support asynchronous command execution."
                return
            elseif (Err%stat>0_IK .and. waitDefault) then
                Err%occurred = .true.
                Err%msg =   PROCEDURE_NAME // ": Unknown error occurred while attempting to execute the command: " // command // &
                            ". The compiler/processor's explanatory message: " // trim(adjustl(Err%msg))
                return
            end if

        else

            call execute_command_line( command, wait=waitDefault, exitstat=exitstatDefault )
            return

        end if

    end subroutine executeCmd

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Fetch the input command-line arguments to the main program.
    !>
    !> \param[inout]    CmdArg : An object of class [CmdArg_type](@ref cmdarg_type) which will contain the command line arguments.
    !>
    !> \remark
    !> This is a method of the class [CmdArg_type](@ref cmdarg_type).
    subroutine queryCmdArg(CmdArg)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: queryCmdArg
#endif
        use String_mod, only: num2str
        use Constants_mod, only: IK, MAX_REC_LEN
        use Err_mod, only: Err_type
        implicit none
        class(CmdArg_type), intent(inout)   :: CmdArg

        integer                             :: i
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@queryCmdArg()"

        CmdArg%Err%occurred = .false.
        CmdArg%Err%msg = ""

        ! first get the full command line
        allocate( character(MAX_REC_LEN) :: CmdArg%cmd )
        call get_command( command=CmdArg%cmd , status = CmdArg%Err%stat )
        if (CmdArg%Err%stat>0) then
            CmdArg%Err%occurred = .true.
            CmdArg%Err%msg = PROCEDURE_NAME // ": Error occurred while fetching the command line."
            return
        elseif (CmdArg%Err%stat==-1) then
            CmdArg%Err%occurred = .true.
            CmdArg%Err%msg = PROCEDURE_NAME // ": Unbelievable error occurred while fetching the command line: &
                           & The length of the command line is longer than " // num2str(MAX_REC_LEN) // "!"
            return
        else
            CmdArg%cmd = trim(adjustl(CmdArg%cmd))
        end if

        ! Now get the command line arguments count
        CmdArg%count = command_argument_count()

        ! Now get the individual command line arguments
        allocate( CmdArg%Arg( 0:CmdArg%count ) )
        do i = 0, CmdArg%count
            allocate( character(MAX_REC_LEN) :: CmdArg%Arg(i)%record )
            call get_command_argument( number=i, value=CmdArg%Arg(i)%record, status=CmdArg%Err%stat )
            if (CmdArg%Err%stat>0) then
                CmdArg%Err%occurred = .true.
                CmdArg%Err%msg = PROCEDURE_NAME // ": Error occurred while fetching the command line."
                return
            elseif (CmdArg%Err%stat==-1) then
                CmdArg%Err%occurred = .true.
                CmdArg%Err%msg = PROCEDURE_NAME // ": Unbelievable error occurred while fetching the command line: &
                               & The length of the command line argument is longer than " // num2str(MAX_REC_LEN) // "!"
                return
            else
                CmdArg%Arg(i)%record = trim(adjustl(CmdArg%Arg(i)%record))
            end if
        end do

    end subroutine queryCmdArg

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Fetch a comprehensive report of the operating system and platform specifications.
    !>
    !> \param[out]  List    :   A list of strings each of which represents one line of information about the system specs.
    !> \param[out]  Err     :   An object of class [Err_type](@ref err_mod::err_type)
    !!                          indicating whether any error has occurred during information collection.
    !> \param[in]   OS      :   An object of class [OS_type](@ref os_type) containing information about the Operating System (optional).
    !> \param[out]  count   :   The count of elements in the output `List` (optional).
    subroutine getSystemInfo(List,Err,OS,count)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getSystemInfo
#endif
        use Err_mod, only: Err_type
        use String_mod, only: num2str
        use Constants_mod, only: IK, RK, MAX_REC_LEN
        use DateTime_mod, only: DateTime_type
        use JaggedArray_mod, only: CharVec_type
        implicit none
        type(CharVec_type)  , intent(out), allocatable  :: List(:)
        type(Err_type)      , intent(out)               :: Err
        type(OS_type)       , intent(in) , optional     :: OS
        integer(IK)         , intent(out), optional     :: count

        type(OS_type)                                   :: OpSy
        character(len=:), allocatable                   :: command,filename,stdErr
        character(len=MAX_REC_LEN)                      :: record
        integer(IK)                                     :: fileUnit,counter,nRecord
        logical                                         :: fileIsOpen

        character(*), parameter                         :: PROCEDURE_NAME = MODULE_NAME // "@getSystemInfo()"

        Err%occurred = .false.
        Err%msg = ""

        ! generate a brand new, non-existing filename
        block
            type(RandomFileName_type)   :: RFN
            RFN = RandomFileName_type(key="getFileList")
            if (RFN%Err%occurred) then
                RFN%Err%msg = PROCEDURE_NAME // RFN%Err%msg
                return
            end if
            filename = RFN%path
        end block
        stdErr = filename // ".stderr"

        ! determine the operating system
        if (present(OS)) then
            OpSy = OS
        else
            call OpSy%query()
            if (OpSy%Err%occurred) then
                Err = OpSy%Err
                Err%msg = PROCEDURE_NAME // Err%msg
                return
            end if
        end if

        if (OpSy%isWindows) then
            command = "systeminfo > " // filename
        elseif (OpSy%isDarwin) then
            command = "uname -a >> " // filename // "; sysctl -a | grep machdep.cpu >> " // filename
        elseif (OpSy%isLinux) then
            !command = "uname -a >> " // filename // "; lshw -short >> " // filename // "; lscpu >> " // filename
            command = "uname -a >> " // filename // "; lscpu >> " // filename
        else ! unknown operating system
            allocate(List(1))
            List(1)%record = "Unknown operating system: " // OpSy%name
            if (present(count)) count = 1_IK
            return
        end if

        call executeCmd( command=command // " 2> " // stdErr, Err=Err )
        if (Err%occurred) then
            Err%msg =   PROCEDURE_NAME // ": Error occurred while attempting to write the system info to external file." // NLC // Err%msg
            return
        end if

        ! now count the number of records in file:

        inquire(file=filename,opened=fileIsOpen,number=fileUnit,iostat=Err%stat)    ! check if the file already exists
        if (Err%stat/=0) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the open status of file = '" // filename // "'."
            return
        end if

        if (fileIsOpen) close(fileUnit,iostat=Err%stat)
        if (Err%stat/=0) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred while attempting to close the open file = '" // filename // "'."
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
            Err%msg = PROCEDURE_NAME // ": Unknown error occurred while opening file = '" // filename // "'."
            return
        end if

        nRecord = 0 ! number of filenames in the file
        do
            read(fileUnit,'(A)',iostat=Err%stat) record
            if ( is_iostat_eor(Err%stat) ) then
                Err%occurred = .true.
                Err%msg  = PROCEDURE_NAME // ": End-Of-Record error condition occurred while attempting to read &
                         & from file = '" // filename // "'."
                return
            elseif ( is_iostat_end(Err%stat) ) then
                exit
            elseif ( Err%stat>0 ) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Unknown error condition occurred while attempting to read &
                        & from file = '" // filename // "'."
                return
            else
                nRecord = nRecord + 1
                cycle
            end if
        end do
        close(fileUnit,iostat=Err%stat)
        if (Err%stat/=0) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred while attempting to close the open file = '" // filename // "'."
            return
        end if

        ! now read the contents of the file

        call sleep(seconds=0.1_RK,Err=Err)
        if (Err%occurred) then
            Err%msg = PROCEDURE_NAME // Err%msg
            return
        end if

        open(newunit=fileUnit,file=filename,status="old",iostat=Err%stat)
        if (Err%stat>0) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Unknown error occurred while opening file = '" // filename // "'."
            return
        end if

        allocate(List(nRecord))
        do counter = 1,nRecord
            read(fileUnit,'(A)',iostat=Err%stat) record
            if ( is_iostat_eor(Err%stat) ) then
                Err%occurred = .true.
                Err%msg  = PROCEDURE_NAME // ": End-Of-Record error condition occurred while attempting to read &
                         & from file = '" // filename // "'."
                return
            elseif ( is_iostat_end(Err%stat) ) then
                exit
            elseif ( Err%stat>0 ) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Unknown error condition occurred while attempting to read &
                        & from file = '" // filename // "'."
                return
            end if
            List(counter)%record = trim(adjustl(record))
        end do

        close(fileUnit,iostat=Err%stat,status="delete")
        if (Err%stat/=0) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred while attempting to close the open file = '" // filename // "'."
            return
        end if

        if (present(count)) count = nRecord

        ! remove the files
        !call removeFile(filename,OpSy%isWindows,Err)
        !if (Err%occurred) then
        !    Err%msg = PROCEDURE_NAME // Err%msg
        !    return
        !end if
        !call removeFile(stdErr,OpSy%isWindows,Err)
        !if (Err%occurred) then
        !    Err%msg = PROCEDURE_NAME // Err%msg
        !    return
        !end if

    end subroutine getSystemInfo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Sleep for the input number of seconds (real number).
    !>
    !> \param[in]   seconds :   The amount of time in seconds to sleep.
    !> \param[out]  Err     :   An object of class [Err_type](@ref err_mod::err_type)
    !!                          indicating whether any error has occurred before, during, or after the sleep.
    subroutine sleep(seconds,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: sleep
#endif

        use, intrinsic :: iso_fortran_env, only: int64
        use Err_mod, only: Err_type
        use Constants_mod, only: RK
        implicit none

        real(RK), intent(in)            :: seconds ! sleep time
        type(Err_type) , intent(out)    :: Err

        integer(int64)                  :: countOld, countNew, countMax
        real(RK)                        :: countRate

        character(*), parameter         :: PROCEDURE_NAME = MODULE_NAME // "@sleep()"

        Err%occurred = .false.
        Err%msg = ""

        call system_clock( count=countOld, count_rate=countRate, count_max=countMax )
        if (countOld==-huge(0) .or. nint(countRate)==0 .or. countMax==0) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred. There is no processor clock."
            return
        end if

        countRate = 1._RK / countRate
        do
            call system_clock( count=countNew )
            if (countNew==countMax) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred. Maximum processor clock count reached."
            end if
            if ( real(countNew-countOld,kind=RK) * countRate > seconds ) exit
            cycle
        end do

    end subroutine sleep

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Copy file from the origin path to the destination path.
    !>
    !> \param[in]   pathOld     :   The original path.
    !> \param[in]   pathNew     :   The destination path.
    !> \param[in]   isWindows   :   Logical value indicating whether the OS and the terminal is Windows shell (CMD or Powershell).
    !> \param[out]  Err         :   An object of class [Err_type](@ref err_mod::err_type)
    !!                              indicating whether any error has occurred the copy.
    subroutine copyFile(pathOld,pathNew,isWindows,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: copyFile
#endif

        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        character(*), intent(in)    :: pathOld, pathNew
        logical     , intent(in)    :: isWindows
        type(Err_type), intent(out) :: Err
        character(:), allocatable   :: cmd
        integer                     :: counter
        logical                     :: fileExists

        character(*), parameter     :: PROCEDURE_NAME = MODULE_NAME // "@copyFile()"

        Err%occurred = .false.

        if (len_trim(adjustl(pathOld))==0) return

        ! First check whether file exists:
        inquire(file=pathNew,exist=fileExists,iostat=Err%stat)    ! check if the file already exists
        if (Err%stat/=0) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file = '" // pathNew // "'."
            return
        end if
        if (fileExists) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": The requested copy file = '" // pathNew // "' already exists."
            return
        end if

        if (isWindows) then
            cmd = 'copy "'  // pathOld // '" "' // pathNew // '" > nul'
        else
            cmd = "cp "     // pathOld // " " // pathNew
        end if

        counter = 0
        do
            counter = counter + 1
            call executeCmd( command=cmd, Err=Err )
            if (Err%occurred) then
                Err%msg = PROCEDURE_NAME // ": Error occurred while executing command "// cmd // "'." // NLC
                return
            end if
            ! ensure file is copied
            inquire(file=pathNew,exist=fileExists,iostat=Err%stat)    ! check if the file already exists
            if (Err%stat/=0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of copied file = '" // pathNew // "'."
                return
            end if
            if (.not. fileExists .and. counter<100) cycle
            exit
        end do
        if (.not. fileExists) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Failed to copy file from '" // pathOld // "' to '" // pathNew // "' after " // num2str(counter) // " attempts."
            return
        end if

    end subroutine copyFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> Remove the requested file.
    !>
    !> \param[in]   path        :   The path to the file to be removed.
    !> \param[in]   isWindows   :   Logical value indicating whether the OS is Windows.
    !> \param[out]  Err         :   An object of class [Err_type](@ref err_mod::err_type)
    !!                              indicating whether any error has occurred before, during, or after the sleep.
    !>
    !> \warning
    !> This subroutine can become extremely dangerous if one does understands the
    !! scopes of the removal of the requested file or pattern. **Use with caution**.
    subroutine removeFile(path,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: removeFile
#endif

        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        character(*), intent(in)                :: path
        type(Err_type), intent(out), optional   :: Err
        !logical     , intent(in), optional      :: isWindows
        logical                                 :: fileExists
        logical                                 :: isPresentErr
        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME // "@removeFile()"

        isPresentErr = present(Err)

        ! First check whether file exists.

        if (isPresentErr) then
            Err%occurred = .false.
            inquire(file=path,exist=fileExists,iostat=Err%stat)    ! check if the file already exists
            if (Err%stat/=0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file = '" // path // "'."
                return
            end if
        else
            inquire(file=path,exist=fileExists) ! check if the file already exists
        end if

        if (.not. fileExists) return

        !if (isPresentErr .and. present(isWindows)) then
        !
        !    blockBrittle: block
        !
        !        character(:), allocatable               :: cmd
        !        integer                                 :: counter
        !
        !        if (isWindows) then
        !            cmd = "del " // path // " > nul"
        !        else
        !            cmd = "rm " // path
        !        end if
        !
        !        counter = 0
        !        do
        !            counter = counter + 1
        !            call executeCmd( command=cmd, Err=Err )
        !            if (Err%occurred) then
        !                Err%msg = PROCEDURE_NAME // ": Error occurred while executing command "// cmd // "'." // NLC
        !                return
        !            end if
        !            ! ensure file is removed
        !            inquire(file=path,exist=fileExists,iostat=Err%stat)    ! check if the file already exists
        !            if (Err%stat/=0) then
        !                Err%occurred = .true.
        !                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of removed file = '" // path // "'."
        !                return
        !            end if
        !            if (fileExists .and. counter<100) cycle
        !            exit
        !        end do
        !        if (fileExists) then
        !            Err%occurred = .true.
        !            Err%msg = PROCEDURE_NAME // ": Failed to remove file = '" // path // "' after " // num2str(counter) // " attempts."
        !            return
        !        end if
        !
        !    end block blockBrittle
        !
        !else

            blockRobust: block
                logical :: isOpen
                integer :: fileUnit
                inquire(file=path,opened=isOpen)
                if (.not. isOpen) open(newunit = fileUnit, file = path, status = "replace")
                if (isPresentErr) then
                    close(fileUnit, status="delete", iostat = Err%stat)
                    Err%occurred = Err%stat > 0_IK
                else
                    close(fileUnit,status="delete")
                end if
            end block blockRobust

        !end if

    end subroutine removeFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module System_mod