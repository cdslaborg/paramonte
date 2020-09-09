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

module System_mod

    use JaggedArray_mod, only: CharVec_type
    use Constants_mod, only: IK, NLC
    use Err_mod, only: Err_type
    implicit none

    character(*), parameter             :: MODULE_NAME = "@System_mod"
    integer(IK) , parameter             :: MAX_OS_NAME_LEN = 63_IK

    type :: RandomFileName_type
        character(:), allocatable       :: path
        character(:), allocatable       :: dir
        character(:), allocatable       :: key
        character(:), allocatable       :: ext
        type(Err_type)                  :: Err
    end type RandomFileName_type

    interface RandomFileName_type
        module procedure                :: getRandomFileName
    end interface RandomFileName_type

    type :: SystemInfo_type
        integer(IK)                     :: nRecord
        type(CharVec_type), allocatable :: List(:)
        type(Err_type)                  :: Err
    contains
        procedure, nopass               :: get => getSystemInfo
    end type SystemInfo_type

    interface SystemInfo_type
        module procedure                :: constructSystemInfo
    end interface SystemInfo_type

    type :: OS_type
        character(:), allocatable       :: name, slash
        logical                         :: isWindows = .false.
        logical                         :: isDarwin = .false.
        logical                         :: isLinux = .false.
        type(Err_type)                  :: Err
    contains
        procedure, pass                 :: query => queryOS
    end type OS_type

    type :: EnvVar_type
        character(:), allocatable       :: name
        character(:), allocatable       :: value
        integer                         :: length
        type(Err_type)                  :: Err
    contains
        procedure, nopass               :: get => getEnvVar
    end type EnvVar_type

    type :: CmdArg_type
        character(:), allocatable       :: cmd
        type(CharVec_type), allocatable :: Arg(:)
        integer                         :: count
        type(Err_type)                  :: Err
    contains
        procedure, pass                 :: query => queryCmdArg
    end type CmdArg_type

    type :: SysCmd_type
        character(:), allocatable       :: cmd
        logical                         :: wait
        integer                         :: exitstat
        type(Err_type)                  :: Err
    contains
        procedure, pass                 :: run => runSysCmd
    end type SysCmd_type

    interface SysCmd_type
        module procedure :: constructSysCmd
    end interface SysCmd_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    ! Queries all of the OS class attributes: name, slash, isWindows, slash, Err
    subroutine queryOS(OS)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: queryOS
#endif

        use String_mod, only: num2str, getLowerCase
        use Constants_mod, only: IK, RK
        use Err_mod, only: Err_type
        implicit none
        class(OS_type)  , intent(out)   :: OS
        character(*)    , parameter     :: PROCEDURE_NAME = MODULE_NAME // "@queryOS()"
        character(:)    , allocatable   :: osname

        OS%Err%occurred = .false.
        OS%Err%msg = ""

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

                    integer                     :: idummy
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

                    open(newunit=idummy,file=RFN%path,status="old",iostat=OS%Err%stat)
                    if (OS%Err%stat>0) then
                        OS%Err%occurred = .true.
                        OS%Err%msg =    PROCEDURE_NAME // ": Unknown error occurred while opening file = '" // RFN%path // "'."
                        OS%name = ""
                        return
                    end if

                    read(idummy,*,iostat=OS%Err%stat) OS%name

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

                    close(idummy)
                    if (OS%Err%stat>0) then
                        OS%Err%occurred = .true.
                        OS%Err%msg = PROCEDURE_NAME // ": Unknown error occurred while closing file = '" // RFN%path // "'."
                        OS%name = ""
                        return
                    end if

                    call sleep( seconds=0.5_RK, Err=OS%Err )
                    if (OS%Err%occurred) then
                        OS%Err%msg = PROCEDURE_NAME // ": Error occurred while calling subroutine sleep() for querying file = '" // RFN%path // "'." // NLC // OS%Err%msg
                        OS%name = ""
                        return
                    end if

                    call removeFile( path=RFN%path, isWindows=.false., Err=OS%Err )
                    if (OS%Err%occurred) then
                        OS%Err%msg = PROCEDURE_NAME // ": Error occurred while removing file = '"// RFN%path // "'." // NLC // OS%Err%msg
                        OS%name = ""
                        return
                    end if

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

    end subroutine queryOS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! generates unique file path in the requested directory for temporary usage
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

    subroutine getEnvVar(name,value,length,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getEnvVar
#endif
        use Constants_mod, only: IK, MAX_REC_LEN
        use Err_mod, only: Err_type
        implicit none
        character(*), intent(in)                    :: name
        character(:), allocatable, intent(out)      :: value
        integer(IK) , optional, intent(out)         :: length
        type(Err_type), optional, intent(out)       :: Err

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

    ! This is the same as runSysCmd, but procedural, rahter than OOP, kept only for legacy usage.
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

        close(fileUnit,iostat=Err%stat)
        if (Err%stat/=0) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred while attempting to close the open file = '" // filename // "'."
            return
        end if

        if (present(count)) count = nRecord

        ! remove the files
        call removeFile(filename,OpSy%isWindows,Err)
        if (Err%occurred) then
            Err%msg = PROCEDURE_NAME // Err%msg
            return
        end if
        call removeFile(stdErr,OpSy%isWindows,Err)
        if (Err%occurred) then
            Err%msg = PROCEDURE_NAME // Err%msg
            return
        end if

    end subroutine getSystemInfo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    subroutine removeFile(path,isWindows,Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: removeFile
#endif

        use Err_mod, only: Err_type
        use String_mod, only: num2str
        implicit none
        character(*), intent(in)    :: path
        logical     , intent(in)    :: isWindows
        type(Err_type), intent(out) :: Err
        character(:), allocatable   :: cmd
        integer                     :: counter
        logical                     :: fileExists

        character(*), parameter     :: PROCEDURE_NAME = MODULE_NAME // "@removeFile()"

        Err%occurred = .false.

        ! First check whether file exists:
        inquire(file=path,exist=fileExists,iostat=Err%stat)    ! check if the file already exists
        if (Err%stat/=0) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file = '" // path // "'."
            return
        end if
        if (.not.fileExists) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": The requested file = '" // path // "' does not exist."
            return
        end if

        if (isWindows) then
            cmd = "del " // path // " > nul"
        else
            cmd = "rm " // path
        end if

        counter = 0
        do
            counter = counter + 1
            call executeCmd( command=cmd, Err=Err )
            if (Err%occurred) then
                Err%msg = PROCEDURE_NAME // ": Error occurred while executing command "// cmd // "'." // NLC
                return
            end if
            ! ensure file is removed
            inquire(file=path,exist=fileExists,iostat=Err%stat)    ! check if the file already exists
            if (Err%stat/=0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of removed file = '" // path // "'."
                return
            end if
            if (fileExists .and. counter<100) cycle
            exit
        end do
        if (fileExists) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Failed to remove file = '" // path // "' after " // num2str(counter) // " attempts."
            return
        end if

    end subroutine removeFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module System_mod