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

module File_mod
    
    use Path_mod, only: Path_type
    use Err_mod, only: Err_type

    implicit none

    character(*), parameter         :: MODULE_NAME = "@File_mod"

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type :: Action_type
        character(:), allocatable   :: value            ! = read, write, readwrite, undefined. Default is processor-dependent.
        logical                     :: isRead           = .false.
        logical                     :: isWrite          = .false.
        logical                     :: isReadWrite      = .false.
        logical                     :: isUndefined      = .false.
        type(Err_type)              :: Err
    end type Action_type

    interface Action_type
        module procedure :: constructAction
    end interface

    type :: Access_type
        character(:), allocatable   :: value            ! = sequential (default), direct, undefined
        logical                     :: isSequential     = .false.
        logical                     :: isDirect         = .false.
        logical                     :: isUndefined      = .false.
        type(Err_type)              :: Err
    end type Access_type

    interface Access_type
        module procedure :: constructAccess
    end interface

    type :: Form_type
        character(:), allocatable   :: value            ! = formatted, unformatted (default depends on ACCESS), undefined.
        logical                     :: isFormatted      = .false.
        logical                     :: isUnformatted    = .false.
        logical                     :: isUndefined      = .false.
        type(Err_type)              :: Err
    end type Form_type

    interface Form_type
        module procedure :: constructForm
    end interface

    type :: Blank_type
        character(:), allocatable   :: value            ! = null (default), zero, undefined.
        logical                     :: isNull           = .false.
        logical                     :: isZero           = .false.
        logical                     :: isUndefined      = .false.
        type(Err_type)              :: Err
    end type Blank_type

    interface Blank_type
        module procedure :: constructBlank
    end interface

    type :: Position_type
        character(:), allocatable   :: value            ! = asis (default), rewind, append, undefined. For ACCESS=sequential.
        logical                     :: isAsis           = .false.
        logical                     :: isRewind         = .false.
        logical                     :: isAppend         = .false.
        logical                     :: isUndefined      = .false.
        type(Err_type)              :: Err
    end type Position_type

    interface Position_type
        module procedure :: constructPosition
    end interface

    type :: Delim_type
        character(:), allocatable   :: value            ! = quote, apostrophe, undefined, or none (default).
        logical                     :: isQuote          = .false.
        logical                     :: isApostrophe     = .false.
        logical                     :: isNone           = .false.
        logical                     :: isUndefined      = .false.
        type(Err_type)              :: Err
    end type Delim_type

    interface Delim_type
        module procedure :: constructDelim
    end interface

    type :: Pad_type
        character(:), allocatable   :: value            ! = yes (default), no, undefined.
        logical                     :: isPadded         = .false.
        logical                     :: isNotPadded      = .false.
        logical                     :: isUndefined      = .false.
        type(Err_type)              :: Err
    end type Pad_type

    interface Pad_type
        module procedure :: constructPad
    end interface

    type :: Round_type
        character(:), allocatable   :: value            ! = up, down, zero, nearest, compatible, processor_defined, or undefined.
        logical                     :: isUp             = .false.
        logical                     :: isDown           = .false.
        logical                     :: isZero           = .false.
        logical                     :: isNearest        = .false.
        logical                     :: isCompatible     = .false.
        logical                     :: isProcessDefined = .false.
        logical                     :: isUndefined      = .false.
        type(Err_type)              :: Err
    end type Round_type

    interface Round_type
        module procedure :: constructRound
    end interface

    type :: Sign_type
        character(:), allocatable   :: value            ! = suppress, plus, processor_defined, or undefined.
        logical                     :: isSuppress       = .false.
        logical                     :: isPlus           = .false.
        logical                     :: isProcessDefined = .false.
        logical                     :: isUndefined      = .false.
        type(Err_type)              :: Err
    end type Sign_type

    interface Sign_type
        module procedure :: constructSign
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type :: File_type
        integer                     :: unit         = -huge(0)
        integer                     :: number       = -huge(0)
        integer                     :: recl         = -huge(0)
        logical                     :: exists       = .false.
        logical                     :: isOpen       = .false.
        logical                     :: isNamed      = .false.
        logical                     :: isInternal   = .false.
        logical                     :: isNumbered   = .false.
        character(:), allocatable   :: status       ! = old, new, replace, scratch, unknown (default).
        character(:), allocatable   :: asynchronous ! = yes, no
        character(:), allocatable   :: format       ! the specific content format statement to be used with read/write statements.
        character(:), allocatable   :: nameByCompiler
        type(Action_type)           :: Action
        type(Access_type)           :: Access
        type(Form_type)             :: Form
        type(Blank_type)            :: Blank
        type(Position_type)         :: Position
        type(Delim_type)            :: Delim
        type(Pad_type)              :: Pad
        type(Round_type)            :: Round
        type(Sign_type)             :: Sign
        type(Path_type)             :: Path
        type(Err_type)              :: Err
    contains
        procedure, pass             :: openFile
        procedure, pass             :: closeFile
        procedure, nopass           :: getNumber
        procedure, nopass           :: getPosition
        procedure, nopass           :: getAction
        procedure, nopass           :: getDelim
        procedure, nopass           :: getRecl
        procedure, nopass           :: getBlank
        procedure, nopass           :: getOpenStatus
        procedure, nopass           :: getExistStatus
        procedure, nopass           :: getInqErr
        procedure, nopass           :: getReadErr
        procedure, nopass           :: getOpenErr
        procedure, nopass           :: getCloseErr
        procedure, nopass           :: getWriteErr
    end type File_type

    interface File_type
        module procedure :: constructFile
    end interface

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructFile( unit, recl, path, status, position, access, form, action, delim &
                          , round, sign,pad, blank, format, asynchronous &
                          , OS &
                          ) result(File)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructFile
#endif

        use String_mod, only: getLowerCase
        use System_mod, only: OS_type
        implicit none
        type(File_type) :: File
        integer     , intent(in), optional  :: unit
        integer     , intent(in), optional  :: recl
        character(*), intent(in), optional  :: status
        character(*), intent(in), optional  :: asynchronous
        character(*), intent(in), optional  :: access
        character(*), intent(in), optional  :: position
        character(*), intent(in), optional  :: form
        character(*), intent(in), optional  :: action
        character(*), intent(in), optional  :: delim
        character(*), intent(in), optional  :: round
        character(*), intent(in), optional  :: sign
        character(*), intent(in), optional  :: pad
        character(*), intent(in), optional  :: blank
        character(*), intent(in), optional  :: path
        character(*), intent(in), optional  :: format
        type(OS_type), intent(in), optional :: OS

        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME // "@constructFile()"

        File%Err%occurred = .false.
        File%Err%stat = -huge(0)
        File%Err%msg = ""

        if (present(unit)) then
            File%unit = unit
        else
            File%unit = -huge(0)
        end if

        if (present(recl)) then
            File%recl = recl
        else
            File%recl = -huge(0)
        end if

        ! set up file path

        if (present(path)) then
!write(*,*) OS%slash
!write(*,*) OS%isWindows
!write(*,*) path
            File%Path = path_type(inputPath=path,OS=OS)
        else
            File%Path = path_type(inputPath="",OS=OS)
        end if
!write(*,*) File%Path%original
!write(*,*) File%Path%modified
        if (File%Path%Err%occurred) then
            File%Err = File%Path%Err
            File%Err%msg = PROCEDURE_NAME // File%Err%msg
            return
        end if

        ! check if file exists

        call File%getExistStatus( exists = File%exists      &
                                , Err = File%err            &
                                , file = File%Path%modified )
        if (File%Err%occurred) then
            File%Err%msg = PROCEDURE_NAME // File%Err%msg
            return
        end if

         ! if it does not exist, try the original file path

        if (.not.File%exists) then
            call File%getExistStatus( exists = File%exists      &
                                    , Err = File%err            &
                                    , file = File%Path%original )
            if (File%exists) File%Path%modified = File%Path%original    ! restore the original path, which is apparently the correct path
        end if
        if (File%Err%occurred) then
            File%Err%msg = PROCEDURE_NAME // File%Err%msg
            return
        end if

        ! set up the rest of attributes

        if (present(format)) then
            File%format = trim(adjustl(format))
        else
            File%format = ""
        end if

        if (present(status)) then
            File%status = getLowerCase(trim(adjustl(status)))
        else
            File%status = "unknown"
        end if

        if (present(asynchronous)) then
            File%asynchronous = getLowerCase(trim(adjustl(asynchronous)))
        else
            File%asynchronous = "no"
        end if

        File%Action = Action_type(action)
        If (File%Action%Err%occurred) then
            File%Err = File%Action%Err
            File%Err%msg = PROCEDURE_NAME // File%Err%msg
            return
        end if

        File%Delim = Delim_type(delim)
        If (File%Delim%Err%occurred) then
            File%Err = File%Delim%Err
            File%Err%msg = PROCEDURE_NAME // File%Err%msg
            return
        end if

        File%Access = Access_type(access)
        If (File%Access%Err%occurred) then
            File%Err = File%Access%Err
            File%Err%msg = PROCEDURE_NAME // File%Err%msg
            return
        end if

        File%Position = Position_type(position)
        If (File%Position%Err%occurred) then
            File%Err = File%Position%Err
            File%Err%msg = PROCEDURE_NAME // File%Err%msg
            return
        end if
        
        if (present(form)) then
            File%form = Form_type(form)
        else
            if ( File%Access%isDirect ) then
                File%Form = Form_type("unformatted")
            else    ! if ( File%Access%isSequential ) then
                File%Form = Form_type("formatted")
            end if
        end if
        If (File%Form%Err%occurred) then
            File%Err = File%Form%Err
            File%Err%msg = PROCEDURE_NAME // File%Err%msg
            return
        end if

        File%Round = Round_type(round)
        If (File%Round%Err%occurred) then
            File%Err = File%Round%Err
            File%Err%msg = PROCEDURE_NAME // File%Err%msg
            return
        end if

        File%Sign = Sign_type(sign)
        If (File%Sign%Err%occurred) then
            File%Err = File%Sign%Err
            File%Err%msg = PROCEDURE_NAME // File%Err%msg
            return
        end if

        File%Pad = Pad_type(pad)
        If (File%Pad%Err%occurred) then
            File%Err = File%Pad%Err
            File%Err%msg = PROCEDURE_NAME // File%Err%msg
            return
        end if

        File%Blank = Blank_type(blank)
        If (File%Blank%Err%occurred) then
            File%Err = File%Blank%Err
            File%Err%msg = PROCEDURE_NAME // File%Err%msg
            return
        end if

        File%nameByCompiler = ""

    end function constructFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! subroutine inquireFile(Self,unit,file)
        ! use Path_mod, only: MAX_FILE_PATH_LEN
        ! use String_mod, only: num2str
        ! implicit none
        ! class(File_type), intent(inout)         :: Self
        ! integer         , intent(in), optional  :: unit
        ! character(*)    , intent(in), optional  :: file
        ! character(*)    , parameter             :: PROCEDURE_NAME = MODULE_NAME // "@inquireFile()"

        ! Self%Err%msg = ""
        ! Self%Err%occurred = .false.

        ! if (allocated(Self%nameByCompiler)) deallocate(Self%nameByCompiler)
        ! allocate( character(MAX_FILE_PATH_LEN) :: Self%nameByCompiler )

        ! if (allocated(Self%Access%value)) deallocate(Self%access)
        ! allocate( character(63) :: Self%access )
        
        ! if (allocated(Self%form)) deallocate(Self%form)
        ! allocate( character(63) :: Self%form )


        ! if (present(unit)) then
            ! inquire( unit   = unit &
                   ! , exist  = Self%exists &
                   ! , opened =  Self%isOpen &
                   ! , number =  Self%number &
                   ! , named  =  Self%isNamed &
                   ! , name   =  Self%nameByCompiler &
                   ! , access =  Self%access &
                   ! , iostat = Err%stat &
                   ! )
            
            ! if (Err%stat>0) then
                ! Err%occurred = .true.
                ! Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with unit=" // num2str(unit) // "."
                ! return
            ! end if
        ! if (present(file)) then
            ! inquire(file=file,exist=exists,iostat=Err%stat)
            ! if (Err%stat>0) then
                ! Err%occurred = .true.
                ! Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with name=" // file // "."
                ! return
            ! end if
        ! else
            ! Err%occurred = .true.
            ! Err%msg = PROCEDURE_NAME // ": At least one of the two input arguments (unit,path) must be provided."
            ! return
        ! end if
        ! if (Self%number==-1) Self%isNumbered = .false.
        ! Self%nameByCompiler = trim(adjustl(Self%nameByCompiler))
    ! end subroutine inquireFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine getExistStatus(exists,Err,unit,file)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getExistStatus
#endif
        use String_mod, only: num2str
        use Err_mod, only: Err_type
        implicit none
        logical, intent(out)                :: exists
        type(Err_type), intent(out)         :: Err
        integer, intent(in), optional       :: unit
        character(*), intent(in), optional  :: file
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@getExistStatus()"
        Err%msg = ""
        Err%occurred = .false.
        if (present(unit)) then
            inquire(unit=unit,exist=exists,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with unit=" // num2str(unit) // "."
                return
            end if
        elseif (present(file)) then
            inquire(file=file,exist=exists,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with name=" // file // "."
                return
            end if
        elseif (present(unit) .and. present(file)) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Only one of the two optional arguments (unit, file) must be provided as input."
            return
        else
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": At least one of the two input arguments (unit,path) must be provided."
            return
        end if
    end subroutine getExistStatus

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine getOpenStatus(isOpen,Err,unit,file)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getOpenStatus
#endif
        use String_mod, only: num2str
        use Err_mod, only: Err_type
        implicit none
        logical, intent(out)                :: isOpen
        type(Err_type), intent(out)         :: Err
        integer, intent(in), optional       :: unit
        character(*), intent(in), optional  :: file
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@getOpenStatus()"
        Err%msg = ""
        Err%occurred = .false.
        if (present(unit)) then
            inquire(unit=unit,opened=isOpen,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with unit=" // num2str(unit) // "."
                return
            end if
        elseif (present(file)) then
            inquire(file=file,opened=isOpen,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with name=" // file // "."
                return
            end if
        else
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": At least one of the two input arguments (unit,path) must be provided."
            return
        end if
    end subroutine getOpenStatus

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine getNumber(isNumbered,number,Err,unit,file)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getNumber
#endif
        use String_mod, only: num2str
        use Err_mod, only: Err_type
        implicit none
        logical, intent(out)                :: isNumbered
        integer, intent(out)                :: number
        type(Err_type), intent(out)         :: Err
        integer, intent(in), optional       :: unit
        character(*), intent(in), optional  :: file
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@getNumber()"
        Err%msg = ""
        Err%occurred = .false.
        isNumbered = .true.
        if (present(unit)) then
            inquire(unit=unit,number=number,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with unit=" // num2str(unit) // "."
                return
            end if
        elseif (present(file)) then
            inquire(file=file,number=number,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with name=" // file // "."
                return
            end if
        else
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": At least one of the two input arguments (unit,path) must be provided."
            return
        end if
        if (number==-1) isNumbered = .false.
    end subroutine getNumber

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine getName(isNamed,nameByCompiler,Err,unit,file)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getName
#endif
        use String_mod, only: num2str
        use Path_mod, only: MAX_FILE_PATH_LEN
        use Err_mod, only: Err_type
        implicit none
        logical, intent(out)                    :: isNamed
        character(:), allocatable, intent(out)  :: nameByCompiler
        type(Err_type), intent(out)             :: Err
        integer, intent(in), optional           :: unit
        character(*), intent(in), optional      :: file
        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME // "@getName()"
        Err%msg = ""
        Err%occurred = .false.
        allocate( character(MAX_FILE_PATH_LEN) :: nameByCompiler )
        if (present(unit)) then
            inquire(unit=unit,named=isNamed,name=nameByCompiler,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with unit=" // num2str(unit) // "."
                return
            end if
        elseif (present(file)) then
            inquire(file=file,named=isNamed,name=nameByCompiler,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with name=" // file // "."
                return
            end if
        else
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": At least one of the two input arguments (unit,path) must be provided."
            return
        end if
        nameByCompiler = trim(adjustl(nameByCompiler))
    end subroutine getName

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine getAccess(access,Err,unit,file)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getAccess
#endif
        use String_mod, only: num2str, getLowerCase
        use Err_mod, only: Err_type
        implicit none
        character(:), allocatable, intent(out)  :: access
        type(Err_type), intent(out)             :: Err
        integer, intent(in), optional           :: unit
        character(*), intent(in), optional      :: file
        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME // "@getAccess()"
        Err%msg = ""
        Err%occurred = .false.
        allocate( character(63) :: access )
        if (present(unit)) then
            inquire(unit=unit,access=access,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with unit=" // num2str(unit) // "."
                return
            end if
        elseif (present(file)) then
            inquire(file=file,access=access,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with name=" // file // "."
                return
            end if
        else
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": At least one of the two input arguments (unit,path) must be provided."
            return
        end if
        access = getLowerCase( trim(adjustl(access)) )
    end subroutine getAccess

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine getForm(form,Err,unit,file)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getForm
#endif
        use String_mod, only: num2str, getLowerCase
        use Err_mod, only: Err_type
        implicit none
        character(:), allocatable, intent(out)  :: form
        type(Err_type), intent(out)             :: Err
        integer, intent(in), optional           :: unit
        character(*), intent(in), optional      :: file
        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME // "@getForm()"
        Err%msg = ""
        Err%occurred = .false.
        allocate( character(63) :: form )
        if (present(unit)) then
            inquire(unit=unit,form=form,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with unit=" // num2str(unit) // "."
                return
            end if
        elseif (present(file)) then
            inquire(file=file,form=form,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with name=" // file // "."
                return
            end if
        else
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": At least one of the two input arguments (unit,path) must be provided."
            return
        end if
        form = getLowerCase( trim(adjustl(form)) )
    end subroutine getForm

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine getRecl(recl,Err,unit,file)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getRecl
#endif
        use String_mod, only: num2str
        use Err_mod, only: Err_type
        implicit none
        integer, intent(out)                :: recl
        type(Err_type), intent(out)         :: Err
        integer, intent(in), optional       :: unit
        character(*), intent(in), optional  :: file
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@getRecl()"
        Err%msg = ""
        Err%occurred = .false.
        if (present(unit)) then
            inquire(unit=unit,recl=recl,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with unit=" // num2str(unit) // "."
                return
            end if
        elseif (present(file)) then
            inquire(file=file,recl=recl,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with name=" // file // "."
                return
            end if
        else
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": At least one of the two input arguments (unit,path) must be provided."
            return
        end if
    end subroutine getRecl

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine getBlank(blank,Err,unit,file)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getBlank
#endif
        use String_mod, only: num2str, getLowerCase
        use Err_mod, only: Err_type
        implicit none
        character(:), allocatable, intent(out)  :: blank
        type(Err_type), intent(out)             :: Err
        integer, intent(in), optional           :: unit
        character(*), intent(in), optional      :: file
        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME // "@getBlank()"
        Err%msg = ""
        Err%occurred = .false.
        allocate( character(63) :: blank )
        if (present(unit)) then
            inquire(unit=unit,blank=blank,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with unit=" // num2str(unit) // "."
                return
            end if
        elseif (present(file)) then
            inquire(file=file,blank=blank,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with name=" // file // "."
                return
            end if
        else
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": At least one of the two input arguments (unit,path) must be provided."
            return
        end if
        blank = getLowerCase( trim(adjustl(blank)) )
    end subroutine getBlank

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine getPosition(position,Err,unit,file)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getPosition
#endif
        use String_mod, only: num2str, getLowerCase
        use Err_mod, only: Err_type
        implicit none
        character(:), allocatable, intent(out)  :: position
        type(Err_type), intent(out)             :: Err
        integer, intent(in), optional           :: unit
        character(*), intent(in), optional      :: file
        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME // "@getPosition()"
        Err%msg = ""
        Err%occurred = .false.
        allocate( character(63) :: position )
        if (present(unit)) then
            inquire(unit=unit,position=position,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with unit=" // num2str(unit) // "."
                return
            end if
        elseif (present(file)) then
            inquire(file=file,position=position,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with name=" // file // "."
                return
            end if
        else
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": At least one of the two input arguments (unit,path) must be provided."
            return
        end if
        position = getLowerCase( trim(adjustl(position)) )
    end subroutine getPosition

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine getAction(action,Err,unit,file)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getAction
#endif
        use String_mod, only: num2str, getLowerCase
        use Err_mod, only: Err_type
        implicit none
        character(:), allocatable, intent(out)  :: action
        type(Err_type), intent(out)             :: Err
        integer, intent(in), optional           :: unit
        character(*), intent(in), optional      :: file
        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME // "@getAction()"
        Err%msg = ""
        Err%occurred = .false.
        allocate( character(63) :: action )
        if (present(unit)) then
            inquire(unit=unit,action=action,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with unit=" // num2str(unit) // "."
                return
            end if
        elseif (present(file)) then
            inquire(file=file,action=action,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with name=" // file // "."
                return
            end if
        else
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": At least one of the two input arguments (unit,path) must be provided."
            return
        end if
        action = getLowerCase( trim(adjustl(action)) )
    end subroutine getAction

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine getDelim(delim,Err,unit,file)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getDelim
#endif
        use String_mod, only: num2str, getLowerCase
        use Err_mod, only: Err_type
        implicit none
        character(:), allocatable, intent(out)  :: delim
        type(Err_type), intent(out)             :: Err
        integer, intent(in), optional           :: unit
        character(*), intent(in), optional      :: file
        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME // "@getDelim()"
        Err%msg = ""
        Err%occurred = .false.
        allocate( character(63) :: delim )
        if (present(unit)) then
            inquire(unit=unit,delim=delim,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with unit=" // num2str(unit) // "."
                return
            end if
        elseif (present(file)) then
            inquire(file=file,delim=delim,iostat=Err%stat)
            if (Err%stat>0) then
                Err%occurred = .true.
                Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file with name=" // file // "."
                return
            end if
        else
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": At least one of the two input arguments (unit,path) must be provided."
            return
        end if
        delim = getLowerCase( trim(adjustl(delim)) )
    end subroutine getDelim

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine closeFile( File )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: closeFile
#endif
        implicit none
        class(File_type), intent(inout) :: File
        character(*)    , parameter     :: PROCEDURE_NAME = "@close()"
        write(*,*) File%Path%modified
        inquire( file   = File%Path%modified    &
               , exist  = File%exists           &
               , opened = File%isOpen           &
               , number = File%number           &
               , iostat = File%Err%stat         &
               )
        if (File%Err%stat/=0) then
            File%Err%occurred = .true.
            File%Err%msg =  PROCEDURE_NAME // &
                            ": Error occurred while inquiring the open status and unit number of &
                            &file='" // File%Path%modified // "'."
            return
        end if
        if (File%exists) then
            if (File%isOpen) close(unit=File%number,iostat=File%Err%stat)
            File%Err = File%getCloseErr(File%Err%stat)
            if (File%Err%occurred) then
                File%Err%msg =    PROCEDURE_NAME // &
                                ": Error occurred while attempting to close the open file='" // File%Path%modified // "'."
            end if
        else
            ! check if the file with the original filename is open, and if so, close it.
            inquire( file   = File%Path%original    &
                   , exist  = File%exists           &
                   , opened = File%isOpen           &
                   , number = File%number           &
                   , iostat = File%Err%stat         &
                   )
            if (File%Err%stat/=0) then
                File%Err%occurred = .true.
                File%Err%msg =  PROCEDURE_NAME // &
                                ": Error occurred while inquiring the open status and unit number of &
                                &file='" // File%Path%original // "'."
                return
            end if
            if (File%exists) then
                if (File%isOpen) close(unit=File%number,iostat=File%Err%stat)
                File%Err = File%getCloseErr(File%Err%stat)
                if (File%Err%occurred) then
                    File%Err%msg =  PROCEDURE_NAME // &
                                    ": Error occurred while attempting to close the open file='" // File%Path%original // "'."
                end if
            end if
        end if
    end subroutine closeFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! sets values for File%unit, File%exists, File%isOpen, File%number, File%Err, and updates File%Path%modified (if needed)
    subroutine openFile( File )
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: openFile
#endif

        implicit none
        class(File_type), intent(inout) :: File
        character(*)    , parameter     :: PROCEDURE_NAME = MODULE_NAME // "@openFile()"

        ! if file is already open, first close it:
        inquire( file   = File%Path%original &
               , exist  = File%exists   &
               , opened = File%isOpen   &
               , number = File%number   &
               , iostat = File%Err%stat &
               )
        if (File%Err%stat/=0) then
            File%Err%occurred = .true.
            File%Err%msg =  PROCEDURE_NAME // &
                            ": Error occurred while inquiring the existence and open status, unit number of &
                            &file='" // File%Path%original // "'."
            return
        end if
        if (File%exists) then
            File%Path%modified = File%Path%original
            if (File%isOpen) then
                File%unit = File%number
                return
            else
                    write(*,*) File%status
                    write(*,*) File%Action%value
                    write(*,*) File%Access%value
                    write(*,*) File%Delim%value
                    write(*,*) File%Form%value
                    write(*,*) File%Position%value
                read(*,*)
                open( newunit  = File%unit              &
                    , file     = File%Path%modified     &
                    , form     = File%Form%value        &
                    , delim    = File%Delim%value       &
                    , status   = File%status            &
                    , action   = File%Action%value      &
                    , access   = File%Access%value      &
                    , iostat   = File%Err%stat          &
                    , position = File%Position%value    &
                    )
            end if
        else
            ! try the modified path file name
            inquire( file   = File%Path%modified    &
                   , exist  = File%exists           &
                   , opened = File%isOpen           &
                   , number = File%number           &
                   , iostat = File%Err%stat         &
                   )
            if (File%Err%stat/=0) then
                File%Err%occurred = .true.
                File%Err%msg =  PROCEDURE_NAME // &
                                ": Error occurred while inquiring the existence and open status, unit number of &
                                &file='" // File%Path%modified // "'."
                return
            end if
            if (File%exists) then
                if (File%isOpen) then
                    File%unit = File%number
                    return
                else
                    write(*,*) File%status
                    write(*,*) File%Action%value
                    write(*,*) File%Access%value
                    write(*,*) File%Delim%value
                    write(*,*) File%Form%value
                    write(*,*) File%Position%value
                    read(*,*)
                    open( newunit  = File%unit              &
                        , form     = File%Form%value        &
                        , delim    = File%Delim%value       &
                        , status   = File%status            &
                        , action   = File%Action%value      &
                        , access   = File%Access%value      &
                        , file     = File%Path%modified     &
                        , iostat   = File%Err%stat          &
                        , position = File%Position%value    &
                        )
                end if
            else
                File%Err%occurred = .true.
                File%Err%msg =  PROCEDURE_NAME // &
                                ": The requested file to open with possible addresses '" // File%Path%original // "' or '" // &
                                File%Path%modified // "' does not exist."
                return
            end if
        end if
    end subroutine openFile

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getWriteErr(stat) result(Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getWriteErr
#endif
        use Err_mod, only: Err_type
        implicit none
        integer, intent(in)     :: stat
        type(Err_type)          :: Err
        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME // "@getWriteErr()"
        Err%occurred = .false.
        Err%stat = stat
        Err%msg = ""
        if ( is_iostat_eor(Err%stat) ) then
            Err%occurred = .true.
            Err%msg  = PROCEDURE_NAME // ": End-Of-Record error condition occurred while attempting to write to file."
            return
        elseif ( is_iostat_end(Err%stat) ) then
            Err%occurred = .true.
            Err%msg  = PROCEDURE_NAME // ": End-Of-File error condition occurred while attempting to write to file."
            return
        elseif ( Err%stat>0 ) then
            Err%occurred = .true.
            Err%msg  = PROCEDURE_NAME // ": Unknown error condition occurred while attempting to write to file."
            return
        end if
    end function getWriteErr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getReadErr(stat,path) result(Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getReadErr
#endif
        use Err_mod, only: Err_type
        implicit none
        integer, intent(in)                 :: stat
        character(*), intent(in), optional  :: path
        type(Err_type)                      :: Err
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@getReadErr()"
        if (stat==0) then
            Err%occurred = .false.
            Err%stat = stat
            Err%msg = ""
            return
        else
            Err%occurred = .true.
            !write(*,*) stat
            Err%stat = stat
            if ( is_iostat_eor(stat) ) then
                Err%msg  = PROCEDURE_NAME // ": End-Of-Record error condition occurred while attempting to read from file."
            elseif ( is_iostat_end(stat) ) then
                Err%msg  = PROCEDURE_NAME // ": End-Of-File error condition occurred while attempting to read from file."
            elseif ( stat>0 ) then
                Err%msg  = PROCEDURE_NAME // ": Unknown error condition occurred while attempting to read from file."
            end if
            if (present(path)) Err%msg = Err%msg(1:len(Err%msg)-1) // "='" // path // "'."
        end if
    end function getReadErr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getCloseErr(stat) result(Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getCloseErr
#endif
        use Err_mod, only: Err_type
        implicit none
        integer, intent(in)     :: stat
        type(Err_type)          :: Err
        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME // "@getCloseErr()"
        Err%occurred = .false.
        Err%stat = stat
        Err%msg = ""
        if (Err%stat/=0) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred while attempting to close the open file."
            return
        end if
    end function getCloseErr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getOpenErr(stat) result(Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getOpenErr
#endif
        use Err_mod, only: Err_type
        implicit none
        integer, intent(in)     :: stat
        type(Err_type)          :: Err
        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME // "@getOpenErr()"
        Err%occurred = .false.
        Err%stat = stat
        Err%msg = ""
        if (Err%stat>0) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Unknown error occurred while opening file."
            return
        end if
    end function getOpenErr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function getInqErr(stat) result(Err)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: getInqErr
#endif
        use Err_mod, only: Err_type
        implicit none
        integer, intent(in)     :: stat
        type(Err_type)          :: Err
        character(*), parameter :: PROCEDURE_NAME = MODULE_NAME // "@getInqErr()"
        Err%occurred = .false.
        Err%stat = stat
        Err%msg = ""
        if (Err%stat/=0) then
            Err%occurred = .true.
            Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the status of file."
            return
        end if
    end function getInqErr

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function constructAction(value) result(Action)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructAction
#endif
        use String_mod, only: getLowerCase
        character(*), intent(in), optional  :: value
        type(Action_type)                   :: Action
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@constructAction()"
        if (present(value)) then
            Action%value = getLowerCase(trim(adjustl(value)))
            if (Action%value=="read") then
                Action%isRead = .true.
            elseif (Action%value=="write") then
                Action%isWrite = .true.
            elseif (Action%value=="readwrite") then
                Action%isReadWrite = .true.
            elseif (Action%value=="undefined") then
                Action%isUndefined = .true.
            else
                Action%value = ""
                Action%Err%occurred = .true.
                Action%Err%msg = PROCEDURE_NAME // ": Invalid requested Action%value='" // Action%value // "'."
            end if
        else
            Action%value = "readwrite"
            Action%isReadWrite = .true.
        end if
    end function constructAction

    function constructAccess(value) result(Access)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructAccess
#endif
        use String_mod, only: getLowerCase
        character(*), intent(in), optional  :: value
        type(Access_type)                   :: Access
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@constructAccess()"
        if (present(value)) then
            Access%value = getLowerCase(trim(adjustl(value)))
            if (Access%value=="sequential") then
                Access%isSequential = .true.
            elseif (Access%value=="direct") then
                Access%isDirect = .true.
            elseif (Access%value=="undefined") then
                Access%isUndefined = .true.
            else
                Access%value = ""
                Access%Err%occurred = .true.
                Access%Err%msg = PROCEDURE_NAME // ": Invalid requested Access%value='" // Access%value // "'."
            end if
        else
            Access%value = "sequential"
            Access%isSequential = .true.
        end if
    end function constructAccess

    function constructForm(value) result(Form)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructForm
#endif
        use String_mod, only: getLowerCase
        character(*), intent(in), optional  :: value
        type(Form_type)                     :: Form
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@constructForm()"
        if (present(value)) then
            Form%value = getLowerCase(trim(adjustl(value)))
            if (Form%value=="formatted") then
                Form%isFormatted = .true.
            elseif (Form%value=="unformatted") then
                Form%isUnformatted = .true.
            elseif (Form%value=="undefined") then
                Form%isUndefined = .true.
            else
                Form%value = ""
                Form%Err%occurred = .true.
                Form%Err%msg = PROCEDURE_NAME // ": Invalid requested Form%value='" // Form%value // "'."
            end if
        else
            Form%value = "formatted"
            Form%isFormatted = .true.
        end if
    end function constructForm

    function constructBlank(value) result(Blank)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructBlank
#endif
        use String_mod, only: getLowerCase
        character(*), intent(in), optional  :: value
        type(Blank_type)                    :: Blank
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@constructBlank()"
        if (present(value)) then
            Blank%value = getLowerCase(trim(adjustl(value)))
            if (Blank%value=="null") then
                Blank%isNull = .true.
            elseif (Blank%value=="zero") then
                Blank%isZero = .true.
            elseif (Blank%value=="undefined") then
                Blank%isUndefined = .true.
            else
                Blank%value = ""
                Blank%Err%occurred = .true.
                Blank%Err%msg = PROCEDURE_NAME // ": Invalid requested Blank%value='" // Blank%value // "'."
            end if
        else
            Blank%value = "null"
            Blank%isNull = .true.
        end if
    end function constructBlank

    function constructPosition(value) result(Position)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructPosition
#endif
        use String_mod, only: getLowerCase
        character(*), intent(in), optional  :: value
        type(Position_type)                 :: Position
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@constructPosition()"
        if (present(value)) then
            Position%value = getLowerCase(trim(adjustl(value)))
            if (Position%value=="asis") then
                Position%isAsis = .true.
            elseif (Position%value=="rewind") then
                Position%isRewind = .true.
            elseif (Position%value=="append") then
                Position%isAppend = .true.
            elseif (Position%value=="undefined") then
                Position%isUndefined = .true.
            else
                Position%value = ""
                Position%Err%occurred = .true.
                Position%Err%msg = PROCEDURE_NAME // ": Invalid requested Position%value='" // Position%value // "'."
            end if
        else
            Position%value = "asis"
            Position%isAsis = .true.
        end if
    end function constructPosition

    function constructDelim(value) result(Delim)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructDelim
#endif
        use String_mod, only: getLowerCase
        character(*), intent(in), optional  :: value
        type(Delim_type)                    :: Delim
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@constructDelim()"
        if (present(value)) then
            Delim%value = getLowerCase(trim(adjustl(value)))
            if (Delim%value=="quote") then
                Delim%isQuote = .true.
            elseif (Delim%value=="apostrophe") then
                Delim%isApostrophe = .true.
            elseif (Delim%value=="none") then
                Delim%isNone = .true.
            elseif (Delim%value=="undefined") then
                Delim%isUndefined = .true.
            else
                Delim%value = ""
                Delim%Err%occurred = .true.
                Delim%Err%msg = PROCEDURE_NAME // ": Invalid requested Delim%value='" // Delim%value // "'."
            end if
        else
            Delim%value = "none"
            Delim%isNone = .true.
        end if
    end function constructDelim

    function constructPad(value) result(Pad)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructPad
#endif
        use String_mod, only: getLowerCase
        character(*), intent(in), optional  :: value
        type(Pad_type)                      :: Pad
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@constructPad()"
        if (present(value)) then
            Pad%value = getLowerCase(trim(adjustl(value)))
            if (Pad%value=="yes") then
                Pad%isPadded = .true.
            elseif (Pad%value=="no") then
                Pad%isNotPadded = .true.
            elseif (Pad%value=="undefined") then
                Pad%isUndefined = .true.
            else
                Pad%value = ""
                Pad%Err%occurred = .true.
                Pad%Err%msg = PROCEDURE_NAME // ": Invalid requested Pad%value='" // Pad%value // "'."
            end if
        else
            Pad%value = "yes"
            Pad%isPadded = .true.
        end if
    end function constructPad

    function constructRound(value) result(Round)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructRound
#endif
        use String_mod, only: getLowerCase
        character(*), intent(in), optional  :: value
        type(Round_type)                    :: Round
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@constructRound()"
        if (present(value)) then
            Round%value = getLowerCase(trim(adjustl(value)))
            if (Round%value=="up") then
                Round%isUp = .true.
            elseif (Round%value=="down") then
                Round%isDown = .true.
            elseif (Round%value=="zero") then
                Round%isZero = .true.
            elseif (Round%value=="nearest") then
                Round%isNearest = .true.
            elseif (Round%value=="compatible") then
                Round%isCompatible = .true.
            elseif (Round%value=="processor_defined") then
                Round%isProcessDefined = .true.
            elseif (Round%value=="undefined") then
                Round%isUndefined = .true.
            else
                Round%value = ""
                Round%Err%occurred = .true.
                Round%Err%msg = PROCEDURE_NAME // ": Invalid requested Round%value='" // Round%value // "'."
            end if
        else
            Round%value = "processor_defined"
            Round%isProcessDefined = .true.
        end if
    end function constructRound

    function constructSign(value) result(Sign)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: constructSign
#endif
        use String_mod, only: getLowerCase
        character(*), intent(in), optional  :: value
        type(Sign_type)                     :: Sign
        character(*), parameter             :: PROCEDURE_NAME = MODULE_NAME // "@constructSign()"
        if (present(value)) then
            Sign%value = getLowerCase(trim(adjustl(value)))
            if (Sign%value=="suppress") then
                Sign%isSuppress = .true.
            elseif (Sign%value=="plus") then
                Sign%isPlus = .true.
            elseif (Sign%value=="processor_defined") then
                Sign%isProcessDefined = .true.
            elseif (Sign%value=="undefined") then
                Sign%isUndefined = .true.
            else
                Sign%value = ""
                Sign%Err%occurred = .true.
                Sign%Err%msg = PROCEDURE_NAME // ": Invalid requested Sign%value='" // Sign%value // "'."
            end if
        else
            Sign%value = "processor_defined"
            Sign%isProcessDefined = .true.
        end if
    end function constructSign

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module File_mod