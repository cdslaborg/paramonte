!***********************************************************************************************************************************
!***********************************************************************************************************************************
!
!   ParaMonte: plain powerful parallel Monte Carlo library.
!
!   Copyright (C) 2012-present, The Computational Data Science Lab
!
!   This file is part of the ParaMonte library.
!
!   ParaMonte is free software: you can redistribute it and/or modify it
!   under the terms of the GNU Lesser General Public License as published
!   by the Free Software Foundation, version 3 of the License.
!
!   ParaMonte is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with the ParaMonte library. If not, see,
!
!       https://github.com/cdslaborg/paramonte/blob/master/LICENSE
!
!   ACKNOWLEDGMENT
!
!   As per the ParaMonte library license agreement terms,
!   if you use any parts of this library for any purposes,
!   we ask you to acknowledge the use of the ParaMonte library
!   in your work (education/research/industry/development/...)
!   by citing the ParaMonte library as described on this page:
!
!       https://github.com/cdslaborg/paramonte/blob/master/ACKNOWLEDGMENT.md
!
!***********************************************************************************************************************************
!***********************************************************************************************************************************

module ParaMonte_mod

    use System_mod, only: SystemInfo_type
    use Decoration_mod, only: Decoration_type
    use Constants_mod, only: RK, IK, CIK, CRK, HUGE_IK, HUGE_RK
    use String_mod, only: IntStr_type
    use System_mod, only: OS_type
    use Timer_mod, only: Timer_type
    use File_mod, only: File_type
    use Path_mod, only: MAX_FILE_PATH_LEN
    use Err_mod, only: Err_type, informUser, note, warn, abort
    use SpecBase_mod, only: SpecBase_type

    implicit none

    public

    character(*), parameter         :: MODULE_NAME = "@ParaMonte_mod"

    type                            :: QuantileProbability_type
        integer(IK)                 :: count = 9_IK
        real(RK)                    :: Value(9) = [0._RK,0.05_RK,0.10_RK,0.25_RK,0.50_RK,0.75_RK,0.90_RK,0.95_RK,1.0_RK]
        character(4)                :: Name(9) = ["  Q0","  Q5"," Q10"," Q25"," Q50"," Q75"," Q90"," Q95","Q100"]
    end type QuantileProbability_type
    type(QuantileProbability_type), parameter :: QPROB = QuantileProbability_type()

    type                            :: Moment_type
        integer(IK)                 :: count = 0_IK
        real(RK), allocatable       :: Mean(:)
        real(RK), allocatable       :: CovMat(:,:)
        real(RK), allocatable       :: CorMat(:,:)
        real(RK), allocatable       :: Quantile(:,:)
    end type Moment_type

    type                            :: ParaMonteNumFunCall_type
        integer(IK)                 :: accepted                         ! accepted in the simulation
        integer(IK)                 :: acceptedRejected                 ! accepted + rejected function calls
    end type ParaMonteNumFunCall_type

    type                            :: ParaMonteLogFuncMode_type
        real(RK)                    :: val
        real(RK), allocatable       :: Crd(:)
    end type ParaMonteLogFuncMode_type

    type                            :: ParaMonteStatistics_type
        real(RK)                    :: avgTimePerFunCalInSec = 0._RK
        real(RK)                    :: avgCommTimePerFunCall = 0._RK
        type(Moment_type)           :: Sample
    end type ParaMonteStatistics_type

    !*******************************************************************************************************************************
    ! ParaMonte IO variables and types
    !*******************************************************************************************************************************

    type                            :: Image_type
        integer(IK)                 :: id, count
        logical                     :: isFirst, isNotFirst, isMaster, isNotMaster
        character(:), allocatable   :: name
    end type Image_type

    type, extends(File_type)        :: LogFile_type
        type(IntStr_type)           :: maxColWidth
        character(6)                :: suffix = "report"
    end type LogFile_type

    type, extends(File_type)        :: TimeFile_type
        character(8)                :: suffix = "progress"
    end type TimeFile_type

    type, extends(File_type)        :: ChainFile_type
        character(5)                :: suffix = "chain"
    end type ChainFile_type

    type, extends(File_type)        :: SampleFile_type
        character(6)                :: suffix = "sample"
    end type SampleFile_type

    type, extends(File_type)        :: RestartFile_type
        integer(IK)                 :: counter = 0_IK
        character(7)                :: suffix = "restart"
    end type RestartFile_type

    !*******************************************************************************************************************************
    ! ParaMonte type
    !*******************************************************************************************************************************

    type                                        :: ParaMonte_type
        type(IntStr_type)                       :: nd
        character(8)                            :: name
        character(16)                           :: brand
        character(:), allocatable               :: date
        character(:), allocatable               :: version
        logical                                 :: isDryRun
        logical                                 :: isFreshRun
        logical                                 :: procArgNeeded
        logical                                 :: procArgHasPriority
        logical                                 :: inputFileArgIsPresent
        type(OS_type)                           :: OS
        type(Err_type)                          :: Err
        type(Image_type)                        :: Image
        type(SpecBase_type)                     :: SpecBase
        !type(ParaMonteStatistics_type)         :: Stats
        type(SystemInfo_type)                   :: SystemInfo
        type(Timer_type)                        :: Timer
        type(File_type)                         :: InputFile
        type(LogFile_type)                      :: LogFile
        type(TimeFile_type)                     :: TimeFile
        type(ChainFile_type)                    :: ChainFile
        type(SampleFile_type)                   :: SampleFile
        type(RestartFile_type)                  :: RestartFile
        type(Decoration_type)                   :: Decor
    contains
        procedure, pass                         :: setupParaMonte
        procedure, pass                         :: addSplashScreen
        procedure, pass                         :: setupOutputFiles
        procedure, pass                         :: noteUserAboutEnvSetup
        procedure, pass                         :: addCompilerPlatformInfo
        procedure, pass                         :: warnUserAboutInputFilePresence
        procedure, pass                         :: setWarnAboutProcArgHasPriority
        procedure, nopass                       :: informUser, note, warn, abort
        procedure, nopass                       :: warnUserAboutMissingNamelist
    end type ParaMonte_type

!***********************************************************************************************************************************
!***********************************************************************************************************************************

contains

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    ! To be called by all images.
    ! Tasks: setup initial variables, as well as construct default and null values for SpecBase. Sets
    ! PM%InputFile%exists = .true. if the input file exists and opens and assigns to it a unit number and sets
    ! and PM%InputFile%isOpen = .true. if the opening process is successful.
    ! If the input file exists, the path used to open it successfully will be also written to InpuFile%Path%modified
    !subroutine setupParaMonte(PM,nd,name,date,version,inputFile)
    subroutine setupParaMonte(PM,nd,name,inputFile)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setupParaMonte
#endif
        use, intrinsic :: iso_fortran_env, only: output_unit
        use Decoration_mod, only: INDENT
        use Constants_mod, only: IK, NLC
        use String_mod, only: getLowerCase, num2str
        use System_mod, only: OS_type
        implicit none
        class(ParaMonte_type), intent(inout)    :: PM
        integer(IK), intent(in)                 :: nd
        character(*), intent(in)                :: name !, date, version
        character(*), intent(in), optional      :: inputFile
        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME // "@setupParaMonte()"

        PM%Timer = Timer_type(PM%Err)
        if (PM%Err%occurred) then
            PM%Err%msg = PROCEDURE_NAME // ": Error occurred while setting up the " // PM%name // "timer."//NLC// PM%Err%msg
            call PM%abort( Err = PM%Err, prefix = PM%brand, newline = NLC, outputUnit = PM%LogFile%unit )
            return
        end if
        PM%nd%val = nd
        PM%LogFile%unit = output_unit   ! temporarily set the report file to stdout.
        PM%Decor = Decoration_type()    ! initialize the TAB character and decoration symbol to the default values.

        PM%Err%occurred = .false.
        PM%Err%msg = ""

        PM%name     = name
        PM%brand    = INDENT // PM%name
#if defined IFORT_ENABLED || __GFORTRAN__
        PM%date     = "Build: " // __TIMESTAMP__
#else
        PM%date     = "Unknown Release Date"
#endif
#if defined PARAMONTE_VERSION
        PM%version  = "Version " // PARAMONTE_VERSION
#else
        PM%version  = "Unknown Version"
#endif

        ! setup general processor / coarray image variables

#if defined CAF_ENABLED
        PM%Image%id             = this_image()
        PM%Image%count          = num_images()
#elif defined MPI_ENABLED
        block
            use mpi
            integer(IK) :: ierrMPI
            logical     :: isInitialized
            call mpi_initialized( isInitialized, ierrMPI )
            if (.not. isInitialized) call mpi_init(ierrMPI)
            call mpi_comm_rank(mpi_comm_world, PM%Image%id, ierrMPI)
            call mpi_comm_size(mpi_comm_world, PM%Image%count, ierrMPI)
            PM%Image%id = PM%Image%id + 1_IK ! make the ranks consistent with Fortran coarray indexing conventions
        end block
#else
        PM%Image%id             = 1_IK
        PM%Image%count          = 1_IK
#endif

        PM%Image%name           = "@process(" // num2str(PM%Image%id) // ")"
        PM%Image%isFirst        = PM%Image%id==1_IK
        PM%Image%isNotFirst     = PM%Image%id/=1_IK
        PM%Image%isMaster       = .false.  ! ATTN: this will have to change later on, depending on the requested type of parallelism
        PM%Image%isNotMaster    = .false.

        ! setup formatting variables

        PM%nd%str = num2str(PM%nd%val)

        ! determine OS

        call PM%OS%query()
        if (PM%OS%Err%occurred) then
            PM%Err = PM%OS%Err
            PM%Err%msg = PROCEDURE_NAME // ": Error occurred while querying OS type."//NLC//PM%Err%msg
            call PM%abort( Err = PM%Err, prefix = PM%brand, newline = NLC, outputUnit = PM%LogFile%unit )
            return
        end if

        ! get system info by all images

        block
            use Constants_mod, only: IK
            use System_mod, only: SystemInfo_type
            integer(IK) :: irecord
            PM%SystemInfo = SystemInfo_type(OS=PM%OS)
            if (PM%SystemInfo%Err%occurred) then
                PM%Err = PM%SystemInfo%Err
                PM%Err%msg = PROCEDURE_NAME//": Error occurred while collecting system info."//NLC//PM%Err%msg
                call PM%abort( Err = PM%Err, prefix = PM%brand, newline = NLC, outputUnit = PM%LogFile%unit )
                return
            end if
        end block

        blockSplashByFirstImage: if (PM%Image%isFirst) then
            call PM%addSplashScreen()
            call PM%noteUserAboutEnvSetup()
        end if blockSplashByFirstImage

        ! check if input file exists by all images

        PM%InputFile%isInternal = .false.
        PM%inputFileArgIsPresent = present(inputFile)
        blockInputFileExistence: if (PM%inputFileArgIsPresent) then
            PM%InputFile = File_type( path=inputFile, status="old", OS=PM%OS )
            if (PM%InputFile%Err%occurred) then
                PM%Err = PM%InputFile%Err
                PM%Err%msg = PROCEDURE_NAME//": Error occurred while attempting to setup the user's input file='"//inputFile//"'."//NLC//PM%Err%msg
                call PM%abort( Err = PM%Err, prefix = PM%brand, newline = NLC, outputUnit = PM%LogFile%unit )
                return
            end if
            ! determine if the file is internal
            PM%InputFile%isInternal = PM%inputFileArgIsPresent .and. .not.PM%InputFile%exists .and. index(getLowerCase(inputFile),"&"//getLowerCase(PM%name)) > 0
            if (.not.(PM%InputFile%isInternal .or. PM%InputFile%exists)) then
                ! file is given, but is neither a path to an external file, nor an internal file containing a namelist
                ! Therefore, there must be an error/mistake by the user.
                PM%Err%msg =    "The user's input file='"//inputFile//"' is neither the path to an existing " // &
                                "external input file nor a string containing the input "//PM%name//" specifications namelist. &
                                &This may be due to the user's mistake, providing a wrong path to the input external file or &
                                &a wrong list of input specifications for the " //PM%name// " simulation. " //PM%name//" will &
                                &assume no input specifications file is given by the user..."
                if (PM%Image%isFirst) call PM%warn(msg=PM%Err%msg, prefix=PM%brand, newline=NLC, outputUnit=PM%LogFile%unit)
                PM%inputFileArgIsPresent = .false.
            end if
        end if blockInputFileExistence

        if (PM%Image%isFirst) call PM%warnUserAboutInputFilePresence()

        ! Set the default and null values for ParaMonte SpecBase

        PM%SpecBase = SpecBase_type(nd=PM%nd%val,methodName=PM%name,imageID=PM%Image%id,imageCount=PM%Image%count)

    end subroutine setupParaMonte

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine addSplashScreen(PM)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: addSplashScreen
#endif
        implicit none
        class(ParaMonte_type), intent(inout) :: PM

        PM%Decor%text = &
        "\n\n"// &
        !PM%name // "\n" // &
        "ParaMonte\n"// &
        "Plain Powerful Parallel\n"// &
        "Monte Carlo Library\n"// &
        "\n"// &
        PM%version // "\n" // &
        "\n"// &
        PM%date // "\n" // &
        "\n"// &
        "Department of Physics\n"// &
        "Computational & Data Science Lab\n"// &
        "Data Science Program, College of Science\n"// &
        "The University of Texas at Arlington\n"// &
        "\n"// &
        "originally developed at\n"// &
        "\n"// &
        "Multiscale Modeling Group\n"// &
        "Center for Computational Oncology (CCO)\n"// &
        "Oden Institute for Computational Engineering and Sciences\n"// &
        "Department of Aerospace Engineering and Engineering Mechanics\n"// &
        "Department of Neurology, Dell-Seton Medical School\n"// &
        "Department of Biomedical Engineering\n"// &
        "The University of Texas at Austin\n"// &
        "\n"// &
        "For questions and further information, please contact:\n"// &
        "\n"// &
        "Amir Shahmoradi\n"// &
        "\n"// &
        "shahmoradi@utexas.edu\n"// &
        "amir.shahmoradi@uta.edu\n"// &
        "ashahmoradi@gmail.com\n"// &
        !"amir@physics.utexas.edu\n"// &
        !"amir@austin.utexas.edu\n"// &
        !"amir@ph.utexas.edu\n"// &
        "\n"// &
        !"https://www.shahmoradi.org\n"// &
        !"shahmoradi.org\n"// &
        "cdslab.org/pm\n"// &
        "\n"// &
        "https://www.cdslab.org/paramonte/\n"// &
        "\n"

        call PM%Decor%writeDecoratedText( text=PM%Decor%text &
                                        , symbol="*" &
                                        , width=132 &
                                        , thicknessHorz=4 &
                                        , thicknessVert=2 &
                                        , marginTop=1 &
                                        , marginBot=2 &
                                        , outputUnit=PM%LogFile%unit &
                                        , newLine="\n" &
                                        )

    end subroutine addSplashScreen

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine addCompilerPlatformInfo(PM)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: addCompilerPlatformInfo
#endif
        use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options
        implicit none
        class(ParaMonte_type), intent(inout)    :: PM
        integer(IK)                             :: i, j

        ! report the interface type to ParaMonte

!#if defined CFI_ENABLED
        call PM%Decor%writeDecoratedText( text="\nParaMonte library interface specifications\n" &
                                        , symbol="*" &
                                        , width=132 &
                                        , thicknessHorz=4 &
                                        , thicknessVert=1 &
                                        , marginTop=2 &
                                        , marginBot=1 &
                                        , outputUnit=PM%LogFile%unit &
                                        , newLine="\n" &
                                        )
        PM%Decor%List = PM%Decor%wrapText( PM%SpecBase%InterfaceType%val , 132 )
        do i = 1,size(PM%Decor%List)
            write(PM%LogFile%unit,"(*(g0))") PM%Decor%List(i)%record
        end do
!#endif

        ! report the ParaMonte compiler version and options

        call PM%Decor%writeDecoratedText( text="\nParaMonte library compiler version\n" &
                                        , symbol="*" &
                                        , width=132 &
                                        , thicknessHorz=4 &
                                        , thicknessVert=1 &
                                        , marginTop=2 &
                                        , marginBot=1 &
                                        , outputUnit=PM%LogFile%unit &
                                        , newLine="\n" &
                                        )
        PM%Decor%List = PM%Decor%wrapText( compiler_version() , 132 )
        do i = 1,size(PM%Decor%List)
            write(PM%LogFile%unit,"(*(g0))") PM%Decor%List(i)%record
        end do

        call PM%Decor%writeDecoratedText( text="\nParaMonte library compiler options\n" &
                                        , symbol="*" &
                                        , width=132 &
                                        , thicknessHorz=4 &
                                        , thicknessVert=1 &
                                        , marginTop=2 &
                                        , marginBot=1 &
                                        , outputUnit=PM%LogFile%unit &
                                        , newLine="\n" &
                                        )
        PM%Decor%List = PM%Decor%wrapText( compiler_options() , 132 )
        do i = 1,size(PM%Decor%List)
            write(PM%LogFile%unit,"(*(g0))") PM%Decor%List(i)%record
        end do

        call PM%Decor%writeDecoratedText( text="\nRuntime platform specifications\n" &
                                        , symbol="*" &
                                        , width=132 &
                                        , thicknessHorz=4 &
                                        , thicknessVert=1 &
                                        , marginTop=2 &
                                        , marginBot=1 &
                                        , outputUnit=PM%LogFile%unit &
                                        , newLine="\n" &
                                        )
        do j = 1, PM%SystemInfo%nRecord
            PM%Decor%List = PM%Decor%wrapText( PM%SystemInfo%List(j)%record , 132 )
            do i = 1,size(PM%Decor%List)
                write(PM%LogFile%unit,"(*(g0))") PM%Decor%List(i)%record
            end do
        end do
        call PM%Decor%write(PM%LogFile%unit)

    end subroutine addCompilerPlatformInfo

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine noteUserAboutEnvSetup(PM)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: noteUserAboutEnvSetup
#endif
        implicit none
        class(ParaMonte_type), intent(inout) :: PM
        call PM%Decor%writeDecoratedText( text = "\nSetting up the " // PM%name // " simulation environment\n" &
                                        , marginTop = 1     &
                                        , marginBot = 1     &
                                        , newline = "\n"    &
                                        , outputUnit = PM%LogFile%unit )
    end subroutine noteUserAboutEnvSetup

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine warnUserAboutMissingNamelist(prefix,name,namelist,outputUnit)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: warnUserAboutMissingNamelist
#endif
        use, intrinsic :: iso_fortran_env, only: output_unit
        use Constants_mod, only: IK
        use Err_mod, only: warn
        implicit none
        character(*), intent(in)    :: prefix, name, namelist
        integer(IK) , intent(in)    :: outputUnit
        character(:), allocatable   :: msg
        msg = "No namelist group of variables named "//namelist//" was detected in user's input file for "//name//" options.\n"//&
              "All " // name // " options will be assigned appropriate default values."
        call warn( prefix = prefix, outputUnit = outputUnit, newline = "\n", msg = msg )
        if (outputUnit/=output_unit) then
            call warn( prefix = prefix, outputUnit = output_unit, newline = "\n", msg = msg )
        end if
    end subroutine warnUserAboutMissingNamelist

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine warnUserAboutInputFilePresence(PM)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: warnUserAboutInputFilePresence
#endif
        use Constants_mod, only: NLC
        implicit none
        class(ParaMonte_type), intent(inout) :: PM
#if defined CFI_ENABLED
        if (PM%SpecBase%InterfaceType%isPython) then
            call PM%note( prefix     = PM%brand &
                        , outputUnit = PM%LogFile%unit &
                        , newline    = NLC &
                        , msg        = "Interfacing Python with "// PM%name //"..." )
        elseif (PM%SpecBase%InterfaceType%isClang) then
#else
            if (PM%inputFileArgIsPresent) then
                if (PM%InputFile%exists) then
                    call PM%note( prefix     = PM%brand &
                                , outputUnit = PM%LogFile%unit &
                                , newline    = NLC &
                                , msg        = "The user's input file for " // PM%name // " options was detected."// NLC // &
                                               "All " // PM%name // " options will be read from the input file."// NLC // &
                                               "Here is " // PM%name // " input options file:"// NLC // PM%inputFile%Path%modified )
                elseif (PM%InputFile%isInternal) then
                    call PM%note( prefix     = PM%brand &
                                , outputUnit = PM%LogFile%unit &
                                , newline    = NLC &
                                , msg        = "No external file corresponding to the user's input file for "//PM%name//" options &
                                               &could be found."//NLC//"The user-provided input file will be processed as an input &
                                               &string of "//PM%name//" options." )
                end if
            else
                call PM%note( prefix     = PM%brand &
                            , outputUnit = PM%LogFile%unit &
                            , newline    = NLC &
                            , msg        = "No " // PM%name  // " input file is provided by the user."//NLC//&
                                           "Variable values from the procedure arguments will be used instead, where provided."//NLC//&
                                           "Otherwise, the default options will be used." )
            end if
#endif
#if defined CFI_ENABLED
        end if
#endif
    end subroutine warnUserAboutInputFilePresence

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setWarnAboutProcArgHasPriority(PM)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setWarnAboutProcArgHasPriority
#endif
        implicit none
        class(ParaMonte_type), intent(inout) :: PM
#if !defined CFI_ENABLED
        character(:), allocatable            :: msg
#endif
        PM%procArgHasPriority = .not. PM%SpecBase%InputFileHasPriority%val
        PM%procArgNeeded = PM%procArgHasPriority .or. (.not.PM%inputFileArgIsPresent)
#if !defined CFI_ENABLED
        if (PM%procArgHasPriority) then
            msg =   "Variable inputFileHasPriority = .false.\n&
                    &All variable values will be overwritten by the corresponding procedure argument values,\n&
                    &only if provided as procedure arguments."
        else
            msg =   "Variable inputFileHasPriority = .true.\n&
                    &All variable values will be read from the user-provided input file"
        end if
        call PM%note( prefix = PM%brand, outputUnit = PM%LogFile%unit, newline = "\n", msg = msg )
#endif
    end subroutine setWarnAboutProcArgHasPriority

!***********************************************************************************************************************************
!***********************************************************************************************************************************

    subroutine setupOutputFiles(PM)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setupOutputFiles
#endif
        use Decoration_mod, only: getGenericFormat, INDENT
        use Constants_mod, only: NLC, FILE_EXT
        use String_mod, only: num2str
        use Path_mod, only: MAX_FILE_PATH_LEN, mkdir

        implicit none

        class(ParaMonte_type), intent(inout)    :: PM
        character(:), allocatable               :: msg, workingOn, currentWorkingDir
        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME // "@setupOutputFiles()"

        if (PM%SpecBase%OutputFileName%original==PM%SpecBase%OutputFileName%def) then
            msg =   "No user-input filename prefix for " // PM%name // " output files detected." // NLC // &
                    "Generating appropriate filenames for " // PM%name // " output files from the current date and time..."
        else
            msg = "Variable outputFileName detected among the input variables to " // PM%name // ":" //NLC// PM%SpecBase%OutputFileName%original
        end if

        call PM%SpecBase%OutputFileName%query(OS=PM%OS)

        if (PM%SpecBase%OutputFileName%Err%occurred) then
            PM%Err = PM%SpecBase%OutputFileName%Err
            PM%Err%msg = PROCEDURE_NAME // ": Error occurred while attempting to construct OutputFileName path type." //NLC// PM%Err%msg
            call PM%abort( Err = PM%Err, prefix = PM%brand, newline = NLC, outputUnit = PM%LogFile%unit )
            return
        end if

        PM%SpecBase%OutputFileName%namePrefix = PM%SpecBase%OutputFileName%name // PM%SpecBase%OutputFileName%ext

        ! get the current working directory

        if (allocated(currentWorkingDir)) deallocate(currentWorkingDir)
        allocate( character(MAX_FILE_PATH_LEN) :: currentWorkingDir )
        block
#if defined IFORT_ENABLED
            use ifport ! only: getcwd
#endif
#if defined IFORT_ENABLED || __GFORTRAN__
            PM%Err%stat = getcwd(currentWorkingDir)
            currentWorkingDir = trim(adjustl(currentWorkingDir))
#else
            PM%Err%stat = 0_IK
            currentWorkingDir = "."
#endif
        end block
        if (PM%Err%stat/=0) then
            PM%Err%msg = PROCEDURE_NAME//": Error occurred while fetching the current working directory via getcwd()."//NLC
            call PM%abort( Err = PM%Err, prefix = PM%brand, newline = NLC, outputUnit = PM%LogFile%unit )
            return
        end if
        msg = msg //NLC//NLC// "Absolute path to the current working directory:"//NLC//currentWorkingDir

        if (len_trim(adjustl(PM%SpecBase%OutputFileName%dir))==0) then
            PM%SpecBase%OutputFileName%dir = trim(adjustl(currentWorkingDir)) // PM%SpecBase%OutputFileName%slashOS
            msg = msg //NLC//NLC// "All output files will be written to the current working directory:"//NLC//PM%SpecBase%OutputFileName%dir
        else
            msg = msg //NLC//NLC// "Generating the requested directory for ParaDRAM output files:"//NLC//PM%SpecBase%OutputFileName%dir
        end if

        ! Generate the output files directory:

        if (PM%Image%isFirst) then
            PM%Err = mkdir( dirPath = PM%SpecBase%OutputFileName%dir, isWindows = PM%OS%isWindows )
            if (PM%Err%occurred) then
                PM%Err%msg = PROCEDURE_NAME//": Error occurred while making directory = '"//PM%SpecBase%OutputFileName%dir//"'."//NLC//PM%Err%msg
                call PM%abort( Err = PM%Err, prefix = PM%brand, newline = NLC, outputUnit = PM%LogFile%unit )
                return
            end if
        end if

        ! in parallel mode, ensure the directory exists before moving on

#if defined CAF_ENABLED
        sync all
#elif defined MPI_ENABLED
        block
            use mpi
            integer :: ierrMPI
            call mpi_barrier(mpi_comm_world,ierrMPI)
        end block
#endif

        if (len_trim(adjustl(PM%SpecBase%OutputFileName%namePrefix))==0) then
            msg = msg //NLC//NLC// "No user-input filename prefix for " // PM%name // " output files detected."//NLC//&
            "Generating appropriate filenames for " // PM%name // " output files from the current date and time..."
            PM%SpecBase%OutputFileName%namePrefix = PM%SpecBase%OutputFileName%def
        end if

        PM%SpecBase%OutputFileName%pathPrefix = PM%SpecBase%OutputFileName%dir // PM%SpecBase%OutputFileName%namePrefix

        ! Variable msg will be used down this subroutine, so it should not be changed beyond this point
        msg  =  msg //NLC//NLC// PM%name // " output files will be prefixed with:"//NLC// PM%SpecBase%OutputFileName%pathPrefix
        if (PM%Image%isFirst) then
            call PM%note( prefix = PM%brand, outputUnit = PM%LogFile%unit, newline = NLC, msg = msg )
        end if

        ! Generate the output filenames, search for pre-existing runs, and open the report file:

        if (PM%Image%isFirst) call PM%note( prefix = PM%brand, outputUnit = PM%LogFile%unit, newline = NLC, msg = "Searching for previous runs of " // PM%name // "..." )

        ! this block could be all executed by only the master images

        PM%LogFile%Path%Ext = FILE_EXT%ascii
        PM%TimeFile%Path%Ext = FILE_EXT%ascii
        PM%ChainFile%Path%Ext = FILE_EXT%ascii
        PM%SampleFile%Path%Ext = FILE_EXT%ascii
        PM%RestartFile%Path%Ext = FILE_EXT%ascii

        PM%RestartFile%Form%value = "formatted"
        if (PM%SpecBase%RestartFileFormat%isBinary) then
            PM%RestartFile%Form%value = "unformatted"
            PM%RestartFile%Path%Ext = FILE_EXT%binary
        end if

        PM%ChainFile%Form%value = "formatted"
        if (PM%SpecBase%ChainFileFormat%isBinary) then
            PM%ChainFile%Form%value = "unformatted"
            PM%ChainFile%Path%Ext = FILE_EXT%binary
        end if

        block
            use Path_mod, only: Path_type
            !use String_mod, only: num2str
            character(:), allocatable :: fullOutputFileName
            integer(IK) :: imageID

#if defined CAF_ENABLED || defined MPI_ENABLED
            if (PM%SpecBase%ParallelizationModel%isMultiChain) then
                imageID = PM%Image%id
            else
#endif
                imageID = 1_IK
#if defined CAF_ENABLED || defined MPI_ENABLED
            end if
#endif
            fullOutputFileName  = PM%SpecBase%OutputFileName%pathPrefix // "_process_" // num2str(imageID) // "_"
            PM%LogFile%Path     = Path_type( inputPath = fullOutputFileName // PM%LogFile%suffix        // PM%LogFile%Path%Ext      , OS = PM%OS )
            PM%TimeFile%Path    = Path_type( inputPath = fullOutputFileName // PM%TimeFile%suffix       // PM%TimeFile%Path%Ext     , OS = PM%OS )
            PM%ChainFile%Path   = Path_type( inputPath = fullOutputFileName // PM%ChainFile%suffix      // PM%ChainFile%Path%Ext    , OS = PM%OS )
            PM%SampleFile%Path  = Path_type( inputPath = fullOutputFileName // PM%SampleFile%suffix     // PM%SampleFile%Path%Ext   , OS = PM%OS )
            PM%RestartFile%Path = Path_type( inputPath = fullOutputFileName // PM%RestartFile%suffix    // PM%RestartFile%Path%Ext  , OS = PM%OS )
        end block

        inquire( file = PM%LogFile%Path%original, exist = PM%LogFile%exists, iostat = PM%LogFile%Err%stat )
        PM%Err = PM%LogFile%getInqErr( PM%LogFile%Err%stat )
        if (PM%Err%occurred) then
            PM%Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file='" // PM%LogFile%Path%original // PM%Err%msg
            call PM%abort( Err = PM%Err, prefix = PM%brand, newline = NLC, outputUnit = PM%LogFile%unit )
            return
        end if

        inquire( file = PM%SampleFile%Path%original, exist = PM%SampleFile%exists, iostat = PM%SampleFile%Err%stat )
        PM%Err = PM%SampleFile%getInqErr( PM%SampleFile%Err%stat )
        if (PM%Err%occurred) then
            PM%Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file='" // PM%SampleFile%Path%original // PM%Err%msg
            call PM%abort( Err = PM%Err, prefix = PM%brand, newline = NLC, outputUnit = PM%LogFile%unit )
            return
        end if

        inquire( file = PM%TimeFile%Path%original, exist = PM%TimeFile%exists, iostat = PM%TimeFile%Err%stat )
        PM%Err = PM%TimeFile%getInqErr( PM%TimeFile%Err%stat )
        if (PM%Err%occurred) then
            PM%Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file='" // PM%TimeFile%Path%original // PM%Err%msg
            call PM%abort( Err = PM%Err, prefix = PM%brand, newline = NLC, outputUnit = PM%LogFile%unit )
            return
        end if

        inquire( file = PM%ChainFile%Path%original, exist = PM%ChainFile%exists, iostat = PM%ChainFile%Err%stat )
        PM%Err = PM%ChainFile%getInqErr( PM%ChainFile%Err%stat )
        if (PM%Err%occurred) then
            PM%Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file='" // PM%ChainFile%Path%original // PM%Err%msg
            call PM%abort( Err = PM%Err, prefix = PM%brand, newline = NLC, outputUnit = PM%LogFile%unit )
            return
        end if

        inquire( file = PM%RestartFile%Path%original, exist = PM%RestartFile%exists, iostat = PM%RestartFile%Err%stat )
        PM%Err = PM%RestartFile%getInqErr( PM%RestartFile%Err%stat )
        if (PM%Err%occurred) then
            PM%Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file='" // PM%RestartFile%Path%original // PM%Err%msg
            call PM%abort( Err = PM%Err, prefix = PM%brand, newline = NLC, outputUnit = PM%LogFile%unit )
            return
        end if

        PM%isDryRun = PM%LogFile%exists .or. PM%TimeFile%exists .or. PM%RestartFile%exists .or. PM%ChainFile%exists .or. PM%SampleFile%exists ! not fresh, if any file exists
        PM%isFreshRun = .not. PM%isDryRun

        if (PM%isFreshRun) then
            if (PM%Image%isFirst) call PM%note( prefix = PM%brand, outputUnit = PM%LogFile%unit, newline = NLC, msg = "No pre-existing "//PM%name//" run detected."//NLC//"Starting a fresh ParaDRAM run..." )
        else
            if (PM%Image%isFirst) call PM%note( prefix = PM%brand, outputUnit = PM%LogFile%unit, newline = NLC, msg = "Previous run of "//PM%name//" detected."//NLC//"Searching for restart files..." )
            if (PM%SampleFile%exists) then ! sampling is already complete
                PM%Err%occurred = .true.
                PM%Err%msg =    PROCEDURE_NAME//": Error occurred. Output sample file detected: "//PM%SampleFile%Path%original//&
                                NLC//PM%name//" cannot overwrite an already-completed simulation."//&
                                NLC//"Please provide an alternative file name for the new simulation outputs."
            elseif (PM%LogFile%exists .and. PM%TimeFile%exists .and. PM%RestartFile%exists .and. PM%ChainFile%exists) then  ! restart mode
                if (PM%SpecBase%SampleSize%val==0_IK) then ! sampling is already complete
                    PM%Err%occurred = .true.
                    PM%Err%msg = PROCEDURE_NAME//": Error occurred. The input variable sampleSize = 0 indicates that the output files belong to a completed simulation."
                else
                    PM%Err%occurred = .false.
                end if
            else
                PM%Err%occurred = .true.
                PM%Err%msg = PROCEDURE_NAME//": Error occurred. For a successful simulation restart, all output files are necessary."//NLC//"List of missing simulation output files:"
                if (.not. PM%LogFile%exists)        PM%Err%msg = PM%Err%msg//NLC//PM%LogFile%Path%original
                if (.not. PM%TimeFile%exists)       PM%Err%msg = PM%Err%msg//NLC//PM%TimeFile%Path%original
                if (.not. PM%ChainFile%exists)      PM%Err%msg = PM%Err%msg//NLC//PM%ChainFile%Path%original
                if (.not. PM%RestartFile%exists)    PM%Err%msg = PM%Err%msg//NLC//PM%RestartFile%Path%original
            end if
            if (PM%Err%occurred) then
                call PM%abort( Err = PM%Err, prefix = PM%brand, newline = NLC, outputUnit = PM%LogFile%unit )
                return
            end if
        end if

        ! open/append the output files:

        if (PM%Image%isMaster) then
            if (PM%isFreshRun) then
                workingOn = "Generating the output "
                PM%LogFile%status = "new"
                PM%TimeFile%status = "new"
                PM%ChainFile%status = "new"
                PM%SampleFile%status = "new"
                PM%RestartFile%status = "new"
                PM%LogFile%Position%value = "asis"
                PM%TimeFile%Position%value = "asis"
                PM%ChainFile%Position%value = "asis"
                PM%SampleFile%Position%value = "asis"
                PM%RestartFile%Position%value = "asis"
            else
                workingOn = "Appending to the existing "
                PM%LogFile%status = "old"
                PM%TimeFile%status = "old"
                PM%ChainFile%status = "old"
                PM%SampleFile%status = "replace"
                PM%RestartFile%status = "old"
                PM%LogFile%Position%value = "append"
                PM%TimeFile%Position%value = "asis"
                PM%ChainFile%Position%value = "asis"
                PM%SampleFile%Position%value = "asis"
                PM%RestartFile%Position%value = "asis"
            end if
        end if

        ! print the stdout message for generating / appending the output report file

        blockLogFileListByFirstImage: if (PM%Image%isFirst) then

            ! print the stdout message for generating / appending the output report file(s)

            call PM%note( prefix = PM%brand             &
                        , outputUnit = PM%LogFile%unit  &
                        , newline = NLC                 &
                        , marginBot = 0_IK              &
                        , msg = workingOn // PM%LogFile%suffix // " file:" )

            ! print the the output report file name of the images

            call PM%note( prefix = PM%brand             &
                        , outputUnit = PM%LogFile%unit  &
                        , newline = NLC                 &
                        , marginTop = 0_IK              &
                        , marginBot = 0_IK              &
                        , msg = PM%LogFile%Path%original )

#if defined CAF_ENABLED || defined MPI_ENABLED
            if (PM%SpecBase%ParallelizationModel%isMultiChain) then
                block
                    use String_mod, only: replaceStr !, num2str
                    integer(IK) :: imageID
                    do imageID = 2, PM%Image%count
                        call PM%note( prefix = PM%brand             &
                                    , outputUnit = PM%LogFile%unit  &
                                    , newline = NLC                 &
                                    , marginTop = 0_IK              &
                                    , marginBot = 0_IK              &
                                    , msg = replaceStr( string = PM%LogFile%Path%original, search = "process_1", substitute = "process_"//num2str(imageID) ) )
                    end do
                end block
            end if
            PM%Err%msg = "Running the simulation in parallel on " // num2str(PM%Image%count) // " processes." // NLC
#else
            PM%Err%msg = "Running the simulation in serial on " // num2str(PM%Image%count) // " process." // NLC
#endif

            call PM%note( prefix = PM%brand             &
                        , outputUnit = PM%LogFile%unit  &
                        , newline = NLC                 &
                        , marginTop = 3_IK              &
                        , marginBot = 3_IK              &
                        , msg = PM%Err%msg // "Please see the output " // PM%LogFile%suffix // " and " // PM%TimeFile%suffix // " files for further realtime simulation details." &
                        )

        end if blockLogFileListByFirstImage

        ! ensure all images sync here to avoid wrong inquire result for the existence of the files

#if defined CAF_ENABLED
        sync all
#elif defined MPI_ENABLED
        block
            use mpi
            integer :: ierrMPI
            call mpi_barrier(mpi_comm_world,ierrMPI)
        end block
#endif

        ! open the output files
        ! Intel ifort SHARED attribute is essential for file unlocking

        blockMasterFileSetup: if (PM%Image%isMaster) then

            PM%LogFile%unit = 1001 + PM%Image%id ! for some unknown reason, if newunit is used, GFortran opens the file as an internal file
            open( unit = PM%LogFile%unit                    &
                , file = PM%LogFile%Path%original           &
                , status = PM%LogFile%status                &
                , iostat = PM%LogFile%Err%stat              &
#if defined IFORT_ENABLED && defined OS_IS_WINDOWS
                , SHARED                                    &
#endif
                , position = PM%LogFile%Position%value      )
            PM%Err = PM%LogFile%getOpenErr(PM%LogFile%Err%stat)
            if (PM%Err%occurred) then
                PM%Err%msg = PROCEDURE_NAME // ": Error occurred while opening the " // PM%name // " " // PM%LogFile%suffix // " file='" // PM%LogFile%Path%original // "'. "
                if (scan(" ",trim(adjustl(PM%LogFile%Path%original)))/=0) then
                    PM%Err%msg = PM%Err%msg // "It appears that absolute path used for the output files contains whitespace characters. " &
                                            // "This could be one potential cause of the simulation failure. " &
                                            // "The whitespace characters are always problematic in paths. " &
                                            // "Ensure the path used for the output files does not contain whitespace characters. "
                end if
                call PM%abort( Err = PM%Err, prefix = PM%brand, newline = NLC, outputUnit = PM%LogFile%unit )
                return
            end if

            ! rewrite the same old stuff to all report files

            if (PM%isFreshRun) then
                call PM%addSplashScreen()
                if (PM%SpecBase%SilentModeRequested%isFalse) call PM%addCompilerPlatformInfo()   ! this takes about 0.75 seconds to execute on Stampede Login nodes.
                call PM%noteUserAboutEnvSetup()
                call PM%warnUserAboutInputFilePresence()
                call PM%setWarnAboutProcArgHasPriority()
                call PM%note( prefix = PM%brand, outputUnit = PM%LogFile%unit, newline = NLC, msg = msg )
            end if

            ! open/append the output files

            if (PM%isFreshRun) call PM%note( prefix = PM%brand, outputUnit = PM%LogFile%unit, newline = NLC, msg = workingOn//PM%TimeFile%suffix//" file:"//NLC//PM%TimeFile%Path%original )

            PM%TimeFile%unit = 1000001  ! for some unknown reason, if newunit is used, GFortran opens the file as an internal file
            open( unit = PM%TimeFile%unit                   &
                , file = PM%TimeFile%Path%original          &
                , status = PM%TimeFile%status               &
                , iostat = PM%TimeFile%Err%stat             &
#if defined IFORT_ENABLED && defined OS_IS_WINDOWS
                , SHARED                                    &
#endif
                , position = PM%TimeFile%Position%value     )
            PM%Err = PM%TimeFile%getOpenErr(PM%TimeFile%Err%stat)
            if (PM%Err%occurred) then
                PM%Err%msg = PROCEDURE_NAME // ": Error occurred while opening the " // PM%name // " " // PM%TimeFile%suffix // " file='" // PM%TimeFile%Path%original // "'. "
                call PM%abort( Err = PM%Err, prefix = PM%brand, newline = NLC, outputUnit = PM%LogFile%unit )
                return
            end if

            if (PM%isFreshRun) call PM%note( prefix = PM%brand, outputUnit = PM%LogFile%unit, newline = NLC, msg = workingOn//PM%ChainFile%suffix//"file:"//NLC//PM%ChainFile%Path%original )

            PM%ChainFile%unit = 2000001  ! for some unknown reason, if newunit is used, GFortran opens the file as an internal file
            open( unit = PM%ChainFile%unit                  &
                , file = PM%ChainFile%Path%original         &
                , form = PM%ChainFile%Form%value            &
                , status = PM%ChainFile%status              &
                , iostat = PM%ChainFile%Err%stat            &
#if defined IFORT_ENABLED && defined OS_IS_WINDOWS
                , SHARED                                    &
#endif
                , position = PM%ChainFile%Position%value    )
            PM%Err = PM%ChainFile%getOpenErr(PM%ChainFile%Err%stat)
            if (PM%Err%occurred) then
                PM%Err%msg = PROCEDURE_NAME // ": Error occurred while opening the " // PM%name // " " // PM%ChainFile%suffix // " file='" // PM%ChainFile%Path%original // "'. "
                call PM%abort( Err = PM%Err, prefix = PM%brand, newline = NLC, outputUnit = PM%LogFile%unit )
                return
            end if

            PM%RestartFile%unit = 3000001  ! for some unknown reason, if newunit is used, GFortran opens the file as an internal file
            open( unit = PM%RestartFile%unit                &
                , file = PM%RestartFile%Path%original       &
                , form = PM%RestartFile%Form%value          &
                , status = PM%RestartFile%status            &
                , iostat = PM%RestartFile%Err%stat          &
#if defined IFORT_ENABLED && defined OS_IS_WINDOWS
                , SHARED                                    &
#endif
                , position = PM%RestartFile%Position%value  )
            PM%Err = PM%RestartFile%getOpenErr(PM%RestartFile%Err%stat)
            if (PM%Err%occurred) then
                PM%Err%msg = PROCEDURE_NAME // ": Error occurred while opening the " // PM%name // " " // PM%RestartFile%suffix // " file='" // PM%RestartFile%Path%original // "'. "
                call PM%abort( Err = PM%Err, prefix = PM%brand, newline = NLC, outputUnit = PM%LogFile%unit )
                return
            end if

            if (PM%isFreshRun) then
                call PM%Decor%writeDecoratedText( text = NLC // PM%name // " simulation specifications" // NLC &
                                                , marginTop = 1     &
                                                , marginBot = 1     &
                                                , newline = NLC     &
                                                , outputUnit = PM%LogFile%unit )
            end if

        end if blockMasterFileSetup

        ! These must be defined for all images, because they may be passed as arguments to the kernel subroutines.

        PM%LogFile%maxColWidth%val = max(PM%SpecBase%OutputRealPrecision%val, PM%SpecBase%OutputColumnWidth%val, PM%SpecBase%VariableNameList%MaxLen%val) + 9_IK
        PM%LogFile%maxColWidth%str = num2str(PM%LogFile%maxColWidth%val)
        PM%LogFile%format = getGenericFormat( width = PM%LogFile%maxColWidth%val &
                                            , precision = PM%SpecBase%OutputRealPrecision%val &
                                            , prefix = INDENT ) ! this is the generic indented format required mostly in postprocessing report
        PM%TimeFile%format = "(*(g" // PM%SpecBase%OutputColumnWidth%str // "." // PM%SpecBase%OutputRealPrecision%str // ",:,'" // PM%SpecBase%OutputDelimiter%val // "'))"
        PM%ChainFile%format = PM%TimeFile%format
        PM%SampleFile%format = PM%ChainFile%format
        PM%RestartFile%format = "(*(g0,:,'"//NLC//"'))"

    end subroutine setupOutputFiles

!***********************************************************************************************************************************
!***********************************************************************************************************************************

end module ParaMonte_mod