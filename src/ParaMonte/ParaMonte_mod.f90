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
    use ParaMonteChainFileContents_mod, only: ChainFileContents_type

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

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! ParaMonte IO variables and types
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! ParaMonte type
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type                                :: ParaMonte_type
        type(IntStr_type)               :: nd
        character(8)                    :: name
        character(16)                   :: brand
        character(:), allocatable       :: date
        character(:), allocatable       :: version
        logical                         :: isDryRun
        logical                         :: isFreshRun
        logical                         :: procArgNeeded
        logical                         :: procArgHasPriority
        logical                         :: inputFileArgIsPresent
        type(OS_type)                   :: OS
        type(Err_type)                  :: Err
        type(Image_type)                :: Image
        type(SpecBase_type)             :: SpecBase
        !type(ParaMonteStatistics_type) :: Stats
        type(SystemInfo_type)           :: SystemInfo
        type(Timer_type)                :: Timer
        type(File_type)                 :: InputFile
        type(LogFile_type)              :: LogFile
        type(TimeFile_type)             :: TimeFile
        type(ChainFile_type)            :: ChainFile
        type(SampleFile_type)           :: SampleFile
        type(RestartFile_type)          :: RestartFile
        type(ChainFileContents_type)    :: Chain
        type(Decoration_type)           :: Decor
    contains    
        procedure, pass                 :: reportDesc
        procedure, pass                 :: setupParaMonte
        procedure, pass                 :: addSplashScreen
        procedure, pass                 :: setupOutputFiles
        procedure, pass                 :: noteUserAboutEnvSetup
        procedure, pass                 :: addCompilerPlatformInfo
        procedure, pass                 :: warnUserAboutInputFilePresence
        procedure, pass                 :: setWarnAboutProcArgHasPriority
        procedure, nopass               :: informUser, note, warn, abort
        procedure, nopass               :: warnUserAboutMissingNamelist
    end type ParaMonte_type

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! To be called by all images.
    ! Tasks: setup initial variables, as well as construct default and null values for SpecBase. Sets
    ! self%InputFile%exists = .true. if the input file exists and opens and assigns to it a unit number and sets
    ! and self%InputFile%isOpen = .true. if the opening process is successful.
    ! If the input file exists, the path used to open it successfully will be also written to InpuFile%Path%modified
    !subroutine setupParaMonte(self,nd,name,date,version,inputFile)
    subroutine setupParaMonte(self,nd,name,inputFile)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setupParaMonte
#endif
        use, intrinsic :: iso_fortran_env, only: output_unit
        use Decoration_mod, only: INDENT
        use Constants_mod, only: IK, NLC
        use String_mod, only: getLowerCase, num2str
        use System_mod, only: OS_type
        implicit none
        class(ParaMonte_type), intent(inout)    :: self
        integer(IK), intent(in)                 :: nd
        character(*), intent(in)                :: name !, date, version
        character(*), intent(in), optional      :: inputFile
        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME // "@setupParaMonte()"

        self%Timer = Timer_type(self%Err)
        if (self%Err%occurred) then
            self%Err%msg = PROCEDURE_NAME // ": Error occurred while setting up the " // self%name // "timer."//NLC// self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if
        self%nd%val = nd
        self%LogFile%unit = output_unit   ! temporarily set the report file to stdout.
        self%Decor = Decoration_type()    ! initialize the TAB character and decoration symbol to the default values.

        self%Err%occurred = .false.
        self%Err%msg = ""

        self%name     = name
        self%brand    = INDENT // self%name
#if defined IFORT_ENABLED || __GFORTRAN__
        self%date     = "Build: " // __TIMESTAMP__
#else
        self%date     = "Unknown Release Date"
#endif
#if defined PARAMONTE_VERSION
        self%version  = "Version " // PARAMONTE_VERSION
#else
        self%version  = "Unknown Version"
#endif

        ! setup general processor / coarray image variables

#if defined CAF_ENABLED
        self%Image%id             = this_image()
        self%Image%count          = num_images()
#elif defined MPI_ENABLED
        block
            use mpi
            integer(IK) :: ierrMPI
            logical     :: isInitialized
            call mpi_initialized( isInitialized, ierrMPI )
            if (.not. isInitialized) call mpi_init(ierrMPI)
            call mpi_comm_rank(mpi_comm_world, self%Image%id, ierrMPI)
            call mpi_comm_size(mpi_comm_world, self%Image%count, ierrMPI)
            self%Image%id = self%Image%id + 1_IK ! make the ranks consistent with Fortran coarray indexing conventions
        end block
#else
        self%Image%id             = 1_IK
        self%Image%count          = 1_IK
#endif

        self%Image%name           = "@process(" // num2str(self%Image%id) // ")"
        self%Image%isFirst        = self%Image%id==1_IK
        self%Image%isNotFirst     = self%Image%id/=1_IK
        self%Image%isMaster       = .false.  ! ATTN: this will have to change later on, depending on the requested type of parallelism
        self%Image%isNotMaster    = .false.

        ! setup formatting variables

        self%nd%str = num2str(self%nd%val)

        ! determine OS. Should be only needed by the Master processes. But apparently not.

        call self%OS%query()
        if (self%OS%Err%occurred) then
            self%Err = self%OS%Err
            self%Err%msg = PROCEDURE_NAME // ": Error occurred while querying OS type."//NLC//self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if

        ! This is where SystemInfo used to live, but not anymore.

        blockSplashByFirstImage: if (self%Image%isFirst) then
            call self%addSplashScreen()
            call self%noteUserAboutEnvSetup()
        end if blockSplashByFirstImage

        ! check if input file exists by all images

        self%InputFile%isInternal = .false.
        self%inputFileArgIsPresent = present(inputFile)
        blockInputFileExistence: if (self%inputFileArgIsPresent) then
            self%InputFile = File_type( path=inputFile, status="old", OS=self%OS )
            if (self%InputFile%Err%occurred) then
                self%Err = self%InputFile%Err
                self%Err%msg = PROCEDURE_NAME//": Error occurred while attempting to setup the user's input file='"//inputFile//"'."//NLC//self%Err%msg
                call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                return
            end if
            ! determine if the file is internal
            self%InputFile%isInternal = self%inputFileArgIsPresent .and. .not.self%InputFile%exists .and. index(getLowerCase(inputFile),"&"//getLowerCase(self%name)) > 0
            if (.not.(self%InputFile%isInternal .or. self%InputFile%exists)) then
                ! file is given, but is neither a path to an external file, nor an internal file containing a namelist
                ! Therefore, there must be an error/mistake by the user.
                self%Err%msg =    "The user's input file='"//inputFile//"' is neither the path to an existing " // &
                                "external input file nor a string containing the input "//self%name//" specifications namelist. &
                                &This may be due to the user's mistake, providing a wrong path to the input external file or &
                                &a wrong list of input specifications for the " //self%name// " simulation. " //self%name//" will &
                                &assume no input specifications file is given by the user..."
                if (self%Image%isFirst) call self%warn(msg=self%Err%msg, prefix=self%brand, newline=NLC, outputUnit=self%LogFile%unit)
                self%inputFileArgIsPresent = .false.
            end if
        end if blockInputFileExistence

        if (self%Image%isFirst) call self%warnUserAboutInputFilePresence()

        ! Set the default and null values for ParaMonte SpecBase

        self%SpecBase = SpecBase_type(nd=self%nd%val,methodName=self%name,imageID=self%Image%id,imageCount=self%Image%count)

    end subroutine setupParaMonte

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine addSplashScreen(self)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: addSplashScreen
#endif
        implicit none
        class(ParaMonte_type), intent(inout) :: self

        self%Decor%text = &
        "\n\n"// &
        !self%name // "\n" // &
        "ParaMonte\n"// &
        "Plain Powerful Parallel\n"// &
        "Monte Carlo Library\n"// &
        "\n"// &
        self%version // "\n" // &
        "\n"// &
        self%date // "\n" // &
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

        call self%Decor%writeDecoratedText  ( text=self%Decor%text &
                                            , symbol="*" &
                                            , width=132 &
                                            , thicknessHorz=4 &
                                            , thicknessVert=2 &
                                            , marginTop=1 &
                                            , marginBot=2 &
                                            , outputUnit=self%LogFile%unit &
                                            , newLine="\n" &
                                            )

    end subroutine addSplashScreen

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine addCompilerPlatformInfo(self)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: addCompilerPlatformInfo
#endif
        use, intrinsic :: iso_fortran_env, only: compiler_version, compiler_options
        use Constants_mod, only: NLC
        implicit none
        class(ParaMonte_type), intent(inout)    :: self
        integer(IK)                             :: i, j

        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME//"@addCompilerPlatformInfo()"

        ! report the interface type to ParaMonte

        call self%Decor%writeDecoratedText  ( text="\nParaMonte library interface specifications\n" &
                                            , symbol="*" &
                                            , width=132 &
                                            , thicknessHorz=4 &
                                            , thicknessVert=1 &
                                            , marginTop=2 &
                                            , marginBot=1 &
                                            , outputUnit=self%LogFile%unit &
                                            , newLine="\n" &
                                            )
        self%Decor%List = self%Decor%wrapText( self%SpecBase%InterfaceType%val , 132 )
        do i = 1, size(self%Decor%List)
            write(self%LogFile%unit,"(*(g0))") self%Decor%List(i)%record
        end do

        ! report the ParaMonte compiler version and options

        call self%Decor%writeDecoratedText  ( text="\nParaMonte library compiler version\n" &
                                            , symbol="*" &
                                            , width=132 &
                                            , thicknessHorz=4 &
                                            , thicknessVert=1 &
                                            , marginTop=2 &
                                            , marginBot=1 &
                                            , outputUnit=self%LogFile%unit &
                                            , newLine="\n" &
                                            )
        self%Decor%List = self%Decor%wrapText( compiler_version() , 132 )
        do i = 1,size(self%Decor%List)
            write(self%LogFile%unit,"(*(g0))") self%Decor%List(i)%record
        end do

        call self%Decor%writeDecoratedText  (text="\nParaMonte library compiler options\n" &
                                            , symbol="*" &
                                            , width=132 &
                                            , thicknessHorz=4 &
                                            , thicknessVert=1 &
                                            , marginTop=2 &
                                            , marginBot=1 &
                                            , outputUnit=self%LogFile%unit &
                                            , newLine="\n" &
                                            )
        self%Decor%List = self%Decor%wrapText( compiler_options() , 132 )
        do i = 1,size(self%Decor%List)
            write(self%LogFile%unit,"(*(g0))") self%Decor%List(i)%record
        end do

        call self%Decor%writeDecoratedText  ( text="\nRuntime platform specifications\n" &
                                            , symbol="*" &
                                            , width=132 &
                                            , thicknessHorz=4 &
                                            , thicknessVert=1 &
                                            , marginTop=2 &
                                            , marginBot=1 &
                                            , outputUnit=self%LogFile%unit &
                                            , newLine="\n" &
                                            )

        ! Get system info by all images. why? Not anymore.
        ! On many parallel processors via singlChain this leads to 
        ! the creation of thousands of files on the system, simultaneously.
        ! this is not needed by any process other than the masters.

        if (allocated(self%SpecBase%SystemInfoFilePath%val)) then
            block
                use FileContents_mod, only: FileContents_type
                type(FileContents_type) :: FileContents
                FileContents = FileContents_type(filePath = self%SpecBase%SystemInfoFilePath%val)
                if (FileContents%Err%occurred) then
                    self%Err = FileContents%Err
                    self%Err%msg = PROCEDURE_NAME//": Error occurred while collecting system info."//NLC//self%Err%msg
                    call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                    return
                else
                    do j = 1, FileContents%numRecord
                        self%Decor%List = self%Decor%wrapText( FileContents%Line(j)%record , 132 )
                        do i = 1,size(self%Decor%List)
                            write(self%LogFile%unit,"(A)") self%Decor%List(i)%record
                        end do
                    end do
                    deallocate(self%SpecBase%SystemInfoFilePath%val)
                end if
            end block
        else
            self%SystemInfo = SystemInfo_type(OS=self%OS)
            if (self%SystemInfo%Err%occurred) then
                self%Err = self%SystemInfo%Err
                self%Err%msg = PROCEDURE_NAME//": Error occurred while collecting system info."//NLC//self%Err%msg
                call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                return
            end if
            do j = 1, self%SystemInfo%nRecord
                self%Decor%List = self%Decor%wrapText( self%SystemInfo%List(j)%record , 132 )
                do i = 1,size(self%Decor%List)
                    write(self%LogFile%unit,"(A)") self%Decor%List(i)%record
                end do
            end do
        end if
        call self%Decor%write(self%LogFile%unit)

    end subroutine addCompilerPlatformInfo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine noteUserAboutEnvSetup(self)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: noteUserAboutEnvSetup
#endif
        implicit none
        class(ParaMonte_type), intent(inout) :: self
        call self%Decor%writeDecoratedText  ( text = "\nSetting up the " // self%name // " simulation environment\n" &
                                            , marginTop = 1     &
                                            , marginBot = 1     &
                                            , newline = "\n"    &
                                            , outputUnit = self%LogFile%unit )
    end subroutine noteUserAboutEnvSetup

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine warnUserAboutInputFilePresence(self)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: warnUserAboutInputFilePresence
#endif
        use Constants_mod, only: NLC
        implicit none
        class(ParaMonte_type), intent(inout) :: self
#if defined CFI_ENABLED
        if (self%SpecBase%InterfaceType%isPython) then
            call self%note  ( prefix     = self%brand &
                            , outputUnit = self%LogFile%unit &
                            , newline    = NLC &
                            , msg        = "Interfacing Python with "// self%name //"..." )
        elseif (self%SpecBase%InterfaceType%isClang) then
#else
            if (self%inputFileArgIsPresent) then
                if (self%InputFile%exists) then
                    call self%note  ( prefix     = self%brand &
                                    , outputUnit = self%LogFile%unit &
                                    , newline    = NLC &
                                    , msg        = "The user's input file for " // self%name // " options was detected."// NLC // &
                                                   "All " // self%name // " options will be read from the input file."// NLC // &
                                                   "Here is " // self%name // " input options file:"// NLC // self%inputFile%Path%modified )
                elseif (self%InputFile%isInternal) then
                    call self%note  ( prefix     = self%brand &
                                    , outputUnit = self%LogFile%unit &
                                    , newline    = NLC &
                                    , msg        = "No external file corresponding to the user's input file for "//self%name//" options &
                                                   &could be found."//NLC//"The user-provided input file will be processed as an input &
                                                   &string of "//self%name//" options." )
                end if
            else
                call self%note  ( prefix     = self%brand &
                                , outputUnit = self%LogFile%unit &
                                , newline    = NLC &
                                , msg        = "No " // self%name  // " input file is provided by the user."//NLC//&
                                               "Variable values from the procedure arguments will be used instead, where provided."//NLC//&
                                               "Otherwise, the default options will be used." )
            end if
#endif
#if defined CFI_ENABLED
        end if
#endif
    end subroutine warnUserAboutInputFilePresence

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setWarnAboutProcArgHasPriority(self)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setWarnAboutProcArgHasPriority
#endif
        implicit none
        class(ParaMonte_type), intent(inout) :: self
#if !defined CFI_ENABLED
        character(:), allocatable            :: msg
#endif
        self%procArgHasPriority = .not. self%SpecBase%InputFileHasPriority%val
        self%procArgNeeded = self%procArgHasPriority .or. (.not.self%inputFileArgIsPresent)

#if defined FORTRAN_ENABLED
        if (self%Image%isFirst) then 
            if (self%procArgHasPriority) then 
                msg =   "Variable inputFileHasPriority = .false.\n&
                        &All variable values will be overwritten by the corresponding procedure argument values,\n&
                        &only if provided as procedure arguments."
            else
                msg =   "Variable inputFileHasPriority = .true.\n&
                        &All variable values will be read from the user-provided input file"
            end if
            call self%note( prefix = self%brand, outputUnit = self%LogFile%unit, newline = "\n", msg = msg )
        end if
#endif
    end subroutine setWarnAboutProcArgHasPriority

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine setupOutputFiles(self)
#if defined DLL_ENABLED && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: setupOutputFiles
#endif
        use Decoration_mod, only: getGenericFormat, INDENT
        use Constants_mod, only: NLC, FILE_EXT
        use String_mod, only: num2str
        use Path_mod, only: MAX_FILE_PATH_LEN, mkdir

        implicit none

        class(ParaMonte_type), intent(inout)    :: self
        character(:), allocatable               :: msg, workingOn, currentWorkingDir
        character(*), parameter                 :: PROCEDURE_NAME = MODULE_NAME // "@setupOutputFiles()"

        if (self%SpecBase%OutputFileName%original==self%SpecBase%OutputFileName%def) then
            msg =   "No user-input filename prefix for " // self%name // " output files detected." // NLC // &
                    "Generating appropriate filenames for " // self%name // " output files from the current date and time..."
        else
            msg = "Variable outputFileName detected among the input variables to " // self%name // ":" //NLC// self%SpecBase%OutputFileName%original
        end if

        call self%SpecBase%OutputFileName%query(OS=self%OS)

        if (self%SpecBase%OutputFileName%Err%occurred) then
            self%Err = self%SpecBase%OutputFileName%Err
            self%Err%msg = PROCEDURE_NAME // ": Error occurred while attempting to construct OutputFileName path type." //NLC// self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if

        self%SpecBase%OutputFileName%namePrefix = self%SpecBase%OutputFileName%name // self%SpecBase%OutputFileName%ext

        ! get the current working directory

        if (allocated(currentWorkingDir)) deallocate(currentWorkingDir)
        allocate( character(MAX_FILE_PATH_LEN) :: currentWorkingDir )
        block
#if defined IFORT_ENABLED
            use ifport ! only: getcwd
#endif
#if defined IFORT_ENABLED || __GFORTRAN__
            self%Err%stat = getcwd(currentWorkingDir)
            currentWorkingDir = trim(adjustl(currentWorkingDir))
#else
            self%Err%stat = 0_IK
            currentWorkingDir = "."
#endif
        end block
        if (self%Err%stat/=0) then
            self%Err%msg = PROCEDURE_NAME//": Error occurred while fetching the current working directory via getcwd()."//NLC
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if
        msg = msg //NLC//NLC// "Absolute path to the current working directory:"//NLC//currentWorkingDir

        if (len_trim(adjustl(self%SpecBase%OutputFileName%dir))==0) then
            self%SpecBase%OutputFileName%dir = trim(adjustl(currentWorkingDir)) // self%SpecBase%OutputFileName%slashOS
            msg = msg //NLC//NLC// "All output files will be written to the current working directory:"//NLC//self%SpecBase%OutputFileName%dir
        else
            msg = msg //NLC//NLC// "Generating the requested directory for the "//self%name//" output files:"//NLC//self%SpecBase%OutputFileName%dir
        end if

        ! Generate the output files directory:

        if (self%Image%isFirst) then
            self%Err = mkdir( dirPath = self%SpecBase%OutputFileName%dir, isWindows = self%OS%isWindows )
            if (self%Err%occurred) then
                self%Err%msg = PROCEDURE_NAME//": Error occurred while making directory = '"//self%SpecBase%OutputFileName%dir//"'."//NLC//self%Err%msg
                call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
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

        if (len_trim(adjustl(self%SpecBase%OutputFileName%namePrefix))==0) then
            msg = msg //NLC//NLC// "No user-input filename prefix for " // self%name // " output files detected."//NLC//&
            "Generating appropriate filenames for " // self%name // " output files from the current date and time..."
            self%SpecBase%OutputFileName%namePrefix = self%SpecBase%OutputFileName%def
        end if

        self%SpecBase%OutputFileName%pathPrefix = self%SpecBase%OutputFileName%dir // self%SpecBase%OutputFileName%namePrefix

        ! Variable msg will be used down this subroutine, so it should not be changed beyond this point
        msg  =  msg //NLC//NLC// self%name // " output files will be prefixed with:"//NLC// self%SpecBase%OutputFileName%pathPrefix
        if (self%Image%isFirst) then
            call self%note( prefix = self%brand, outputUnit = self%LogFile%unit, newline = NLC, msg = msg )
        end if

        ! Generate the output filenames, search for pre-existing runs, and open the report file:

        if (self%Image%isFirst) call self%note( prefix = self%brand, outputUnit = self%LogFile%unit, newline = NLC, msg = "Searching for previous runs of " // self%name // "..." )

        ! this block could be all executed by only the master images

        self%LogFile%Path%Ext = FILE_EXT%ascii
        self%TimeFile%Path%Ext = FILE_EXT%ascii
        self%ChainFile%Path%Ext = FILE_EXT%ascii
        self%SampleFile%Path%Ext = FILE_EXT%ascii
        self%RestartFile%Path%Ext = FILE_EXT%ascii

        self%RestartFile%Form%value = "formatted"
        if (self%SpecBase%RestartFileFormat%isBinary) then
            self%RestartFile%Form%value = "unformatted"
            self%RestartFile%Path%Ext = FILE_EXT%binary
        end if

        self%ChainFile%Form%value = "formatted"
        if (self%SpecBase%ChainFileFormat%isBinary) then
            self%ChainFile%Form%value = "unformatted"
            self%ChainFile%Path%Ext = FILE_EXT%binary
        end if

        block
            use Path_mod, only: Path_type
            integer(IK) :: imageID
#if defined CAF_ENABLED || defined MPI_ENABLED
            if (self%SpecBase%ParallelizationModel%isMultiChain) then
                imageID = self%Image%id
            else
#endif
                imageID = 1_IK
#if defined CAF_ENABLED || defined MPI_ENABLED
            end if
#endif
            self%SpecBase%OutputFileName%pathPrefix = self%SpecBase%OutputFileName%pathPrefix // "_process_" // num2str(imageID) // "_"
            self%LogFile%Path       = Path_type( inputPath = self%SpecBase%OutputFileName%pathPrefix // self%LogFile%suffix        // self%LogFile%Path%Ext      , OS = self%OS )
            self%TimeFile%Path      = Path_type( inputPath = self%SpecBase%OutputFileName%pathPrefix // self%TimeFile%suffix       // self%TimeFile%Path%Ext     , OS = self%OS )
            self%ChainFile%Path     = Path_type( inputPath = self%SpecBase%OutputFileName%pathPrefix // self%ChainFile%suffix      // self%ChainFile%Path%Ext    , OS = self%OS )
            self%SampleFile%Path    = Path_type( inputPath = self%SpecBase%OutputFileName%pathPrefix // self%SampleFile%suffix     // self%SampleFile%Path%Ext   , OS = self%OS )
            self%RestartFile%Path   = Path_type( inputPath = self%SpecBase%OutputFileName%pathPrefix // self%RestartFile%suffix    // self%RestartFile%Path%Ext  , OS = self%OS )
        end block

        inquire( file = self%LogFile%Path%original, exist = self%LogFile%exists, iostat = self%LogFile%Err%stat )
        self%Err = self%LogFile%getInqErr( self%LogFile%Err%stat )
        if (self%Err%occurred) then
            self%Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file='" // self%LogFile%Path%original // self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if

        inquire( file = self%SampleFile%Path%original, exist = self%SampleFile%exists, iostat = self%SampleFile%Err%stat )
        self%Err = self%SampleFile%getInqErr( self%SampleFile%Err%stat )
        if (self%Err%occurred) then
            self%Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file='" // self%SampleFile%Path%original // self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if

        inquire( file = self%TimeFile%Path%original, exist = self%TimeFile%exists, iostat = self%TimeFile%Err%stat )
        self%Err = self%TimeFile%getInqErr( self%TimeFile%Err%stat )
        if (self%Err%occurred) then
            self%Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file='" // self%TimeFile%Path%original // self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if

        inquire( file = self%ChainFile%Path%original, exist = self%ChainFile%exists, iostat = self%ChainFile%Err%stat )
        self%Err = self%ChainFile%getInqErr( self%ChainFile%Err%stat )
        if (self%Err%occurred) then
            self%Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file='" // self%ChainFile%Path%original // self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if

        inquire( file = self%RestartFile%Path%original, exist = self%RestartFile%exists, iostat = self%RestartFile%Err%stat )
        self%Err = self%RestartFile%getInqErr( self%RestartFile%Err%stat )
        if (self%Err%occurred) then
            self%Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file='" // self%RestartFile%Path%original // self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if

        self%isDryRun = (.not. self%SpecBase%OverwriteRequested%val) .and. & ! not fresh, if any file exists
                        (self%LogFile%exists .or. self%TimeFile%exists .or. self%RestartFile%exists .or. self%ChainFile%exists .or. self%SampleFile%exists)
        self%isFreshRun = .not. self%isDryRun

        if (self%isFreshRun) then
            if (self%Image%isFirst) call self%note( prefix = self%brand, outputUnit = self%LogFile%unit, newline = NLC, msg = "No pre-existing "//self%name//" run detected."//NLC//"Starting a fresh "//self%name//" run..." )
        else
            if (self%Image%isFirst) call self%note( prefix = self%brand, outputUnit = self%LogFile%unit, newline = NLC, msg = "Previous run of "//self%name//" detected."//NLC//"Searching for restart files..." )
            if (self%SampleFile%exists) then ! sampling is already complete
                self%Err%occurred = .true.
                self%Err%msg =    PROCEDURE_NAME//": Error occurred. Output sample file detected: "//self%SampleFile%Path%original//&
                                NLC//self%name//" cannot overwrite an already-completed simulation."//&
                                NLC//"Please provide an alternative file name for the new simulation outputs."
            elseif (self%LogFile%exists .and. self%TimeFile%exists .and. self%RestartFile%exists .and. self%ChainFile%exists) then  ! restart mode
                if (self%SpecBase%SampleSize%val==0_IK) then ! sampling is already complete
                    self%Err%occurred = .true.
                    self%Err%msg = PROCEDURE_NAME//": Error occurred. The input variable sampleSize = 0 indicates that the output files belong to a completed simulation."
                else
                    self%Err%occurred = .false.
                end if
            else
                self%Err%occurred = .true.
                self%Err%msg = PROCEDURE_NAME//": Error occurred. For a successful simulation restart, all output files are necessary."//NLC//"List of missing simulation output files:"
                if (.not. self%LogFile%exists)        self%Err%msg = self%Err%msg//NLC//self%LogFile%Path%original
                if (.not. self%TimeFile%exists)       self%Err%msg = self%Err%msg//NLC//self%TimeFile%Path%original
                if (.not. self%ChainFile%exists)      self%Err%msg = self%Err%msg//NLC//self%ChainFile%Path%original
                if (.not. self%RestartFile%exists)    self%Err%msg = self%Err%msg//NLC//self%RestartFile%Path%original
            end if
            if (self%Err%occurred) then
                call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                return
            end if
        end if

        ! open/append the output files:

        if (self%Image%isMaster) then
            if (self%isFreshRun) then
                workingOn = "Generating the output "
                self%LogFile%status = "replace"
                self%TimeFile%status = "replace"
                self%ChainFile%status = "replace"
                self%SampleFile%status = "replace"
                self%RestartFile%status = "replace"
                self%LogFile%Position%value = "asis"
                self%TimeFile%Position%value = "asis"
                self%ChainFile%Position%value = "asis"
                self%SampleFile%Position%value = "asis"
                self%RestartFile%Position%value = "asis"
            else
                workingOn = "Appending to the existing "
                self%LogFile%status = "old"
                self%TimeFile%status = "old"
                self%ChainFile%status = "old"
                self%SampleFile%status = "replace"
                self%RestartFile%status = "old"
                self%LogFile%Position%value = "append"
                self%TimeFile%Position%value = "asis"
                self%ChainFile%Position%value = "asis"
                self%SampleFile%Position%value = "asis"
                self%RestartFile%Position%value = "asis"
            end if
        end if

        ! print the stdout message for generating / appending the output report file

        blockLogFileListByFirstImage: if (self%Image%isFirst) then

            ! print the stdout message for generating / appending the output report file(s)

            call self%note  ( prefix = self%brand               &
                            , outputUnit = self%LogFile%unit    &
                            , newline = NLC                     &
                            , marginBot = 0_IK                  &
                            , msg = workingOn // self%LogFile%suffix // " file:" )

            ! print the the output report file name of the images

            call self%note  ( prefix = self%brand               &
                            , outputUnit = self%LogFile%unit    &
                            , newline = NLC                     &
                            , marginTop = 0_IK                  &
                            , marginBot = 0_IK                  &
                            , msg = self%LogFile%Path%original  )

#if defined CAF_ENABLED || defined MPI_ENABLED
            if (self%SpecBase%ParallelizationModel%isMultiChain) then
                block
                    use String_mod, only: replaceStr !, num2str
                    integer(IK) :: imageID
                    do imageID = 2, self%Image%count
                        call self%note  ( prefix = self%brand               &
                                        , outputUnit = self%LogFile%unit    &
                                        , newline = NLC                     &
                                        , marginTop = 0_IK                  &
                                        , marginBot = 0_IK                  &
                                        , msg = replaceStr( string = self%LogFile%Path%original, search = "process_1", substitute = "process_"//num2str(imageID) ) )
                    end do
                end block
            end if
            self%Err%msg = "Running the simulation in parallel on " // num2str(self%Image%count) // " processes." // NLC
#else
            self%Err%msg = "Running the simulation in serial on " // num2str(self%Image%count) // " process." // NLC
#endif

            call self%note  ( prefix = self%brand               &
                            , outputUnit = self%LogFile%unit    &
                            , newline = NLC                     &
                            , marginTop = 3_IK                  &
                            , marginBot = 3_IK                  &
                            , msg = self%Err%msg // "Please see the output " // self%LogFile%suffix // " and " // self%TimeFile%suffix // " files for further realtime simulation details." &
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

        blockMasterFileSetup: if (self%Image%isMaster) then

            self%LogFile%unit = 1001 + self%Image%id ! for some unknown reason, if newunit is used, GFortran opens the file as an internal file
            open( unit = self%LogFile%unit              &
                , file = self%LogFile%Path%original     &
                , status = self%LogFile%status          &
                , iostat = self%LogFile%Err%stat        &
#if defined IFORT_ENABLED && defined OS_IS_WINDOWS
                , SHARED                                &
#endif
                , position = self%LogFile%Position%value)
            self%Err = self%LogFile%getOpenErr(self%LogFile%Err%stat)
            if (self%Err%occurred) then
                self%Err%msg = PROCEDURE_NAME // ": Error occurred while opening the " // self%name // " " // self%LogFile%suffix // " file='" // self%LogFile%Path%original // "'. "
                if (scan(" ",trim(adjustl(self%LogFile%Path%original)))/=0) then
                    self%Err%msg = self%Err%msg // "It appears that absolute path used for the output files contains whitespace characters. " &
                                            // "This could be one potential cause of the simulation failure. " &
                                            // "The whitespace characters are always problematic in paths. " &
                                            // "Ensure the path used for the output files does not contain whitespace characters. "
                end if
                call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                return
            end if

            ! rewrite the same old stuff to all report files

            if (self%isFreshRun) then
                call self%addSplashScreen()
                if (self%SpecBase%SilentModeRequested%isFalse) call self%addCompilerPlatformInfo()   ! this takes about 0.75 seconds to execute on Stampede Login nodes.
                call self%noteUserAboutEnvSetup()
                call self%warnUserAboutInputFilePresence()
                call self%setWarnAboutProcArgHasPriority()
                call self%note( prefix = self%brand, outputUnit = self%LogFile%unit, newline = NLC, msg = msg )
            end if

            ! open/append the output files

            if (self%isFreshRun) call self%note( prefix = self%brand, outputUnit = self%LogFile%unit, newline = NLC, msg = workingOn//self%TimeFile%suffix//" file:"//NLC//self%TimeFile%Path%original )

            self%TimeFile%unit = 1000001  ! for some unknown reason, if newunit is used, GFortran opens the file as an internal file
            open( unit = self%TimeFile%unit                 &
                , file = self%TimeFile%Path%original        &
                , status = self%TimeFile%status             &
                , iostat = self%TimeFile%Err%stat           &
#if defined IFORT_ENABLED && defined OS_IS_WINDOWS
                , SHARED                                    &
#endif
                , position = self%TimeFile%Position%value   )
            self%Err = self%TimeFile%getOpenErr(self%TimeFile%Err%stat)
            if (self%Err%occurred) then
                self%Err%msg = PROCEDURE_NAME // ": Error occurred while opening the " // self%name // " " // self%TimeFile%suffix // " file='" // self%TimeFile%Path%original // "'. "
                call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                return
            end if

            if (self%isFreshRun) call self%note( prefix = self%brand, outputUnit = self%LogFile%unit, newline = NLC, msg = workingOn//self%ChainFile%suffix//"file:"//NLC//self%ChainFile%Path%original )

            self%ChainFile%unit = 2000001  ! for some unknown reason, if newunit is used, GFortran opens the file as an internal file
            open( unit = self%ChainFile%unit                &
                , file = self%ChainFile%Path%original       &
                , form = self%ChainFile%Form%value          &
                , status = self%ChainFile%status            &
                , iostat = self%ChainFile%Err%stat          &
#if defined IFORT_ENABLED && defined OS_IS_WINDOWS
                , SHARED                                    &
#endif
                , position = self%ChainFile%Position%value  )
            self%Err = self%ChainFile%getOpenErr(self%ChainFile%Err%stat)
            if (self%Err%occurred) then
                self%Err%msg = PROCEDURE_NAME // ": Error occurred while opening the " // self%name // " " // self%ChainFile%suffix // " file='" // self%ChainFile%Path%original // "'. "
                call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                return
            end if

            self%RestartFile%unit = 3000001  ! for some unknown reason, if newunit is used, GFortran opens the file as an internal file
            open( unit = self%RestartFile%unit              &
                , file = self%RestartFile%Path%original     &
                , form = self%RestartFile%Form%value        &
                , status = self%RestartFile%status          &
                , iostat = self%RestartFile%Err%stat        &
#if defined IFORT_ENABLED && defined OS_IS_WINDOWS
                , SHARED                                    &
#endif
                , position = self%RestartFile%Position%value)
            self%Err = self%RestartFile%getOpenErr(self%RestartFile%Err%stat)
            if (self%Err%occurred) then
                self%Err%msg = PROCEDURE_NAME // ": Error occurred while opening the " // self%name // " " // self%RestartFile%suffix // " file='" // self%RestartFile%Path%original // "'. "
                call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                return
            end if

            if (self%isFreshRun) then
                call self%Decor%writeDecoratedText  ( text = NLC // self%name // " simulation specifications" // NLC &
                                                    , marginTop = 1     &
                                                    , marginBot = 1     &
                                                    , newline = NLC     &
                                                    , outputUnit = self%LogFile%unit )
            end if

        end if blockMasterFileSetup

        ! These must be defined for all images, because they may be passed as arguments to the kernel subroutines.

        self%LogFile%maxColWidth%val = max(self%SpecBase%OutputRealPrecision%val, self%SpecBase%OutputColumnWidth%val, self%SpecBase%VariableNameList%MaxLen%val) + 9_IK
        self%LogFile%maxColWidth%str = num2str(self%LogFile%maxColWidth%val)
        self%LogFile%format = getGenericFormat( width = self%LogFile%maxColWidth%val &
                                            , precision = self%SpecBase%OutputRealPrecision%val &
                                            , prefix = INDENT ) ! this is the generic indented format required mostly in postprocessing report
        self%TimeFile%format = "(*(g" // self%SpecBase%OutputColumnWidth%str // "." // self%SpecBase%OutputRealPrecision%str // ",:,'" // self%SpecBase%OutputDelimiter%val // "'))"
        self%ChainFile%format = self%TimeFile%format
        self%SampleFile%format = self%ChainFile%format
        self%RestartFile%format = "(*(g0,:,'"//NLC//"'))"

    end subroutine setupOutputFiles

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subroutine reportDesc(self, msg) !, marginTop, marginBot)
#if defined DLL_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: reportDesc
#endif
        use Constants_mod, only: IK, NLC
        implicit none
        class(ParaMonte_type), intent(inout)    :: self
        character(*), intent(in)                :: msg
        !integer(IK) , intent(in)                :: marginTop, marginBot
        !integer(IK)                             :: marginTopDef, marginBotDef
        call self%note  ( prefix     = self%brand           &
                        , outputUnit = self%LogFile%unit    &
                        , newline    = NLC                  &
                        , marginTop  = 1_IK                 &
                        , marginBot  = 1_IK                 &
                        , msg        = msg                  )
    end subroutine reportDesc

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module ParaMonte_mod