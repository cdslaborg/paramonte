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

!>  \mainpage ParaMonte: Plain Powerful <b>Para</b>llel <b>Monte</b> Carlo Library
!>
!>  This is the <b>`ParaMonte::Kernel`</b> developer documentation website.
!>
!>  What is ParaMonte?
!>  ==================
!>
!>  ParaMonte is a serial / parallel library of Monte Carlo routines for sampling
!>  mathematical objective functions of arbitrary-dimensions, in particular,
!>  the posterior distributions of Bayesian models in data science,
!>  Machine Learning, and scientific inference, with the design
!>  goal of unifying the
!>
!>  +    **automation** of Monte Carlo simulations,
!>  +    **user-friendliness** of the library,
!>  +    **accessibility** from multiple programming environments,
!>  +    **high-performance** at runtime, and,
!>  +    **scalability** across many parallel processors.
!>
!>  ### ParaMonte project's repository
!>
!>  The ParaMonte library is open-source is permanently located and maintained on **GitHub** at:
!>
!>  &nbsp;&nbsp;&nbsp;&nbsp;[**https://github.com/cdslaborg/paramonte**](https://github.com/cdslaborg/paramonte)
!>
!>  ### ParaMonte usage and examples website
!>
!>  For information about the usage and examples visit **the ParaMonte documentation and examples website** at:
!>
!>  &nbsp;&nbsp;&nbsp;&nbsp;[**https://www.cdslab.org/paramonte**](https://www.cdslab.org/paramonte)
!>
!>  ### ParaMonte API documentation website
!>
!>  For the API developer documentation, visit:
!>
!>  &nbsp;&nbsp;&nbsp;&nbsp;[**https://www.cdslab.org/paramonte/notes/api/kernel**](https://www.cdslab.org/paramonte/notes/api/kernel)
!>
!>  ParaMonte samplers
!>  ==================
!>
!>  The routines currently supported by the ParaMonte kernel library include:
!>
!>  ### ParaDRAM
!>
!>  Parallel Delayed-Rejection Adaptive Metropolis-Hastings Markov
!>  Chain Monte Carlo Sampler. For a quick start, example scripts,
!>  and instructions on how to use he ParaDRAM sampler in your
!>  language of choice, visit:
!>
!>  &nbsp;&nbsp;&nbsp;&nbsp;[**https://www.cdslab.org/paramonte/notes/usage/paradram/interface**](https://www.cdslab.org/paramonte/notes/usage/paradram/interface)
!>
!>  Naming conventions
!>  ==================
!>
!>  +   The CamelCase naming style is used throughout the entire ParaMonte
!>      kernel library.
!>
!>  +   Although the Fortran language is case-insensitive, by convention,
!>      all scalar variable names begin with a lower case, whereas all vectors,
!>      arrays, types, and module names begin with an upper-case letter.
!>
!>  +   The name of any variable that represents a vector of values is normally
!>      suffixed with `Vec` or `Vector`, for example: `StartPointVec`, ...
!>
!>  +   The name of any variable that represents a matrix of values is normally
!>      suffixed with `Mat`, for example: `proposalStartCorMat`, ...
!>
!>  +   The name of any variable that represents a list of varying-size values
!>      is normally suffixed with `List`, like: `variableNameList`, ...
!>
!>  +   All static functions or methods of classes begin with a lowercase verb.
!>
!>  +   Significant attempt has been made to end all boolean variables with a
!>      passive verb, such that the full variable name virtually forms a
!>      proposition, that is, an English-language statement that should
!>      be either `.true.` or `.false.`, set by the user.
!>
!-------------------------------------------------------------------------------

!>  \brief This module contains the base class of all ParaMonte samplers and its associated methods.
!>  @author Amir Shahmoradi

module ParaMonte_mod

    use System_mod, only: SystemInfo_type
    use Parallelism_mod, only: Image_type
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

    !> The Quantile derived type containing the distribution quantiles.
    type                            :: QuantileProbability_type
        integer(IK)                 :: count = 9_IK
        real(RK)                    :: Value(9) = [0._RK,0.05_RK,0.10_RK,0.25_RK,0.50_RK,0.75_RK,0.90_RK,0.95_RK,1.0_RK]
        character(4)                :: Name(9) = ["  Q0","  Q5"," Q10"," Q25"," Q50"," Q75"," Q90"," Q95","Q100"]
    end type QuantileProbability_type
    type(QuantileProbability_type), parameter :: QPROB = QuantileProbability_type()

    !> The derived type containing the statistical moments of the objective function.
    type                            :: Moment_type
        integer(IK)                 :: count = 0_IK
        real(RK), allocatable       :: Mean(:)
        real(RK), allocatable       :: CovMat(:,:)
        real(RK), allocatable       :: CorMat(:,:)
        real(RK), allocatable       :: Quantile(:,:)
    end type Moment_type

    !> The derived type containing the number of function calls.
    type                            :: ParaMonteNumFunCall_type
        integer(IK)                 :: accepted         !< The number of objective function calls accepted in the simulation.
        integer(IK)                 :: acceptedRejected !< The accepted + rejected function calls.
    end type ParaMonteNumFunCall_type

    !> The derived type containing information about the function mode.
    type                            :: ParaMonteLogFuncMode_type
        real(RK)                    :: val      !< The function mode.
        real(RK), allocatable       :: Crd(:)   !< The location of the mode in the domain of the the objective function.
    end type ParaMonteLogFuncMode_type

    !> The derived type containing information about the statistics of the objective function and the runtime performance of the simulation.
    type                            :: ParaMonteStatistics_type
        real(RK)                    :: avgTimePerFunCalInSec = 0._RK    !< Average time per objective function call.
        real(RK)                    :: avgCommTimePerFunCall = 0._RK    !< Average inter-process communication time per function call.
        type(Moment_type)           :: Sample                           !< The statistical moments of the objective function.
    end type ParaMonteStatistics_type

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! ParaMonte IO variables and types
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

    !> The `ParaMonte_type` sampler base class.
    type                                :: ParaMonte_type
        type(IntStr_type)               :: nd                       !< The number of dimensions of the domain of the objective function.
        character(8)                    :: name                     !< The ParaMonte sampler method name.
        character(16)                   :: brand                    !< The ParaMonte sampler brand (The decorated tabbed version of the sampler name).
        character(:), allocatable       :: date                     !< The date of the simulation.
        character(:), allocatable       :: version                  !< The ParaMonte version.
        logical                         :: isDryRun                 !< The logical flag that, if `.true.`, indicates the simulation is in restart mode.
        logical                         :: isFreshRun               !< The logical flag that, if `.false.`, indicates the simulation is in new mode.
        logical                         :: procArgNeeded            !< The logical flag that, if `.true.`, requires reading the simulation specification
                                                                    !< from the Object-Oriented interface of the sampler and prioritizing them over the
                                                                    !< corresponding values in the external input file.
        logical                         :: procArgHasPriority       !< The logical flag that, if `.true.`, indicates that the simulation specifications
                                                                    !< have priority over the corresponding values from the external input file.
        logical                         :: inputFileArgIsPresent    !< The logical flag that indicates whether the external input file has been provided by the user.
        type(OS_type)                   :: OS                       !< An object of class [OS_type](@ref system_mod::os_type) containing information about the Operating System.
        type(Err_type)                  :: Err                      !< An object of class [Err_type](@ref err_mod::err_type) containing error-handling information.
                                                                    !< about error occurrence and message during the simulation setup and runtime.
        type(Image_type)                :: Image                    !< An object of type [Image_type](@ref image_type) containing information about
                                                                    !< the processor count and types in the simulation.
        type(SpecBase_type)             :: SpecBase                 !< An object of class [SpecBase_type](@ref specbase_mod::specbase_type) containing information
                                                                    !< about the basic simulation specification properties.
        !type(ParaMonteStatistics_type) :: Stats
        type(SystemInfo_type)           :: SystemInfo               !< An object of class [SystemInfo_type](@ref system_mod::systeminfo_type) containing
                                                                    !< information about the operating system and platform.
        type(Timer_type)                :: Timer                    !< An object of class [Timer_type](@ref timer_mod::timer_type) used for timing of the simulation.
        type(File_type)                 :: InputFile                !< An object of class [File_type](@ref file_mod::file_type) containing information about the simulation input file.
        type(LogFile_type)              :: LogFile                  !< An object of class [LogFile_type](@ref logfile_type) containing information about the simulation report file.
        type(TimeFile_type)             :: TimeFile                 !< An object of class [TimeFile_type](@ref timefile_type) containing information about the simulation timing.
        type(ChainFile_type)            :: ChainFile                !< An object of class [ChainFile_type](@ref chainfile_type) containing information about the simulation output chain.
        type(SampleFile_type)           :: SampleFile               !< An object of class [SampleFile_type](@ref samplefile_type) containing information about the simulation output sample.
        type(RestartFile_type)          :: RestartFile              !< An object of class [RestartFile_type](@ref restartfile_type) containing information about the simulation output restart.
        type(ChainFileContents_type)    :: Chain                    !< An object of class [ChainFileContents_type](@ref paramontechainfilecontents_mod::chainfilecontents_type) containing information and methods for chain IO.
        type(Decoration_type)           :: Decor                    !< An object of class [Decoration_type](@ref decoration_mod::decoration_type) containing IO decoration tools.
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

    !> \brief
    !> This procedure is a method of the [ParaMonte_type](@ref paramonte_type) class.
    !> Set up the ParaMonte sampler with the requested input specifications. This method,
    !> + sets up initial variables,
    !> + and constructs the default and null values for `SpecBase`.
    !> + `self%InputFile%exists = .true.` if the input file exists and opens and assigns to it a unit number and sets
    !>    and `self%InputFile%isOpen = .true.` if the opening process is successful.
    !> + If the input file exists, the path used to open it successfully will be also written to `InpuFile%Path%modified`.
    !>
    !> @param[inout]    self        :   An object of class [ParaMonte_type](@ref paramonte_type).
    !> @param[in]       nd          :   The number of dimensions of the domain of the objective function.
    !> @param[in]       name        :   The name of the sampler. Example: `ParaDRAM`.
    !> @param[in]       inputFile   :   The path to the input file, or the contents of an input file.
    !>
    !> \warning
    !> This routine has to be called by all images (processes).
    subroutine setupParaMonte(self,nd,name,inputFile)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
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

        self%LogFile%unit = output_unit ! temporarily set the report file to stdout.
        self%Err%occurred = .false.
        self%Err%msg = ""

        self%Timer = Timer_type(self%Err)
        if (self%Err%occurred) then
        ! LCOV_EXCL_START
            self%Err%msg = PROCEDURE_NAME // ": Error occurred while setting up the " // self%name // "timer."//NLC// self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if
        ! LCOV_EXCL_STOP

        self%nd%val = nd
        self%Decor = Decoration_type()    ! initialize the TAB character and decoration symbol to the default values.

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
        ! WARNING: The following ParaMonte library version tag as it will be replaced by the version 
        ! WARNING: that is generated by the preprocessing build script of the ParaMonte library.
        ! WARNING: This is superior to the above method of using the compiler Preprocessor
        ! WARNING: since CMAKE triggers a complete lengthy rebuild of the library when the
        ! WARNING: only the library version has changed.
        self%version  = "Unknown ParaMonte Version"
#include "ParaMonte_mod@version@kernel.inc.f90"
#endif

        ! setup general processor / coarray image variables

        call self%Image%query()

        ! setup formatting variables

        self%nd%str = num2str(self%nd%val)

        ! determine OS. Should be only needed by the Leader processes. But apparently not.

        call self%OS%query()
        if (self%OS%Err%occurred) then
        ! LCOV_EXCL_START
            self%Err = self%OS%Err
            self%Err%msg = PROCEDURE_NAME // ": Error occurred while querying OS type."//NLC//self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if
        ! LCOV_EXCL_STOP

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
            ! LCOV_EXCL_START
                self%Err = self%InputFile%Err
                self%Err%msg = PROCEDURE_NAME//": Error occurred while attempting to setup the user's input file='"//inputFile//"'."//NLC//self%Err%msg
                call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                return
            end if
            ! LCOV_EXCL_STOP
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

    !> \brief
    !> This procedure is a method of the [ParaMonte_type](@ref paramonte_type) class.
    !> Add a splash screen to the output report file.
    !>
    !> @param[inout]    self    :   An object of class [ParaMonte_type](@ref paramonte_type).
    !>
    !> \remark
    !> This routine has to be called by all master images (processes).
    subroutine addSplashScreen(self)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
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

    !> \brief
    !> Add information about the compiler and the platform/OS to the output report file(s).
    !>
    !> @param[inout]    self    :   An object of class [ParaMonte_type](@ref paramonte_type).
    !>
    !> \remark
    !> This procedure is a method of the [ParaMonte_type](@ref paramonte_type) class.
    !>
    !> \remark
    !> This routine has to be called by all leader images (processes).
    subroutine addCompilerPlatformInfo(self)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
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

        self%SystemInfo = SystemInfo_type(OS = self%OS, path = self%SpecBase%SystemInfoFilePath%val)
        if (self%SystemInfo%Err%occurred) then
        ! LCOV_EXCL_START
            self%Err = self%SystemInfo%Err
            self%Err%msg = PROCEDURE_NAME//": Error occurred while collecting system info."//NLC//self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if
        ! LCOV_EXCL_STOP

        ! write the system info to the output file

        do j = 1, self%SystemInfo%nRecord
            self%Decor%List = self%Decor%wrapText( self%SystemInfo%Records(j)%record , 132 )
            do i = 1,size(self%Decor%List)
                write(self%LogFile%unit,"(A)") self%Decor%List(i)%record
            end do
        end do
        call self%Decor%write(self%LogFile%unit)

    end subroutine addCompilerPlatformInfo

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a method of the [ParaMonte_type](@ref paramonte_type) class.
    !> Write to the output report file, the relevant platform setup messages.
    !>
    !> @param[inout]    self    :   An object of class [ParaMonte_type](@ref paramonte_type).
    !>
    !> \remark
    !> This routine has to be called by all master images (processes).
    subroutine noteUserAboutEnvSetup(self)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
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

    !> \brief
    !> If the relevant method name is missing in the namelist input file, then warn the user about this issue.
    !>
    !> @param[inout]    prefix      :   The prefix of the warning message.
    !> @param[inout]    name        :   The sampler method name.
    !> @param[inout]    namelist    :   The name of the missing namelist.
    !> @param[inout]    outputUnit  :   The file unit to which the message must be output.
    subroutine warnUserAboutMissingNamelist(prefix,name,namelist,outputUnit)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
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
        ! LCOV_EXCL_START
            call warn( prefix = prefix, outputUnit = output_unit, newline = "\n", msg = msg )
        end if
        ! LCOV_EXCL_STOP
    end subroutine warnUserAboutMissingNamelist

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a method of the [ParaMonte_type](@ref paramonte_type) class.
    !> Warn the user about whether the input file is missing, or is present, and other input file activities.
    !>
    !> @param[inout]    self    :   An object of class [ParaMonte_type](@ref paramonte_type).
    subroutine warnUserAboutInputFilePresence(self)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
        !DEC$ ATTRIBUTES DLLEXPORT :: warnUserAboutInputFilePresence
#endif
        use Constants_mod, only: NLC ! LCOV_EXCL_LINE
        implicit none
        class(ParaMonte_type), intent(inout)    :: self
        character(:), allocatable               :: msg
#if defined JULIA_ENABLED
        msg = "Interfacing Julia with "// self%name //"..."
#elif defined MATLAB_ENABLED
        msg = "Interfacing MATLAB with "// self%name //"..."
#elif defined MATTHEMATICA_ENABLED
        msg = "Interfacing Mathematica with "// self%name //"..."
#elif defined PYTHON_ENABLED
        msg = "Interfacing Python with "// self%name //"..."
#elif defined R_ENABLED
        msg = "Interfacing R with "// self%name //"..."
#elif defined C_ENABLED || defined CPP_ENABLED || defined FORTRAN_ENABLED
        if (self%inputFileArgIsPresent) then
            if (self%InputFile%exists) then
                msg =   "The user's input file for " // self%name // " options was detected."// NLC // &
                        "All " // self%name // " specifications will be read from the input file."// NLC // &
                        "Here is " // self%name // " input specifications file:"// NLC // self%inputFile%Path%modified
            elseif (self%InputFile%isInternal) then
                msg =   "No external file corresponding to the user's input file for "//self%name//" options could be found."//NLC// &
                        "The user-provided input file will be processed as an input string of "//self%name//" options."
            end if
        else
                msg =   "No " // self%name  // " input file is provided by the user."//NLC// &
#if defined FORTRAN_ENABLED
                        "Variable values from the procedure arguments will be used instead, where provided."//NLC//"Otherwise, "// &
#else
                        "Where needed, "// &
#endif
                        "the default options will be used."
        end if
#endif
        call self%note  ( prefix     = self%brand &
                        , outputUnit = self%LogFile%unit &
                        , newline    = NLC &
                        , msg        = msg )
    end subroutine warnUserAboutInputFilePresence

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !> \brief
    !> This procedure is a method of the [ParaMonte_type](@ref paramonte_type) class.
    !> Warn the user about whether the specifications setup from within the program are allowed or not.
    !>
    !> @param[inout]    self    :   An object of class [RefinedChain_type](@ref paramcmcrefinedchain_mod::refinedchain_type).
    subroutine setWarnAboutProcArgHasPriority(self)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
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

    !> \brief
    !> This procedure is a method of the [ParaMonte_type](@ref paramonte_type) class.
    !> Set up the output files of the simulation.
    !>
    !> @param[inout]    self    :   An object of class [RefinedChain_type](@ref paramcmcrefinedchain_mod::refinedchain_type).
    subroutine setupOutputFiles(self)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
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
        ! LCOV_EXCL_START
            self%Err = self%SpecBase%OutputFileName%Err
            self%Err%msg = PROCEDURE_NAME // ": Error occurred while attempting to construct OutputFileName path type." //NLC// self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if
        ! LCOV_EXCL_STOP

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
        ! LCOV_EXCL_START
            self%Err%msg = PROCEDURE_NAME//": Error occurred while fetching the current working directory via getcwd()."//NLC
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if
        ! LCOV_EXCL_STOP
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
            ! LCOV_EXCL_START
                self%Err%msg = PROCEDURE_NAME//": Error occurred while making directory = '"//self%SpecBase%OutputFileName%dir//"'."//NLC//self%Err%msg
                call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                return
            end if
            ! LCOV_EXCL_STOP
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
        ! LCOV_EXCL_START
            self%Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file='" // self%LogFile%Path%original // self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if
        ! LCOV_EXCL_STOP

        inquire( file = self%SampleFile%Path%original, exist = self%SampleFile%exists, iostat = self%SampleFile%Err%stat )
        self%Err = self%SampleFile%getInqErr( self%SampleFile%Err%stat )
        if (self%Err%occurred) then
        ! LCOV_EXCL_START
            self%Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file='" // self%SampleFile%Path%original // self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if
        ! LCOV_EXCL_STOP

        inquire( file = self%TimeFile%Path%original, exist = self%TimeFile%exists, iostat = self%TimeFile%Err%stat )
        self%Err = self%TimeFile%getInqErr( self%TimeFile%Err%stat )
        if (self%Err%occurred) then
        ! LCOV_EXCL_START
            self%Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file='" // self%TimeFile%Path%original // self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if
        ! LCOV_EXCL_STOP

        inquire( file = self%ChainFile%Path%original, exist = self%ChainFile%exists, iostat = self%ChainFile%Err%stat )
        self%Err = self%ChainFile%getInqErr( self%ChainFile%Err%stat )
        if (self%Err%occurred) then
        ! LCOV_EXCL_START
            self%Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file='" // self%ChainFile%Path%original // self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if
        ! LCOV_EXCL_STOP

        inquire( file = self%RestartFile%Path%original, exist = self%RestartFile%exists, iostat = self%RestartFile%Err%stat )
        self%Err = self%RestartFile%getInqErr( self%RestartFile%Err%stat )
        if (self%Err%occurred) then
        ! LCOV_EXCL_START
            self%Err%msg = PROCEDURE_NAME // ": Error occurred while inquiring the existence of file='" // self%RestartFile%Path%original // self%Err%msg
            call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
            return
        end if
        ! LCOV_EXCL_STOP

        self%isDryRun = (.not. self%SpecBase%OverwriteRequested%val) .and. & ! not fresh, if any file exists
                        (self%LogFile%exists .or. self%TimeFile%exists .or. self%RestartFile%exists .or. self%ChainFile%exists .or. self%SampleFile%exists)
        self%isFreshRun = .not. self%isDryRun

        if (self%isFreshRun) then
            if (self%Image%isFirst) call self%note( prefix = self%brand, outputUnit = self%LogFile%unit, newline = NLC, msg = "No pre-existing "//self%name//" run detected."//NLC//"Starting a fresh "//self%name//" run..." )
        else
            if (self%Image%isFirst) call self%note( prefix = self%brand, outputUnit = self%LogFile%unit, newline = NLC, msg = "Previous run of "//self%name//" detected."//NLC//"Searching for restart files..." )
            if (self%SampleFile%exists) then ! sampling is already complete
                self%Err%occurred = .true.
                self%Err%msg =  PROCEDURE_NAME//": Error occurred. Output sample file detected: "//self%SampleFile%Path%original//&
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

        if (self%Image%isLeader) then
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

        blockLeaderFileSetup: if (self%Image%isLeader) then

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
            ! LCOV_EXCL_START
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
            ! LCOV_EXCL_STOP

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
            ! LCOV_EXCL_START
                self%Err%msg = PROCEDURE_NAME // ": Error occurred while opening the " // self%name // " " // self%TimeFile%suffix // " file='" // self%TimeFile%Path%original // "'. "
                call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                return
            end if
            ! LCOV_EXCL_STOP

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
            ! LCOV_EXCL_START
                self%Err%msg = PROCEDURE_NAME // ": Error occurred while opening the " // self%name // " " // self%ChainFile%suffix // " file='" // self%ChainFile%Path%original // "'. "
                call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                return
            end if
            ! LCOV_EXCL_STOP

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
            ! LCOV_EXCL_START
                self%Err%msg = PROCEDURE_NAME // ": Error occurred while opening the " // self%name // " " // self%RestartFile%suffix // " file='" // self%RestartFile%Path%original // "'. "
                call self%abort( Err = self%Err, prefix = self%brand, newline = NLC, outputUnit = self%LogFile%unit )
                return
            end if
            ! LCOV_EXCL_STOP

            if (self%isFreshRun) then
                call self%Decor%writeDecoratedText  ( text = NLC // self%name // " simulation specifications" // NLC &
                                                    , marginTop = 1     &
                                                    , marginBot = 1     &
                                                    , newline = NLC     &
                                                    , outputUnit = self%LogFile%unit )
            end if

        end if blockLeaderFileSetup

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

    !> \brief
    !> This procedure is a method of the [ParaMonte_type](@ref paramonte_type) class.
    !> Output the relevant description.
    !>
    !> @param[inout]    self    :   An object of class [ParaMonte_type](@ref paramonte_type).
    !> @param[inout]    msg     :   The message to be output.
    subroutine reportDesc(self, msg) !, marginTop, marginBot)
#if IFORT_ENABLED && defined DLL_ENABLED && (OS_IS_WINDOWS || defined OS_IS_DARWIN) && !defined CFI_ENABLED
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

end module ParaMonte_mod ! LCOV_EXCL_LINE